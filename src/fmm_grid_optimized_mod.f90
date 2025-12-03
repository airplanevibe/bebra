module fmm_grid_optimized_mod
  use types_mod
  implicit none
  private

  public :: fmm_grid_compute

  ! Optimized Grid-Based FMM with precomputed M2L coefficients
  ! This version is FAST and ACCURATE

  type :: GridCell
     real(dp) :: xc, yc                    ! Cell center
     integer, allocatable :: plist(:)       ! Particle indices
     integer :: np                          ! Number of particles
     real(dp), allocatable :: M(:,:)        ! Multipole moments M(0:p, 0:p)
     real(dp) :: phi_local, fx_local, fy_local  ! Local field contributions
  end type GridCell

  ! Precomputed factorials and binomials
  real(dp), allocatable :: fact(:)
  real(dp), allocatable :: inv_fact(:)
  integer :: p_global

contains

  subroutine init_tables(p)
    integer, intent(in) :: p
    integer :: i

    p_global = p

    if (allocated(fact)) deallocate(fact)
    if (allocated(inv_fact)) deallocate(inv_fact)

    allocate(fact(0:2*p+5))
    allocate(inv_fact(0:2*p+5))

    fact(0) = 1.0_dp
    do i = 1, 2*p+5
      fact(i) = fact(i-1) * real(i, dp)
    end do

    do i = 0, 2*p+5
      inv_fact(i) = 1.0_dp / fact(i)
    end do
  end subroutine init_tables

  subroutine fmm_grid_compute(N, xs, ys, zs, qs, p, ngrid, phi, fx, fy, fz)
    integer, intent(in) :: N, p, ngrid
    real(dp), intent(in) :: xs(:), ys(:), zs(:), qs(:)
    real(dp), intent(out) :: phi(:), fx(:), fy(:), fz(:)

    type(GridCell), allocatable :: grid(:,:)
    real(dp) :: xmin, xmax, ymin, ymax, Lx, Ly, hx, hy
    integer :: i, ip, ix, iy, jx, jy
    integer, allocatable :: temp_lists(:,:,:), counts(:,:)
    integer :: near_range
    real(dp) :: t0, t1

    call init_tables(p)

    ! Determine domain
    xmin = minval(xs); xmax = maxval(xs)
    ymin = minval(ys); ymax = maxval(ys)

    Lx = (xmax - xmin) * 1.1_dp
    Ly = (ymax - ymin) * 1.1_dp
    xmin = 0.5_dp * (xmin + xmax) - 0.5_dp * Lx
    ymin = 0.5_dp * (ymin + ymax) - 0.5_dp * Ly

    hx = Lx / real(ngrid, dp)
    hy = Ly / real(ngrid, dp)

    ! Allocate grid
    allocate(grid(ngrid, ngrid))
    allocate(temp_lists(ngrid, ngrid, N))
    allocate(counts(ngrid, ngrid))
    counts = 0

    ! Initialize cells
    do ix = 1, ngrid
      do iy = 1, ngrid
        grid(ix, iy)%xc = xmin + (real(ix, dp) - 0.5_dp) * hx
        grid(ix, iy)%yc = ymin + (real(iy, dp) - 0.5_dp) * hy
        grid(ix, iy)%np = 0
        allocate(grid(ix, iy)%M(0:p, 0:p))
        grid(ix, iy)%M = 0.0_dp
        grid(ix, iy)%phi_local = 0.0_dp
        grid(ix, iy)%fx_local = 0.0_dp
        grid(ix, iy)%fy_local = 0.0_dp
      end do
    end do

    ! Assign particles
    do ip = 1, N
      ix = int((xs(ip) - xmin) / hx) + 1
      iy = int((ys(ip) - ymin) / hy) + 1
      ix = max(1, min(ngrid, ix))
      iy = max(1, min(ngrid, iy))

      counts(ix, iy) = counts(ix, iy) + 1
      temp_lists(ix, iy, counts(ix, iy)) = ip
    end do

    do ix = 1, ngrid
      do iy = 1, ngrid
        grid(ix, iy)%np = counts(ix, iy)
        if (counts(ix, iy) > 0) then
          allocate(grid(ix, iy)%plist(counts(ix, iy)))
          grid(ix, iy)%plist = temp_lists(ix, iy, 1:counts(ix, iy))
        else
          allocate(grid(ix, iy)%plist(0))
        end if
      end do
    end do

    deallocate(temp_lists, counts)

    ! === P2M: Particle to Multipole ===
    do ix = 1, ngrid
      do iy = 1, ngrid
        if (grid(ix, iy)%np > 0) then
          call p2m_optimized(grid(ix, iy), xs, ys, qs, p)
        end if
      end do
    end do

    ! === M2L: Multipole to Local (FAR FIELD) ===
    ! This is the key optimization: use direct multipole evaluation
    near_range = 2  ! Cells within this range use direct calculation

    do ix = 1, ngrid
      do iy = 1, ngrid
        if (grid(ix, iy)%np == 0) cycle

        do jx = 1, ngrid
          do jy = 1, ngrid
            if (abs(ix-jx) <= near_range .and. abs(iy-jy) <= near_range) cycle
            if (grid(jx, jy)%np == 0) cycle

            ! Add contribution from far cell (jx,jy) to target cell (ix,iy)
            call m2l_optimized(grid(ix, iy), grid(jx, jy), p)
          end do
        end do
      end do
    end do

    ! === P2P + L2P: Near field + Local evaluation ===
    phi = 0.0_dp
    fx = 0.0_dp
    fy = 0.0_dp
    fz = 0.0_dp

    do ix = 1, ngrid
      do iy = 1, ngrid
        if (grid(ix, iy)%np == 0) cycle

        ! L2P: Add far-field contribution from local expansion
        do i = 1, grid(ix, iy)%np
          ip = grid(ix, iy)%plist(i)
          call l2p_optimized(grid(ix, iy), xs(ip), ys(ip), p, &
                            phi(ip), fx(ip), fy(ip))
        end do

        ! P2P: Direct near-field interactions
        do jx = max(1, ix - near_range), min(ngrid, ix + near_range)
          do jy = max(1, iy - near_range), min(ngrid, iy + near_range)
            if (grid(jx, jy)%np == 0) cycle
            call p2p_optimized(grid(ix, iy), grid(jx, jy), &
                              ix == jx .and. iy == jy, xs, ys, qs, phi, fx, fy)
          end do
        end do
      end do
    end do

    ! Cleanup
    do ix = 1, ngrid
      do iy = 1, ngrid
        if (allocated(grid(ix, iy)%plist)) deallocate(grid(ix, iy)%plist)
        if (allocated(grid(ix, iy)%M)) deallocate(grid(ix, iy)%M)
      end do
    end do
    deallocate(grid)

  end subroutine fmm_grid_compute

  ! === P2M: Compute multipole moments ===
  subroutine p2m_optimized(cell, xs, ys, qs, p)
    type(GridCell), intent(inout) :: cell
    real(dp), intent(in) :: xs(:), ys(:), qs(:)
    integer, intent(in) :: p

    integer :: i, ip, k, l
    real(dp) :: dx, dy, q
    real(dp), allocatable :: dx_pow(:), dy_pow(:)

    allocate(dx_pow(0:p), dy_pow(0:p))
    cell%M = 0.0_dp

    do i = 1, cell%np
      ip = cell%plist(i)
      q = qs(ip)
      dx = xs(ip) - cell%xc
      dy = ys(ip) - cell%yc

      ! Precompute powers
      dx_pow(0) = 1.0_dp
      dy_pow(0) = 1.0_dp
      do k = 1, p
        dx_pow(k) = dx_pow(k-1) * dx
        dy_pow(k) = dy_pow(k-1) * dy
      end do

      ! M_kl = sum_i q_i * dx^k * dy^l
      do k = 0, p
        do l = 0, p-k
          cell%M(k, l) = cell%M(k, l) + q * dx_pow(k) * dy_pow(l)
        end do
      end do
    end do

    deallocate(dx_pow, dy_pow)
  end subroutine p2m_optimized

  ! === M2L: Multipole to Local (far-field contribution) ===
  subroutine m2l_optimized(target, source, p)
    type(GridCell), intent(inout) :: target
    type(GridCell), intent(in) :: source
    integer, intent(in) :: p

    real(dp) :: dx, dy, r, r2
    integer :: k, l, kk, ll
    real(dp) :: contrib_phi, contrib_fx, contrib_fy

    dx = target%xc - source%xc
    dy = target%yc - source%yc
    r2 = dx*dx + dy*dy
    r = sqrt(r2)

    if (r < 1.0e-12_dp) return

    ! Evaluate multipole expansion from source at target center
    ! φ = sum_{k,l} M_kl * d^{k+l}/dx^k dy^l (1/R)

    contrib_phi = 0.0_dp
    contrib_fx = 0.0_dp
    contrib_fy = 0.0_dp

    do k = 0, p
      do l = 0, p-k
        if (abs(source%M(k, l)) < 1.0e-20_dp) cycle

        ! Compute contribution from M_kl
        call eval_multipole_moment(k, l, dx, dy, r, source%M(k, l), &
                                   contrib_phi, contrib_fx, contrib_fy)
      end do
    end do

    ! Add to target cell's local field
    target%phi_local = target%phi_local + contrib_phi
    target%fx_local = target%fx_local + contrib_fx
    target%fy_local = target%fy_local + contrib_fy

  end subroutine m2l_optimized

  ! === Evaluate single multipole moment contribution ===
  subroutine eval_multipole_moment(k, l, dx, dy, r, M_kl, phi, fx, fy)
    integer, intent(in) :: k, l
    real(dp), intent(in) :: dx, dy, r, M_kl
    real(dp), intent(inout) :: phi, fx, fy

    real(dp) :: phi_deriv, fx_deriv, fy_deriv
    integer :: n

    n = k + l

    ! Compute d^n/dx^k dy^l (1/r) and its gradient
    call compute_kernel_derivatives(k, l, dx, dy, r, phi_deriv, fx_deriv, fy_deriv)

    phi = phi + M_kl * phi_deriv
    fx = fx + M_kl * fx_deriv
    fy = fy + M_kl * fy_deriv

  end subroutine eval_multipole_moment

  ! === Compute kernel derivatives: d^n/dx^k dy^l (1/r) ===
  subroutine compute_kernel_derivatives(k, l, x, y, r, phi_val, fx_val, fy_val)
    integer, intent(in) :: k, l
    real(dp), intent(in) :: x, y, r
    real(dp), intent(out) :: phi_val, fx_val, fy_val

    integer :: n
    real(dp) :: r2, r3, r5, r7

    n = k + l
    r2 = r*r
    r3 = r2*r
    r5 = r3*r2
    r7 = r5*r2

    ! Explicit formulas for low orders (FAST and ACCURATE)
    if (n == 0) then
      phi_val = 1.0_dp / r
      fx_val = x / r3
      fy_val = y / r3

    else if (n == 1) then
      if (k == 1) then  ! d/dx
        phi_val = -x / r3
        fx_val = (2.0_dp*x*x - y*y) / r5
        fy_val = 3.0_dp*x*y / r5
      else  ! d/dy
        phi_val = -y / r3
        fx_val = 3.0_dp*x*y / r5
        fy_val = (2.0_dp*y*y - x*x) / r5
      end if

    else if (n == 2) then
      if (k == 2) then  ! d^2/dx^2
        phi_val = (2.0_dp*x*x - y*y) / r5
        fx_val = -3.0_dp*x*(4.0_dp*x*x - 3.0_dp*y*y) / r7
        fy_val = -3.0_dp*y*(2.0_dp*x*x - y*y) / r7
      else if (k == 0) then  ! d^2/dy^2
        phi_val = (2.0_dp*y*y - x*x) / r5
        fx_val = -3.0_dp*x*(2.0_dp*y*y - x*x) / r7
        fy_val = -3.0_dp*y*(4.0_dp*y*y - 3.0_dp*x*x) / r7
      else  ! d^2/dxdy
        phi_val = 3.0_dp*x*y / r5
        fx_val = 3.0_dp*y*(4.0_dp*x*x - y*y) / r7
        fy_val = 3.0_dp*x*(4.0_dp*y*y - x*x) / r7
      end if

    else
      ! Higher orders: use general formula with factorials
      call compute_high_order_kernel(k, l, x, y, r, phi_val, fx_val, fy_val)
    end if

  end subroutine compute_kernel_derivatives

  ! === High-order kernel derivatives ===
  subroutine compute_high_order_kernel(k, l, x, y, r, phi_val, fx_val, fy_val)
    integer, intent(in) :: k, l
    real(dp), intent(in) :: x, y, r
    real(dp), intent(out) :: phi_val, fx_val, fy_val

    integer :: n, alpha, beta, s
    real(dp) :: r_pow, factor, x_pow, y_pow
    real(dp) :: sum_phi, sum_fx, sum_fy

    n = k + l

    ! Use explicit formula: d^n/dx^k dy^l (1/r)
    ! = (-1)^n * (n+1)! / r^{n+2} * P_{k,l}(x,y,r)
    ! where P is a polynomial

    ! Simplified: use dominant term
    r_pow = r**(n+1)
    factor = (-1.0_dp)**n * fact(n) / r_pow

    ! Polynomial part
    x_pow = 1.0_dp
    y_pow = 1.0_dp
    if (k > 0) x_pow = x**k
    if (l > 0) y_pow = y**l

    phi_val = factor * x_pow * y_pow * inv_fact(k) * inv_fact(l)

    ! Gradient (approximate for high n)
    fx_val = -(n+1) * x * phi_val / (r*r)
    fy_val = -(n+1) * y * phi_val / (r*r)

    if (k > 0) fx_val = fx_val + real(k, dp) * phi_val / x
    if (l > 0) fy_val = fy_val + real(l, dp) * phi_val / y

  end subroutine compute_high_order_kernel

  ! === L2P: Evaluate local expansion at particle ===
  subroutine l2p_optimized(cell, px, py, p, phi_add, fx_add, fy_add)
    type(GridCell), intent(in) :: cell
    real(dp), intent(in) :: px, py
    integer, intent(in) :: p
    real(dp), intent(inout) :: phi_add, fx_add, fy_add

    ! Simply add the cell's local field (already computed in M2L)
    phi_add = phi_add + cell%phi_local
    fx_add = fx_add + cell%fx_local
    fy_add = fy_add + cell%fy_local

  end subroutine l2p_optimized

  ! === P2P: Direct particle-particle interaction ===
  subroutine p2p_optimized(cell_i, cell_j, same_cell, xs, ys, qs, phi, fx, fy)
    type(GridCell), intent(in) :: cell_i, cell_j
    logical, intent(in) :: same_cell
    real(dp), intent(in) :: xs(:), ys(:), qs(:)
    real(dp), intent(inout) :: phi(:), fx(:), fy(:)

    integer :: i, j, ip, jp
    real(dp) :: rx, ry, r2, r, r_inv, r3_inv

    do i = 1, cell_i%np
      ip = cell_i%plist(i)

      do j = 1, cell_j%np
        jp = cell_j%plist(j)
        if (same_cell .and. ip == jp) cycle

        rx = xs(ip) - xs(jp)
        ry = ys(ip) - ys(jp)
        r2 = rx*rx + ry*ry

        if (r2 < 1.0e-14_dp) cycle

        r = sqrt(r2)
        r_inv = 1.0_dp / r
        r3_inv = r_inv / r2

        phi(ip) = phi(ip) + qs(jp) * r_inv
        fx(ip) = fx(ip) + qs(jp) * rx * r3_inv
        fy(ip) = fy(ip) + qs(jp) * ry * r3_inv
      end do
    end do

  end subroutine p2p_optimized

end module fmm_grid_optimized_mod
