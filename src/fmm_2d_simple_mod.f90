module fmm_2d_simple_mod
  use types_mod
  implicit none
  private

  public :: fmm_2d_compute

  ! Simple 2D Cartesian multipole expansion for Coulomb potential
  ! Using monopole, dipole, quadrupole, octupole, etc.
  ! Much simpler and faster than full spherical harmonic approach

  type :: GridCell2D
     real(dp) :: xc, yc                   ! Cell center
     integer, allocatable :: plist(:)      ! Particle indices
     integer :: np                         ! Number of particles
     ! Multipole moments up to order p
     real(dp), allocatable :: moments(:,:) ! (0:p, 0:p) - x^k * y^l moments
  end type GridCell2D

  ! Precomputed factorials
  real(dp), allocatable :: fact(:)
  integer :: p_global

contains

  subroutine init_factorials(pmax)
    integer, intent(in) :: pmax
    integer :: i

    if (allocated(fact)) deallocate(fact)
    allocate(fact(0:2*pmax+10))

    fact(0) = 1.0_dp
    do i = 1, 2*pmax+10
      fact(i) = fact(i-1) * real(i, dp)
    end do
    p_global = pmax
  end subroutine init_factorials

  subroutine fmm_2d_compute(N, xs, ys, zs, qs, p, ngrid, phi, fx, fy, fz)
    integer, intent(in) :: N, p, ngrid
    real(dp), intent(in) :: xs(:), ys(:), zs(:), qs(:)
    real(dp), intent(out) :: phi(:), fx(:), fy(:), fz(:)

    type(GridCell2D), allocatable :: grid(:,:)
    real(dp) :: xmin, xmax, ymin, ymax, Lx, Ly, hx, hy
    integer :: i, ip, ix, iy, jx, jy
    integer, allocatable :: temp_lists(:,:,:)
    integer, allocatable :: counts(:,:)
    integer :: near_range

    call init_factorials(p)

    ! Domain bounds
    xmin = minval(xs)
    xmax = maxval(xs)
    ymin = minval(ys)
    ymax = maxval(ys)

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

    do ix = 1, ngrid
      do iy = 1, ngrid
        grid(ix, iy)%xc = xmin + (real(ix, dp) - 0.5_dp) * hx
        grid(ix, iy)%yc = ymin + (real(iy, dp) - 0.5_dp) * hy
        grid(ix, iy)%np = 0
        allocate(grid(ix, iy)%moments(0:p, 0:p))
        grid(ix, iy)%moments = 0.0_dp
      end do
    end do

    ! Assign particles to cells
    do ip = 1, N
      ix = int((xs(ip) - xmin) / hx) + 1
      iy = int((ys(ip) - ymin) / hy) + 1

      ix = max(1, min(ngrid, ix))
      iy = max(1, min(ngrid, iy))

      counts(ix, iy) = counts(ix, iy) + 1
      temp_lists(ix, iy, counts(ix, iy)) = ip
    end do

    ! Copy particle lists
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

    ! P2M: Compute multipole moments
    !$OMP PARALLEL DO PRIVATE(ix, iy) SCHEDULE(dynamic)
    do ix = 1, ngrid
      do iy = 1, ngrid
        if (grid(ix, iy)%np > 0) then
          call compute_multipole_moments(grid(ix, iy), xs, ys, qs, p)
        end if
      end do
    end do
    !$OMP END PARALLEL DO

    ! Evaluate at particles
    phi = 0.0_dp
    fx = 0.0_dp
    fy = 0.0_dp
    fz = 0.0_dp

    near_range = max(2, int(sqrt(real(p, dp)) + 1))

    !$OMP PARALLEL DO PRIVATE(ix, iy, i, ip, jx, jy) SCHEDULE(dynamic)
    do ix = 1, ngrid
      do iy = 1, ngrid
        if (grid(ix, iy)%np == 0) cycle

        ! Far-field: multipole expansion from distant cells
        do i = 1, grid(ix, iy)%np
          ip = grid(ix, iy)%plist(i)

          do jx = 1, ngrid
            do jy = 1, ngrid
              if (ix == jx .and. iy == jy) cycle
              if (abs(ix - jx) <= near_range .and. abs(iy - jy) <= near_range) cycle
              if (grid(jx, jy)%np == 0) cycle

              call evaluate_multipole_at_point(xs(ip), ys(ip), grid(jx, jy), p, &
                                              phi(ip), fx(ip), fy(ip))
            end do
          end do
        end do

        ! Near-field: direct calculation
        do jx = max(1, ix - near_range), min(ngrid, ix + near_range)
          do jy = max(1, iy - near_range), min(ngrid, iy + near_range)
            if (grid(jx, jy)%np == 0) cycle

            call direct_interaction(grid(ix, iy), grid(jx, jy), &
                                   ix == jx .and. iy == jy, &
                                   xs, ys, qs, phi, fx, fy)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ! Cleanup
    do ix = 1, ngrid
      do iy = 1, ngrid
        if (allocated(grid(ix, iy)%plist)) deallocate(grid(ix, iy)%plist)
        if (allocated(grid(ix, iy)%moments)) deallocate(grid(ix, iy)%moments)
      end do
    end do
    deallocate(grid)
    if (allocated(fact)) deallocate(fact)

  end subroutine fmm_2d_compute

  ! Compute multipole moments up to order p
  subroutine compute_multipole_moments(cell, xs, ys, qs, p)
    type(GridCell2D), intent(inout) :: cell
    real(dp), intent(in) :: xs(:), ys(:), qs(:)
    integer, intent(in) :: p

    integer :: i, ip, k, l
    real(dp) :: dx, dy, q

    cell%moments = 0.0_dp

    do i = 1, cell%np
      ip = cell%plist(i)
      q = qs(ip)
      dx = xs(ip) - cell%xc
      dy = ys(ip) - cell%yc

      ! M_kl = sum_i q_i * dx^k * dy^l
      do k = 0, p
        do l = 0, p - k  ! Limit total order
          cell%moments(k, l) = cell%moments(k, l) + q * (dx**k) * (dy**l)
        end do
      end do
    end do

  end subroutine compute_multipole_moments

  ! Evaluate multipole expansion at a point using exact formulas
  subroutine evaluate_multipole_at_point(px, py, source, p, phi_add, fx_add, fy_add)
    real(dp), intent(in) :: px, py
    type(GridCell2D), intent(in) :: source
    integer, intent(in) :: p
    real(dp), intent(inout) :: phi_add, fx_add, fy_add

    real(dp) :: dx, dy, r, r2
    integer :: k, l, n, alpha, beta
    real(dp) :: coeff, r_pow
    real(dp) :: phi_contrib, fx_contrib, fy_contrib

    dx = px - source%xc
    dy = py - source%yc
    r2 = dx*dx + dy*dy
    r = sqrt(r2)

    if (r < 1.0e-12_dp) return

    ! Use Taylor expansion of 1/|r-r'|
    ! phi(r) = sum_{k,l} M_kl * d^{k+l}/dx^k dy^l (1/R)|_{R=(dx,dy)}

    do k = 0, p
      do l = 0, p - k
        if (abs(source%moments(k, l)) < 1.0e-20_dp) cycle

        n = k + l
        call evaluate_derivative_1_over_r(k, l, dx, dy, r, &
                                          phi_contrib, fx_contrib, fy_contrib)

        phi_add = phi_add + source%moments(k, l) * phi_contrib
        fx_add = fx_add + source%moments(k, l) * fx_contrib
        fy_add = fy_add + source%moments(k, l) * fy_contrib
      end do
    end do

  end subroutine evaluate_multipole_at_point

  ! Compute d^{k+l}/dx^k dy^l (1/r) and its gradient
  ! Using explicit formulas for 2D
  subroutine evaluate_derivative_1_over_r(k, l, dx, dy, r, phi_val, fx_val, fy_val)
    integer, intent(in) :: k, l
    real(dp), intent(in) :: dx, dy, r
    real(dp), intent(out) :: phi_val, fx_val, fy_val

    integer :: n, i, j, s
    real(dp) :: r2, factor, sign_factor
    real(dp) :: phi_temp, dx_pow, dy_pow

    n = k + l
    r2 = dx*dx + dy*dy

    if (n == 0) then
      ! 1/r
      phi_val = 1.0_dp / r
      fx_val = dx / (r*r2)
      fy_val = dy / (r*r2)

    else if (n == 1) then
      if (k == 1) then
        ! d/dx (1/r) = -x/r^3
        phi_val = -dx / (r*r2)
        fx_val = (2.0_dp*dx*dx - dy*dy) / (r2*r*r2)
        fy_val = 3.0_dp*dx*dy / (r2*r*r2)
      else
        ! d/dy (1/r) = -y/r^3
        phi_val = -dy / (r*r2)
        fx_val = 3.0_dp*dx*dy / (r2*r*r2)
        fy_val = (2.0_dp*dy*dy - dx*dx) / (r2*r*r2)
      end if

    else if (n == 2) then
      if (k == 2) then
        ! d^2/dx^2 (1/r)
        phi_val = (2.0_dp*dx*dx - dy*dy) / (r2*r*r2)
        fx_val = -6.0_dp*dx*(2.0_dp*dx*dx - dy*dy) / (r2*r2*r*r2) + 4.0_dp*dx/(r2*r*r2)
        fy_val = -6.0_dp*dy*(2.0_dp*dx*dx - dy*dy) / (r2*r2*r*r2) - 2.0_dp*dy/(r2*r*r2)
      else if (k == 0) then
        ! d^2/dy^2 (1/r)
        phi_val = (2.0_dp*dy*dy - dx*dx) / (r2*r*r2)
        fx_val = -6.0_dp*dx*(2.0_dp*dy*dy - dx*dx) / (r2*r2*r*r2) - 2.0_dp*dx/(r2*r*r2)
        fy_val = -6.0_dp*dy*(2.0_dp*dy*dy - dx*dx) / (r2*r2*r*r2) + 4.0_dp*dy/(r2*r*r2)
      else ! k == 1, l == 1
        ! d^2/dxdy (1/r) = 3xy/r^5
        phi_val = 3.0_dp*dx*dy / (r2*r*r2)
        fx_val = 3.0_dp*dy*(2.0_dp*dy*dy - 3.0_dp*dx*dx) / (r2*r2*r*r2)
        fy_val = 3.0_dp*dx*(2.0_dp*dx*dx - 3.0_dp*dy*dy) / (r2*r2*r*r2)
      end if

    else
      ! Higher orders: use general formula
      ! d^n/dx^k dy^l (1/r) = sum of terms involving powers of x, y, 1/r
      call compute_high_order_derivative(k, l, dx, dy, r, phi_val, fx_val, fy_val)
    end if

  end subroutine evaluate_derivative_1_over_r

  ! General high-order derivative computation
  subroutine compute_high_order_derivative(k, l, dx, dy, r, phi_val, fx_val, fy_val)
    integer, intent(in) :: k, l
    real(dp), intent(in) :: dx, dy, r
    real(dp), intent(out) :: phi_val, fx_val, fy_val

    integer :: n, m, alpha, beta
    real(dp) :: r2, factor
    real(dp) :: sum_phi, sum_fx, sum_fy
    real(dp) :: term

    n = k + l
    r2 = dx*dx + dy*dy

    ! Use finite difference approximation for simplicity
    ! Or implement recursive formula

    ! For now, use approximate formula
    factor = (-1.0_dp)**n * fact(n) / (r**(n+1))

    ! Approximate using dominant term
    if (k > 0 .and. l > 0) then
      phi_val = factor * (dx**k) * (dy**l) / (fact(k) * fact(l))
    else if (k > 0) then
      phi_val = factor * (dx**k) / fact(k)
    else if (l > 0) then
      phi_val = factor * (dy**l) / fact(l)
    else
      phi_val = 1.0_dp / r
    end if

    ! Gradients (approximate)
    fx_val = -(n+1) * dx * phi_val / r2
    fy_val = -(n+1) * dy * phi_val / r2

    if (k > 0) fx_val = fx_val + real(k, dp) * phi_val / dx
    if (l > 0) fy_val = fy_val + real(l, dp) * phi_val / dy

  end subroutine compute_high_order_derivative

  ! Direct pairwise interaction
  subroutine direct_interaction(cell_i, cell_j, same_cell, xs, ys, qs, phi, fx, fy)
    type(GridCell2D), intent(in) :: cell_i, cell_j
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

  end subroutine direct_interaction

end module fmm_2d_simple_mod
