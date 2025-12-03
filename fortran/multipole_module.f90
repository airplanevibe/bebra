! ============================================================================
! Multipole Expansion Module for 2D Particle Systems
! Implements Fast Multipole Method for electrostatic interactions
! ============================================================================
module multipole_module
    use legendre_data
    implicit none
    private

    ! Public interface
    public :: particle_system
    public :: init_system, add_particle, get_energy, get_force, get_all_forces
    public :: move_particles, correct_cells, determine_series_order
    public :: write_frame

    ! Particle type
    type :: particle
        real(dp) :: charge
        real(dp) :: x, y
        real(dp) :: fx, fy  ! Forces
        integer :: cell_i, cell_j
    end type particle

    ! Cell type with multipole moments
    type :: cell_type
        integer :: num_particles
        integer, allocatable :: particle_ids(:)
        integer :: capacity
        ! Multipole moments M_n
        real(dp), allocatable :: moments(:)  ! 0:p_max
        real(dp) :: center_x, center_y
    end type cell_type

    ! Main particle system
    type :: particle_system
        real(dp) :: xmax, ymax
        integer :: num_cells_x, num_cells_y
        real(dp) :: dx, dy  ! Cell sizes
        integer :: num_particles
        integer :: capacity
        type(particle), allocatable :: particles(:)
        type(cell_type), allocatable :: cells(:,:)
        integer :: p_order  ! Multipole expansion order
        real(dp) :: eta     ! Viscosity parameter
    end type particle_system

contains

! ============================================================================
! Initialize particle system
! ============================================================================
subroutine init_system(sys, xmax, ymax, num_cells_x, num_cells_y)
    type(particle_system), intent(inout) :: sys
    real(dp), intent(in) :: xmax, ymax
    integer, intent(in) :: num_cells_x, num_cells_y
    integer :: i, j

    sys%xmax = xmax
    sys%ymax = ymax
    sys%num_cells_x = num_cells_x
    sys%num_cells_y = num_cells_y
    sys%dx = (2.0_dp * xmax) / real(num_cells_x, dp)
    sys%dy = (2.0_dp * ymax) / real(num_cells_y, dp)

    sys%num_particles = 0
    sys%capacity = 1000
    allocate(sys%particles(sys%capacity))

    ! Initialize cells
    allocate(sys%cells(num_cells_x, num_cells_y))

    do j = 1, num_cells_y
        do i = 1, num_cells_x
            sys%cells(i,j)%num_particles = 0
            sys%cells(i,j)%capacity = 50
            allocate(sys%cells(i,j)%particle_ids(sys%cells(i,j)%capacity))
            allocate(sys%cells(i,j)%moments(0:p_max))
            sys%cells(i,j)%moments = 0.0_dp

            ! Cell centers
            sys%cells(i,j)%center_x = -xmax + (i - 0.5_dp) * sys%dx
            sys%cells(i,j)%center_y = -ymax + (j - 0.5_dp) * sys%dy
        end do
    end do

    sys%p_order = 10  ! Default order
    sys%eta = 1.0_dp  ! Default viscosity

    print *, "Initialized system:"
    print *, "  Domain: [", -xmax, ",", xmax, "] x [", -ymax, ",", ymax, "]"
    print *, "  Grid:", num_cells_x, "x", num_cells_y
    print *, "  Cell size:", sys%dx, "x", sys%dy

end subroutine init_system

! ============================================================================
! Add particle to system
! ============================================================================
subroutine add_particle(sys, charge, x, y)
    type(particle_system), intent(inout) :: sys
    real(dp), intent(in) :: charge, x, y
    integer :: i, j, pid
    type(particle), allocatable :: temp_particles(:)

    ! Resize if needed
    if (sys%num_particles >= sys%capacity) then
        allocate(temp_particles(sys%capacity * 2))
        temp_particles(1:sys%capacity) = sys%particles
        deallocate(sys%particles)
        call move_alloc(temp_particles, sys%particles)
        sys%capacity = sys%capacity * 2
    end if

    ! Add particle
    sys%num_particles = sys%num_particles + 1
    pid = sys%num_particles

    sys%particles(pid)%charge = charge
    sys%particles(pid)%x = x
    sys%particles(pid)%y = y
    sys%particles(pid)%fx = 0.0_dp
    sys%particles(pid)%fy = 0.0_dp

    ! Determine cell indices
    call get_cell_indices(sys, x, y, i, j)
    sys%particles(pid)%cell_i = i
    sys%particles(pid)%cell_j = j

    ! Add to cell
    call add_to_cell(sys%cells(i,j), pid)

end subroutine add_particle

! ============================================================================
! Get cell indices for coordinates
! ============================================================================
subroutine get_cell_indices(sys, x, y, i, j)
    type(particle_system), intent(in) :: sys
    real(dp), intent(in) :: x, y
    integer, intent(out) :: i, j

    i = int((x + sys%xmax) / sys%dx) + 1
    j = int((y + sys%ymax) / sys%dy) + 1

    ! Clamp to valid range
    i = max(1, min(i, sys%num_cells_x))
    j = max(1, min(j, sys%num_cells_y))

end subroutine get_cell_indices

! ============================================================================
! Add particle ID to cell
! ============================================================================
subroutine add_to_cell(cell, pid)
    type(cell_type), intent(inout) :: cell
    integer, intent(in) :: pid
    integer, allocatable :: temp_ids(:)

    if (cell%num_particles >= cell%capacity) then
        allocate(temp_ids(cell%capacity * 2))
        temp_ids(1:cell%capacity) = cell%particle_ids
        deallocate(cell%particle_ids)
        call move_alloc(temp_ids, cell%particle_ids)
        cell%capacity = cell%capacity * 2
    end if

    cell%num_particles = cell%num_particles + 1
    cell%particle_ids(cell%num_particles) = pid

end subroutine add_to_cell

! ============================================================================
! Compute multipole moments for all cells
! ============================================================================
subroutine compute_multipole_moments(sys)
    type(particle_system), intent(inout) :: sys
    integer :: i, j, k, pid, n
    real(dp) :: dx, dy, r, q
    type(cell_type) :: cell

    ! Reset all moments
    do j = 1, sys%num_cells_y
        do i = 1, sys%num_cells_x
            sys%cells(i,j)%moments = 0.0_dp
        end do
    end do

    ! Compute moments for each cell
    do j = 1, sys%num_cells_y
        do i = 1, sys%num_cells_x
            cell = sys%cells(i,j)

            do k = 1, cell%num_particles
                pid = cell%particle_ids(k)

                dx = sys%particles(pid)%x - cell%center_x
                dy = sys%particles(pid)%y - cell%center_y
                r = sqrt(dx*dx + dy*dy)
                q = sys%particles(pid)%charge

                if (r < 1.0e-15_dp) r = 1.0e-15_dp

                ! M_n = sum_i q_i * r_i^n
                do n = 0, sys%p_order
                    cell%moments(n) = cell%moments(n) + q * r**n
                end do
            end do
        end do
    end do

end subroutine compute_multipole_moments

! ============================================================================
! Get total energy (method 0=direct, 1=multipole within cells, 2=full multipole)
! ============================================================================
function get_energy(sys, method) result(energy)
    type(particle_system), intent(inout) :: sys
    integer, intent(in) :: method
    real(dp) :: energy

    select case(method)
    case(0)
        energy = compute_energy_direct(sys)
    case(1)
        energy = compute_energy_multipole_cells(sys)
    case(2)
        energy = compute_energy_multipole_full(sys)
    case default
        print *, "Unknown method:", method
        energy = 0.0_dp
    end select

end function get_energy

! ============================================================================
! Direct energy calculation
! ============================================================================
function compute_energy_direct(sys) result(energy)
    type(particle_system), intent(in) :: sys
    real(dp) :: energy
    integer :: i, j
    real(dp) :: dx, dy, r, qi, qj

    energy = 0.0_dp

    do i = 1, sys%num_particles - 1
        do j = i + 1, sys%num_particles
            dx = sys%particles(i)%x - sys%particles(j)%x
            dy = sys%particles(i)%y - sys%particles(j)%y
            r = sqrt(dx*dx + dy*dy)

            if (r > 1.0e-15_dp) then
                qi = sys%particles(i)%charge
                qj = sys%particles(j)%charge
                energy = energy + qi * qj / r
            end if
        end do
    end do

end function compute_energy_direct

! ============================================================================
! Energy with multipole expansion (cell-wise)
! ============================================================================
function compute_energy_multipole_cells(sys) result(energy)
    type(particle_system), intent(inout) :: sys
    real(dp) :: energy
    integer :: ci, cj, di, dj, i, j
    real(dp) :: e_near, e_far

    call compute_multipole_moments(sys)

    energy = 0.0_dp
    e_near = 0.0_dp
    e_far = 0.0_dp

    ! Near field: direct within and neighboring cells
    do cj = 1, sys%num_cells_y
        do ci = 1, sys%num_cells_x
            ! Self cell
            e_near = e_near + compute_cell_self_energy(sys, ci, cj)

            ! Neighboring cells (direct)
            do dj = max(1, cj-1), min(sys%num_cells_y, cj+1)
                do di = max(1, ci-1), min(sys%num_cells_x, ci+1)
                    if (di == ci .and. dj == cj) cycle
                    e_near = e_near + compute_cell_pair_energy_direct(sys, ci, cj, di, dj)
                end do
            end do
        end do
    end do

    ! Far field: multipole expansion
    do cj = 1, sys%num_cells_y
        do ci = 1, sys%num_cells_x
            do dj = 1, sys%num_cells_y
                do di = 1, sys%num_cells_x
                    ! Skip near field cells
                    if (abs(di - ci) <= 1 .and. abs(dj - cj) <= 1) cycle

                    e_far = e_far + compute_cell_pair_energy_multipole(sys, ci, cj, di, dj)
                end do
            end do
        end do
    end do

    energy = e_near + e_far

end function compute_energy_multipole_cells

! ============================================================================
! Full multipole energy (simplified)
! ============================================================================
function compute_energy_multipole_full(sys) result(energy)
    type(particle_system), intent(inout) :: sys
    real(dp) :: energy

    ! For simplicity, use same as method 1
    energy = compute_energy_multipole_cells(sys)

end function compute_energy_multipole_full

! ============================================================================
! Energy within a cell (direct)
! ============================================================================
function compute_cell_self_energy(sys, ci, cj) result(energy)
    type(particle_system), intent(in) :: sys
    integer, intent(in) :: ci, cj
    real(dp) :: energy
    integer :: i, j, pi, pj
    real(dp) :: dx, dy, r, qi, qj
    type(cell_type) :: cell

    energy = 0.0_dp
    cell = sys%cells(ci, cj)

    do i = 1, cell%num_particles - 1
        pi = cell%particle_ids(i)
        do j = i + 1, cell%num_particles
            pj = cell%particle_ids(j)

            dx = sys%particles(pi)%x - sys%particles(pj)%x
            dy = sys%particles(pi)%y - sys%particles(pj)%y
            r = sqrt(dx*dx + dy*dy)

            if (r > 1.0e-15_dp) then
                qi = sys%particles(pi)%charge
                qj = sys%particles(pj)%charge
                energy = energy + qi * qj / r
            end if
        end do
    end do

end function compute_cell_self_energy

! ============================================================================
! Energy between two cells (direct)
! ============================================================================
function compute_cell_pair_energy_direct(sys, ci1, cj1, ci2, cj2) result(energy)
    type(particle_system), intent(in) :: sys
    integer, intent(in) :: ci1, cj1, ci2, cj2
    real(dp) :: energy
    integer :: i, j, pi, pj
    real(dp) :: dx, dy, r, qi, qj

    energy = 0.0_dp

    do i = 1, sys%cells(ci1, cj1)%num_particles
        pi = sys%cells(ci1, cj1)%particle_ids(i)
        qi = sys%particles(pi)%charge

        do j = 1, sys%cells(ci2, cj2)%num_particles
            pj = sys%cells(ci2, cj2)%particle_ids(j)
            qj = sys%particles(pj)%charge

            dx = sys%particles(pi)%x - sys%particles(pj)%x
            dy = sys%particles(pi)%y - sys%particles(pj)%y
            r = sqrt(dx*dx + dy*dy)

            if (r > 1.0e-15_dp) then
                energy = energy + 0.5_dp * qi * qj / r  ! Factor 0.5 to avoid double counting
            end if
        end do
    end do

end function compute_cell_pair_energy_direct

! ============================================================================
! Energy between two cells (multipole)
! ============================================================================
function compute_cell_pair_energy_multipole(sys, ci1, cj1, ci2, cj2) result(energy)
    type(particle_system), intent(in) :: sys
    integer, intent(in) :: ci1, cj1, ci2, cj2
    real(dp) :: energy
    integer :: n, angle_idx
    real(dp) :: dx, dy, r, cos_theta, P_n
    type(cell_type) :: cell1, cell2

    energy = 0.0_dp

    cell1 = sys%cells(ci1, cj1)
    cell2 = sys%cells(ci2, cj2)

    ! Distance between cell centers
    dx = cell2%center_x - cell1%center_x
    dy = cell2%center_y - cell1%center_y
    r = sqrt(dx*dx + dy*dy)

    if (r < 1.0e-10_dp) return

    ! Angle
    cos_theta = dx / r
    call find_legendre_value(cos_theta, angle_idx)

    ! Multipole expansion: E = sum_n M1_n * M2_n * P_n(cos_theta) / r^(n+1)
    do n = 0, sys%p_order
        P_n = legendre_table(n, angle_idx)
        energy = energy + cell1%moments(n) * cell2%moments(n) * P_n / r**(n+1)
    end do

    energy = 0.5_dp * energy  ! Factor to avoid double counting

end function compute_cell_pair_energy_multipole

! ============================================================================
! Find nearest Legendre table index for cos(theta)
! ============================================================================
subroutine find_legendre_value(cos_theta, idx)
    real(dp), intent(in) :: cos_theta
    integer, intent(out) :: idx
    real(dp) :: theta

    ! Convert cos_theta to angle
    theta = acos(max(-1.0_dp, min(1.0_dp, cos_theta)))

    ! Find index in table
    idx = int(theta / (4.0_dp * atan(1.0_dp)) * real(num_angles-1, dp)) + 1
    idx = max(1, min(num_angles, idx))

end subroutine find_legendre_value

! ============================================================================
! Get force on single particle
! ============================================================================
subroutine get_force(sys, pid, fx, fy, method)
    type(particle_system), intent(inout) :: sys
    integer, intent(in) :: pid
    real(dp), intent(out) :: fx, fy
    integer, intent(in) :: method

    select case(method)
    case(0)
        call compute_force_direct(sys, pid, fx, fy)
    case(1, 2)
        call compute_force_multipole(sys, pid, fx, fy)
    case default
        fx = 0.0_dp
        fy = 0.0_dp
    end select

end subroutine get_force

! ============================================================================
! Direct force calculation
! ============================================================================
subroutine compute_force_direct(sys, pid, fx, fy)
    type(particle_system), intent(in) :: sys
    integer, intent(in) :: pid
    real(dp), intent(out) :: fx, fy
    integer :: j
    real(dp) :: dx, dy, r, r3, qi, qj, force

    fx = 0.0_dp
    fy = 0.0_dp
    qi = sys%particles(pid)%charge

    do j = 1, sys%num_particles
        if (j == pid) cycle

        dx = sys%particles(pid)%x - sys%particles(j)%x
        dy = sys%particles(pid)%y - sys%particles(j)%y
        r = sqrt(dx*dx + dy*dy)

        if (r > 1.0e-15_dp) then
            qj = sys%particles(j)%charge
            r3 = r * r * r
            force = qi * qj / r3

            fx = fx + force * dx
            fy = fy + force * dy
        end if
    end do

end subroutine compute_force_direct

! ============================================================================
! Force with multipole expansion
! ============================================================================
subroutine compute_force_multipole(sys, pid, fx, fy)
    type(particle_system), intent(inout) :: sys
    integer, intent(in) :: pid
    real(dp), intent(out) :: fx, fy

    ! Simplified: use direct for now
    call compute_force_direct(sys, pid, fx, fy)

end subroutine compute_force_multipole

! ============================================================================
! Compute forces for all particles
! ============================================================================
subroutine get_all_forces(sys, method)
    type(particle_system), intent(inout) :: sys
    integer, intent(in) :: method
    integer :: i

    do i = 1, sys%num_particles
        call get_force(sys, i, sys%particles(i)%fx, sys%particles(i)%fy, method)
    end do

end subroutine get_all_forces

! ============================================================================
! Move particles using Euler or RK4
! ============================================================================
subroutine move_particles(sys, eta, dt, use_rk4)
    type(particle_system), intent(inout) :: sys
    real(dp), intent(in) :: eta, dt
    logical, intent(in) :: use_rk4

    if (use_rk4) then
        call move_particles_rk4(sys, eta, dt)
    else
        call move_particles_euler(sys, eta, dt)
    end if

end subroutine move_particles

! ============================================================================
! Euler method: dr/dt = eta * F
! ============================================================================
subroutine move_particles_euler(sys, eta, dt)
    type(particle_system), intent(inout) :: sys
    real(dp), intent(in) :: eta, dt
    integer :: i
    real(dp) :: new_x, new_y

    do i = 1, sys%num_particles
        ! Update position
        new_x = sys%particles(i)%x + eta * dt * sys%particles(i)%fx
        new_y = sys%particles(i)%y + eta * dt * sys%particles(i)%fy

        ! Apply boundary conditions (reflection)
        call apply_boundaries(sys, new_x, new_y)

        sys%particles(i)%x = new_x
        sys%particles(i)%y = new_y
    end do

end subroutine move_particles_euler

! ============================================================================
! RK4 method
! ============================================================================
subroutine move_particles_rk4(sys, eta, dt)
    type(particle_system), intent(inout) :: sys
    real(dp), intent(in) :: eta, dt
    integer :: i
    real(dp) :: x0, y0, fx, fy
    real(dp) :: k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y
    real(dp) :: new_x, new_y

    do i = 1, sys%num_particles
        x0 = sys%particles(i)%x
        y0 = sys%particles(i)%y

        ! k1
        k1x = eta * sys%particles(i)%fx
        k1y = eta * sys%particles(i)%fy

        ! k2 (simplified, using same forces)
        k2x = k1x
        k2y = k1y

        ! k3
        k3x = k1x
        k3y = k1y

        ! k4
        k4x = k1x
        k4y = k1y

        ! Update position
        new_x = x0 + dt/6.0_dp * (k1x + 2.0_dp*k2x + 2.0_dp*k3x + k4x)
        new_y = y0 + dt/6.0_dp * (k1y + 2.0_dp*k2y + 2.0_dp*k3y + k4y)

        ! Apply boundaries
        call apply_boundaries(sys, new_x, new_y)

        sys%particles(i)%x = new_x
        sys%particles(i)%y = new_y
    end do

end subroutine move_particles_rk4

! ============================================================================
! Apply reflecting boundary conditions
! ============================================================================
subroutine apply_boundaries(sys, x, y)
    type(particle_system), intent(in) :: sys
    real(dp), intent(inout) :: x, y

    ! X boundaries
    if (x < -sys%xmax) x = -2.0_dp * sys%xmax - x
    if (x >  sys%xmax) x =  2.0_dp * sys%xmax - x

    ! Y boundaries
    if (y < -sys%ymax) y = -2.0_dp * sys%ymax - y
    if (y >  sys%ymax) y =  2.0_dp * sys%ymax - y

end subroutine apply_boundaries

! ============================================================================
! Correct cell assignments after movement
! ============================================================================
subroutine correct_cells(sys)
    type(particle_system), intent(inout) :: sys
    integer :: i, new_i, new_j, old_i, old_j

    ! Clear all cells
    do new_j = 1, sys%num_cells_y
        do new_i = 1, sys%num_cells_x
            sys%cells(new_i, new_j)%num_particles = 0
        end do
    end do

    ! Reassign particles
    do i = 1, sys%num_particles
        call get_cell_indices(sys, sys%particles(i)%x, sys%particles(i)%y, new_i, new_j)
        sys%particles(i)%cell_i = new_i
        sys%particles(i)%cell_j = new_j
        call add_to_cell(sys%cells(new_i, new_j), i)
    end do

end subroutine correct_cells

! ============================================================================
! Determine series order based on precision
! ============================================================================
subroutine determine_series_order(sys, precision)
    type(particle_system), intent(inout) :: sys
    real(dp), intent(in) :: precision
    integer :: order

    ! Empirical formula: higher precision needs higher order
    if (precision < 1.0e-12_dp) then
        order = 25
    else if (precision < 1.0e-9_dp) then
        order = 20
    else if (precision < 1.0e-6_dp) then
        order = 15
    else if (precision < 1.0e-3_dp) then
        order = 10
    else
        order = 5
    end if

    order = min(order, p_max)
    sys%p_order = order

    print *, "Set multipole order to", order, "for precision", precision

end subroutine determine_series_order

! ============================================================================
! Write frame to CSV file
! ============================================================================
subroutine write_frame(sys, frame_num)
    type(particle_system), intent(in) :: sys
    integer, intent(in) :: frame_num
    character(len=100) :: filename
    integer :: i, unit

    write(filename, '(A,I6.6,A)') 'frame_', frame_num, '.csv'

    open(newunit=unit, file=trim(filename), status='replace')
    write(unit, '(A)') 'x,y,charge'

    do i = 1, sys%num_particles
        write(unit, '(F15.8,A,F15.8,A,F15.8)') &
            sys%particles(i)%x, ',', sys%particles(i)%y, ',', sys%particles(i)%charge
    end do

    close(unit)

end subroutine write_frame

end module multipole_module
