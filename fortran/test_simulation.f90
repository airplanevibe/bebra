! ============================================================================
! Test program for multipole particle simulation
! ============================================================================
program test_simulation
    use legendre_data, only: dp
    use multipole_module
    implicit none

    type(particle_system) :: sys
    integer :: i, step, num_steps, frame_interval
    real(dp) :: energy, dt, eta, precision
    integer :: method, num_particles
    real(dp) :: x, y, charge
    logical :: use_rk4
    integer :: unit

    print *, "============================================"
    print *, "  Multipole Particle Simulation Test"
    print *, "============================================"
    print *

    ! Initialize system
    call init_system(sys, 1.0_dp, 1.0_dp, 10, 10)

    ! Set precision and determine order
    precision = 1.0e-6_dp
    call determine_series_order(sys, precision)

    ! Add particles in a grid pattern
    num_particles = 50
    print *, "Adding", num_particles, "particles..."

    do i = 1, num_particles
        ! Random position
        call random_number(x)
        call random_number(y)
        x = 1.6_dp * x - 0.8_dp  ! [-0.8, 0.8]
        y = 1.6_dp * y - 0.8_dp

        ! Alternating charges
        if (mod(i, 2) == 0) then
            charge = 1.0_dp
        else
            charge = -1.0_dp
        end if

        call add_particle(sys, charge, x, y)
    end do

    print *, "Total particles:", sys%num_particles
    print *

    ! Simulation parameters
    dt = 0.01_dp
    eta = 0.1_dp
    num_steps = 200
    frame_interval = 5
    method = 1  ! 0=direct, 1=multipole cells, 2=full multipole
    use_rk4 = .false.  ! Start with Euler

    print *, "Simulation parameters:"
    print *, "  dt =", dt
    print *, "  eta =", eta
    print *, "  num_steps =", num_steps
    print *, "  method =", method
    if (use_rk4) then
        print *, "  integrator = RK4"
    else
        print *, "  integrator = Euler"
    end if
    print *

    ! Open energy log
    open(newunit=unit, file='energy.log', status='replace')
    write(unit, '(A)') '# Step  Time  Energy'

    ! Write initial frame
    call write_frame(sys, 0)

    ! Main simulation loop
    print *, "Running simulation..."
    print *, "Step    Time      Energy      Method"
    print *, "----------------------------------------"

    do step = 1, num_steps
        ! Compute forces
        call get_all_forces(sys, method)

        ! Move particles
        call move_particles(sys, eta, dt, use_rk4)

        ! Correct cell assignments
        call correct_cells(sys)

        ! Compute energy
        energy = get_energy(sys, method)

        ! Log energy
        write(unit, '(I6,2X,F10.4,2X,ES15.8)') step, step*dt, energy

        ! Print progress
        if (mod(step, 10) == 0) then
            print '(I6,2X,F10.4,2X,ES15.8,2X,I3)', step, step*dt, energy, method
        end if

        ! Write frame
        if (mod(step, frame_interval) == 0) then
            call write_frame(sys, step / frame_interval)
        end if

        ! Switch to RK4 after some steps
        if (step == 50 .and. .not. use_rk4) then
            use_rk4 = .true.
            print *, "  -> Switched to RK4 integrator"
        end if
    end do

    close(unit)

    ! Write final frame
    call write_frame(sys, num_steps / frame_interval + 1)

    print *
    print *, "============================================"
    print *, "  Simulation completed!"
    print *, "============================================"
    print *, "Output files:"
    print *, "  - frame_*.csv (", num_steps/frame_interval + 2, "frames )"
    print *, "  - energy.log"
    print *
    print *, "Visualize with:"
    print *, "  python3 visualize.py"
    print *, "  python3 visualize.py --energy"
    print *

    ! Test different methods
    print *, "Comparing methods:"
    print *, "  Method 0 (direct):", get_energy(sys, 0)
    print *, "  Method 1 (multipole cells):", get_energy(sys, 1)
    print *, "  Method 2 (full multipole):", get_energy(sys, 2)
    print *

end program test_simulation
