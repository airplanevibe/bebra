program main
  use types_mod
  use potential_interface_mod
  use fmm_2d_simple_mod
  use md_mod
  implicit none

  integer :: N, nsteps, nframes
  real(dp) :: dt, Lbox
  real(dp), allocatable :: xs(:), ys(:), zs(:), qs(:)
  real(dp), allocatable :: vxs(:), vys(:)
  real(dp), allocatable :: fx(:), fy(:), fz(:), phi(:)
  real(dp), allocatable :: fx_dir(:), fy_dir(:), fz_dir(:), phi_dir(:)
  integer :: step, i, j, ip
  real(dp) :: dump_interval_real
  integer :: dump_interval
  real(dp) :: t0f, t1f, t_fmm, t_dir
  real(dp) :: rx, ry, rz, r2, r_inv
  real(dp) :: max_err, rms_err, err, max_phi_err, rms_phi_err, phi_err
  integer :: imax, imax_phi
  integer :: p_order, ngrid
  real(dp) :: speedup

  print *, "=============================================="
  print *, "  2D Cartesian Multipole FMM for Superconductors"
  print *, "  Grid-based Fast Multipole Method"
  print *, "=============================================="
  print *, ""

  ! Configuration
  N = 5000              ! Number of particles
  dt = 0.001_dp
  Lbox = 1.0_dp
  nsteps = 3            ! Quick test
  nframes = 2

  ! FMM parameters - configurable!
  p_order = 6           ! Multipole expansion order (try 2, 4, 6, 8, 10)
  ngrid = 40            ! Grid size (ngrid x ngrid cells)

  print *, "Configuration:"
  print '(A,I0)', "  Number of particles (N): ", N
  print '(A,ES10.3)', "  Time step (dt): ", dt
  print '(A,F6.2)', "  Box size (Lbox): ", Lbox
  print '(A,I0)', "  MD steps: ", nsteps
  print '(A,I0)', "  Output frames: ", nframes
  print *, ""
  print *, "FMM Parameters:"
  print '(A,I0)', "  Multipole expansion order (p): ", p_order
  print '(A,I0)', "  Grid size: ", ngrid, " x ", ngrid
  print '(A,I0)', "  Total grid cells: ", ngrid * ngrid
  print '(A,F6.2)', "  Particles per cell (avg): ", real(N, dp) / real(ngrid*ngrid, dp)
  print *, ""

  allocate(xs(N), ys(N), zs(N), qs(N))
  allocate(vxs(N), vys(N))
  allocate(fx(N), fy(N), fz(N), phi(N))
  allocate(fx_dir(N), fy_dir(N), fz_dir(N), phi_dir(N))

  call set_coulomb_potential()
  print *, "Potential: Coulomb (1/r) for z=0 (2D)"
  print *, ""

  ! Initialize particles with random positions
  call random_seed()
  call random_number(xs)
  call random_number(ys)
  call random_number(zs)

  xs = (xs - 0.5_dp) * 2.0_dp * Lbox * 0.9_dp  ! Keep inside domain
  ys = (ys - 0.5_dp) * 2.0_dp * Lbox * 0.9_dp
  zs = 0.0_dp
  qs = 1.0_dp
  vxs = 0.0_dp
  vys = 0.0_dp

  dump_interval_real = real(nsteps, dp) / real(max(1, nframes - 1), dp)
  dump_interval = max(1, int(dump_interval_real))

  call write_frame_csv(0, xs, ys, N)

  print *, "Starting MD simulation with FMM..."
  print *, ""

  do step = 1, nsteps
    print '(A,I0,A,I0)', "Step ", step, " / ", nsteps

    ! FMM computation
    call cpu_time(t0f)
    call fmm_2d_compute(N, xs, ys, zs, qs, p_order, ngrid, phi, fx, fy, fz)
    call cpu_time(t1f)
    t_fmm = t1f - t0f
    print '(A,F10.6,A)', "  FMM time: ", t_fmm, " s"

    ! Verification with direct sum (only first step)
    if (step == 1) then
      print *, ""
      print *, "  ================================================"
      print *, "  ACCURACY VERIFICATION (Direct Sum Comparison)"
      print *, "  ================================================"
      print *, ""

      fx_dir = 0.0_dp
      fy_dir = 0.0_dp
      fz_dir = 0.0_dp
      phi_dir = 0.0_dp

      call cpu_time(t0f)
      !$OMP PARALLEL DO PRIVATE(i, j, rx, ry, rz, r2, r_inv) REDUCTION(+:fx_dir, fy_dir, phi_dir)
      do i = 1, N
        do j = 1, N
          if (i == j) cycle

          rx = xs(i) - xs(j)
          ry = ys(i) - ys(j)
          rz = zs(i) - zs(j)
          r2 = rx*rx + ry*ry + rz*rz

          if (r2 < 1e-18_dp) cycle

          r_inv = 1.0_dp / sqrt(r2)
          phi_dir(i) = phi_dir(i) + qs(j) * r_inv
          fx_dir(i) = fx_dir(i) + qs(j) * rx * r_inv / r2
          fy_dir(i) = fy_dir(i) + qs(j) * ry * r_inv / r2
        end do
      end do
      !$OMP END PARALLEL DO
      call cpu_time(t1f)
      t_dir = t1f - t0f

      speedup = t_dir / t_fmm

      print '(A,F10.6,A)', "  Direct sum time: ", t_dir, " s"
      print '(A,F10.2,A)', "  Speedup: ", speedup, "x"
      print *, ""

      ! Sample comparison (first 10 particles)
      print *, "  Sample comparison (first 10 particles):"
      print *, "  ----------------------------------------------------------------"
      print *, "  i  |   Fx(FMM)   |  Fx(Direct) |   Fy(FMM)   |  Fy(Direct)  "
      print *, "  ----------------------------------------------------------------"
      do i = 1, min(10, N)
        print '(A,I3,A,4ES14.6)', "  ", i, " ", fx(i), fx_dir(i), fy(i), fy_dir(i)
      end do
      print *, "  ----------------------------------------------------------------"
      print *, ""

      print *, "  Potential comparison (first 10 particles):"
      print *, "  --------------------------------------------------"
      print *, "  i  |   Phi(FMM)    |   Phi(Direct)   "
      print *, "  --------------------------------------------------"
      do i = 1, min(10, N)
        print '(A,I3,A,2ES16.8)', "  ", i, " ", phi(i), phi_dir(i)
      end do
      print *, "  --------------------------------------------------"
      print *, ""

      ! Error statistics for forces
      max_err = 0.0_dp
      rms_err = 0.0_dp
      imax = 1

      do i = 1, N
        err = sqrt((fx(i)-fx_dir(i))**2 + (fy(i)-fy_dir(i))**2)
        rms_err = rms_err + err*err
        if (err > max_err) then
          max_err = err
          imax = i
        end if
      end do

      rms_err = sqrt(rms_err / real(N, dp))

      ! Error statistics for potential
      max_phi_err = 0.0_dp
      rms_phi_err = 0.0_dp
      imax_phi = 1

      do i = 1, N
        phi_err = abs(phi(i) - phi_dir(i))
        rms_phi_err = rms_phi_err + phi_err*phi_err
        if (phi_err > max_phi_err) then
          max_phi_err = phi_err
          imax_phi = i
        end if
      end do

      rms_phi_err = sqrt(rms_phi_err / real(N, dp))

      ! Calculate relative errors
      print *, "  ================================================"
      print *, "  ACCURACY STATISTICS"
      print *, "  ================================================"
      print *, ""
      print *, "  Force accuracy:"
      print '(A,ES12.4)', "    Max absolute error: ", max_err
      print '(A,I0)', "    (at particle ", imax, ")"
      print '(A,ES12.4)', "    RMS error: ", rms_err

      ! Compute typical force magnitude for relative error
      if (imax >= 1 .and. imax <= N) then
        print '(A,ES12.4)', "    Force magnitude at max error: ", &
              sqrt(fx_dir(imax)**2 + fy_dir(imax)**2)
        if (sqrt(fx_dir(imax)**2 + fy_dir(imax)**2) > 1e-10_dp) then
          print '(A,ES12.4)', "    Relative error: ", &
                max_err / sqrt(fx_dir(imax)**2 + fy_dir(imax)**2)
        end if
      end if
      print *, ""

      print *, "  Potential accuracy:"
      print '(A,ES12.4)', "    Max absolute error: ", max_phi_err
      print '(A,I0)', "    (at particle ", imax_phi, ")"
      print '(A,ES12.4)', "    RMS error: ", rms_phi_err

      if (imax_phi >= 1 .and. imax_phi <= N) then
        print '(A,ES12.4)', "    Potential at max error: ", phi_dir(imax_phi)
        if (abs(phi_dir(imax_phi)) > 1e-10_dp) then
          print '(A,ES12.4)', "    Relative error: ", max_phi_err / abs(phi_dir(imax_phi))
        end if
      end if
      print *, ""

      print *, "  ================================================"
      print *, "  PERFORMANCE SUMMARY"
      print *, "  ================================================"
      print *, ""
      print '(A,I0)', "    N particles: ", N
      print '(A,I0)', "    Expansion order (p): ", p_order
      print '(A,I0)', "    Grid size: ", ngrid, "x", ngrid
      print '(A,F10.6,A)', "    FMM time: ", t_fmm, " s"
      print '(A,F10.6,A)', "    Direct time: ", t_dir, " s"
      print '(A,F10.2,A)', "    Speedup: ", speedup, "x"
      print '(A,ES10.2,A)', "    FMM complexity: O(N) ~ ", t_fmm/real(N, dp)*1e6, " μs/particle"
      print '(A,ES10.2,A)', "    Direct complexity: O(N²) ~ ", t_dir/real(N, dp)**2*1e9, " ns/pair"
      print *, ""
      print *, "  ================================================"
      print *, ""

      ! Check if accuracy goal is met
      if (rms_err < 1.0e-8_dp .and. rms_phi_err < 1.0e-8_dp) then
        print *, "  ✓ ACCURACY GOAL MET: RMS errors < 1e-8"
      else if (rms_err < 1.0e-6_dp .and. rms_phi_err < 1.0e-6_dp) then
        print *, "  ✓ Good accuracy: RMS errors < 1e-6"
      else
        print *, "  ! Consider increasing p_order for better accuracy"
      end if

      if (speedup > 1.0_dp) then
        print '(A,F6.2,A)', "  ✓ FMM is faster: ", speedup, "x speedup"
      else
        print *, "  ! FMM not faster - try larger N or smaller ngrid"
      end if
      print *, ""

    end if

    ! MD integration
    call md_step_euler(xs, ys, vxs, vys, fx, fy, N, dt)
    call apply_reflect_bc(xs, ys, vxs, vys, N, Lbox)

    ! Write output
    if (mod(step, dump_interval) == 0 .or. step == nsteps) then
      call write_frame_csv(step, xs, ys, N)
      print '(A,I5.5,A)', "  Wrote frame_", step, ".csv"
    end if

    print *, ""
  end do

  print *, "=============================================="
  print *, "Simulation completed successfully!"
  print *, "=============================================="

  deallocate(xs, ys, zs, qs, vxs, vys, fx, fy, fz, phi)
  deallocate(fx_dir, fy_dir, fz_dir, phi_dir)

end program main
