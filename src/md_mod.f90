module md_mod
  use types_mod
  implicit none

  public :: md_step_euler, apply_reflect_bc, write_frame_csv

contains

  subroutine md_step_euler(xs, ys, vxs, vys, fx, fy, n, dt)
    real(dp), intent(inout) :: xs(:), ys(:), vxs(:), vys(:)
    real(dp), intent(in) :: fx(:), fy(:)
    integer, intent(in) :: n
    real(dp), intent(in) :: dt
    integer :: i

    do i = 1, n
      vxs(i) = vxs(i) + fx(i) * dt
      vys(i) = vys(i) + fy(i) * dt
      xs(i) = xs(i) + vxs(i) * dt
      ys(i) = ys(i) + vys(i) * dt
    end do
  end subroutine md_step_euler

  subroutine apply_reflect_bc(xs, ys, vxs, vys, n, L)
    real(dp), intent(inout) :: xs(:), ys(:), vxs(:), vys(:)
    integer, intent(in) :: n
    real(dp), intent(in) :: L
    integer :: i

    do i = 1, n
      if (xs(i) > L) then
        xs(i) = 2.0_dp*L - xs(i)
        vxs(i) = -vxs(i)
      else if (xs(i) < -L) then
        xs(i) = -2.0_dp*L - xs(i)
        vxs(i) = -vxs(i)
      end if

      if (ys(i) > L) then
        ys(i) = 2.0_dp*L - ys(i)
        vys(i) = -vys(i)
      else if (ys(i) < -L) then
        ys(i) = -2.0_dp*L - ys(i)
        vys(i) = -vys(i)
      end if
    end do
  end subroutine apply_reflect_bc

  subroutine write_frame_csv(step, xs, ys, n)
    integer, intent(in) :: step, n
    real(dp), intent(in) :: xs(:), ys(:)
    character(len=64) :: fname
    integer :: i, unit

    write(fname,'(A,I5.5,A)') 'frame_', step, '.csv'
    open(newunit=unit, file=fname, status='replace', action='write')
    write(unit, '(A)') 'x,y'
    do i = 1, n
      write(unit,'(F12.6,",",F12.6)') xs(i), ys(i)
    end do
    close(unit)
  end subroutine write_frame_csv

end module md_mod
