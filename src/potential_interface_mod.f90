module potential_interface_mod
  use types_mod
  implicit none

  ! Abstract interface for different potentials (Coulomb, Yukawa, etc.)
  abstract interface
     pure function potential_func(r, q) result(phi)
       import :: dp
       real(dp), intent(in) :: r, q
       real(dp) :: phi
     end function potential_func

     subroutine force_func(rx, ry, rz, q, fx, fy, fz)
       import :: dp
       real(dp), intent(in) :: rx, ry, rz, q
       real(dp), intent(out) :: fx, fy, fz
     end subroutine force_func
  end interface

  procedure(potential_func), pointer :: active_potential => null()
  procedure(force_func), pointer :: active_force => null()

contains

  subroutine set_coulomb_potential()
    active_potential => coulomb_potential
    active_force => coulomb_force
  end subroutine set_coulomb_potential

  pure function coulomb_potential(r, q) result(phi)
    real(dp), intent(in) :: r, q
    real(dp) :: phi

    if (r < 1e-14_dp) then
      phi = 0.0_dp
    else
      phi = q / r
    endif
  end function coulomb_potential

  subroutine coulomb_force(rx, ry, rz, q, fx, fy, fz)
    real(dp), intent(in) :: rx, ry, rz, q
    real(dp), intent(out) :: fx, fy, fz
    real(dp) :: r2, rinv

    r2 = rx*rx + ry*ry + rz*rz
    if (r2 < 1e-16_dp) then
      fx = 0.0_dp; fy = 0.0_dp; fz = 0.0_dp
      return
    endif

    rinv = 1.0_dp / sqrt(r2)
    fx = q * rx * rinv / r2
    fy = q * ry * rinv / r2
    fz = q * rz * rinv / r2
  end subroutine coulomb_force

end module potential_interface_mod
