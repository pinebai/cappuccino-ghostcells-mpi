subroutine update_velocity_rk4(uc,vc,wc,istage)
  use types
  use parameters
  use variables
  use coef
  use time_mod

  implicit none

  integer, intent(in) :: istage
  real(prec), dimension(nxyza), intent(inout) :: Uc,Vc,Wc

  ! Forth-order Runge-Kutta coefficients
  real(prec), parameter :: a1_rk = 1./6.0_dp
  real(prec), parameter :: a2_rk = 1./3.0_dp
  real(prec), parameter :: a3_rk = 1./3.0_dp
  real(prec), parameter :: a4_rk = 1./6.0_dp

  real(prec) :: dt

  dt = timestep

  if(istage == 1) then
 
    ! U = Uo + dt/den * su
    ! V = Vo + dt/den * sv
    ! W = Wo + dt/den * sw

    Uc = Uc + a1_rk * dt/den * su
    Vc = Vc + a1_rk * dt/den * sv
    Wc = Wc + a1_rk * dt/den * sw

    U = Uo + 0.5 * dt/den * su
    V = Vo + 0.5 * dt/den * sv
    W = Wo + 0.5 * dt/den * sw

  elseif(istage == 2) then

    Uc = Uc + a2_rk * dt/den * su
    Vc = Vc + a2_rk * dt/den * sv
    Wc = Wc + a2_rk * dt/den * sw

    U = Uo + 0.5 * dt/den * su
    V = Vo + 0.5 * dt/den * sv
    W = Wo + 0.5 * dt/den * sw

  elseif(istage == 3) then

    Uc = Uc + a3_rk * dt/den * su
    Vc = Vc + a3_rk * dt/den * sv
    Wc = Wc + a3_rk * dt/den * sw

    U = Uo + dt/den * su
    V = Vo + dt/den * sv
    W = Wo + dt/den * sw

  else ! if(istage == 4) then

    Uc = Uc + a4_rk * dt/den * su
    Vc = Vc + a4_rk * dt/den * sv
    Wc = Wc + a4_rk * dt/den * sw

    U = Uc
    V = Vc
    W = Wc
    
  endif

end subroutine