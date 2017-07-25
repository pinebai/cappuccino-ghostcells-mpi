!***********************************************************************
!
subroutine update_values_at_ghostcells
!
! Ghost cells are inner cells here - the first layer of inner cells and
! last layer of inner cells. By inner cells I mean those that participate
! in usual loops over cells.
!
!
!  Periodic BC with ghost cells:
!      _________________________________________
!     |                                         |
!     v                                         |
!   ingcw    inp                               ine    ingce
! |xxxoxxx|---o---|---o---|-- ... --|---o---|---o---|xxxoxxx|
!             |                                         ^
!             !_________________________________________|
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use variables
  use boundc

  implicit none

!
! Local
!
  integer :: i, j, k, inp, ine, ingce, ingcw,ingc
  real(dp) :: fxp,fxe
!
!***********************************************************************
!


! WEST-EAST direction
  do k=3,nkmm
  do j=3,njmm
   
  ingcw=lk(k)+li(2)+j   
  inp=lk(k)+li(3)+j

  ine=lk(k)+li(nimm)+j
  ingce=lk(k)+li(nim)+j


  u(ingcw)   = u(ine)
  v(ingcw)   = v(ine)
  w(ingcw)   = w(ine)
  te(ingcw)  = te(ine)
  ed(ingcw)  = ed(ine)
  vis(ingcw) = vis(ine)

  u(ingce)   = u(inp)
  v(ingce)   = v(inp)
  w(ingce)   = w(inp)
  te(ingce)  = te(inp)
  ed(ingce)  = ed(inp)
  vis(ingce) = vis(inp)

  enddo
  enddo

! SOUTH-NORTH direction
  do k=3,nkmm
  do i=3,nimm

  ingcw=lk(k)+li(i)+2   
  inp=lk(k)+li(i)+3

  ine=lk(k)+li(i)+njmm
  ingce=lk(k)+li(i)+njm


  u(ingcw)   = u(ine)
  v(ingcw)   = v(ine)
  w(ingcw)   = w(ine)
  te(ingcw)  = te(ine)
  ed(ingcw)  = ed(ine)
  vis(ingcw) = vis(ine)


  u(ingce)   = u(inp)
  v(ingce)   = v(inp)
  w(ingce)   = w(inp)
  te(ingce)  = te(inp)
  ed(ingce)  = ed(inp)
  vis(ingce) = vis(inp)

  enddo
  enddo  

! BOTTOM plane
  do i=3,nimm
  do j=3,njmm

  ingc=lk(2)+li(i)+j   
  inp=lk(3)+li(i)+j

  ! No slip at the wall
  ine=inp-nij
  fxp=fz(ine) 
  fxe=1.0_dp-fxp

  u(ingc) = - u(inp)*fxp/fxe
  v(ingc) = - v(inp)*fxp/fxe
  w(ingc) = - w(inp)*fxp/fxe
  te(ingc) = - te(inp)*fxp/fxe

  ed(ingc)  = ed(inp)

  vis(ingc) = ( vis(ingc) - vis(inp)*fxp ) / fxe

  enddo
  enddo 




! TOP plane
  do i=3,nimm
  do j=3,njmm

  ingc=lk(nkm)+li(i)+j   
  inp=lk(nkmm)+li(i)+j

  ! No slip at the wall
  ! ine=inp+nij 
  fxe=fz(inp) 
  fxp=1.0_dp-fxe

  u(ingc) = - u(inp)*fxp/fxe
  v(ingc) = - v(inp)*fxp/fxe
  w(ingc) = - w(inp)*fxp/fxe
  te(ingc) = - te(inp)*fxp/fxe

  ed(ingc)  = ed(inp)

  vis(ingc) = ( vis(ingc) - vis(inp)*fxp ) / fxe

  enddo
  enddo 



end subroutine