!***********************************************************************
!
subroutine update_values_at_ghostcells_nonperiodic
!
! Ghost cells are inner cells here - the first layer of inner cells and
! last layer of inner cells. By inner cells I mean those that participate
! in usual loops over cells.
!
!  Periodic BC with ghost cells:
!
!        ingc    inp                                            
!     |xxxoxxx|---o---|---o---|-- ... --|---o---|---o---|xxxoxxx|
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
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use variables
  use bc
  use inlet

  implicit none

!
! Local
!
  integer :: i, j, k, inp, ingc, inbc, iface
  real(dp) :: fxp,fxe,fxer
  real(dp) :: Unmag
!
!***********************************************************************
!

!
! west
  do k=3,nkmm
  do j=3,njmm

  ! Index of currecnt inner cell  
  inp=lk(k)+li(3)+j

  ! Index of adjecent ghost cell
  ingc = inp-nj

  ! Index of a face dividing this cell with ghost cell
  iface = inp-nj

  ! Index for boundary condition identifier array
  inbc=(j-1)*nk+k

  ! Boundary condition identifier lbhlp
  lbhlp=lbw(inbc)

  ! Interpolation factor for current inner cell value - fxp
  fxp = fx(inp-nj)

  ! Interpolation factor for adjecent ghost cell and its reciprocal value
  fxe = 1.0_dp-fxp
  fxer = 1.0_dp / fxe

  if(lbhlp.eq.1) then
  ! inlet

    u(ingc)  = (uin  - u(inp)*fxp) * fxer
    v(ingc)  = (vin  - v(inp)*fxp) * fxer
    w(ingc)  = (win  - w(inp)*fxp) * fxer

    te(ingc) = (tein - te(inp)*fxp) * fxer
    ed(ingc) = (edin - ed(inp)*fxp) * fxer
    t(ingc)  = (tin  -  t(inp)*fxp) * fxer

    vis(ingc) = 0.09_dp*te(ingc)**2/ed(ingc)

  endif

  if(lbhlp.eq.2) then
  ! outlet - extrapolate from two inside cells

    u(ingc) = u(inp)-(u(inp+nj)-u(inp))*fx(inp)
    v(ingc) = v(inp)-(v(inp+nj)-v(inp))*fx(inp)
    w(ingc) = w(inp)-(w(inp+nj)-w(inp))*fx(inp)

    te(ingc) = te(inp)-(te(inp+nj)-te(inp))*fx(inp) 
    ed(ingc) = ed(inp)-(ed(inp+nj)-ed(inp))*fx(inp)
    t(ingc) = t(inp)-(t(inp+nj)-t(inp))*fx(inp) 

    vis(ingc) = vis(inp)-(vis(inp+nj)-vis(inp))*fx(inp)

  endif

  if(lbhlp.eq.3) then
  ! symmetry
  
    Unmag = u(inp)*ar1x(iface)+v(inp)*ar1y(iface)+w(inp)*ar1z(iface)

    u(ingc) = u(inp)-Unmag*ar1x(iface)
    v(ingc) = v(inp)-Unmag*ar1y(iface)
    w(ingc) = w(inp)-Unmag*ar1z(iface)

    te(ingc) = te(inp)
    ed(ingc) = ed(inp)
    t(ingc) = t(inp)
    vis(ingc) = vis(inp)

  endif

  if(lbhlp.eq.4) then
  ! wall

    ! No-slip for velocity (zero at wall)
    u(ingc)  = ( - u(inp)*fxp) * fxer
    v(ingc)  = ( - v(inp)*fxp) * fxer
    w(ingc)  = ( - w(inp)*fxp) * fxer
    
    ! Zero for turbulence kinetic energy
    te(ingc) = ( - te(inp)*fxp) * fxer

    ed(ingc) = ed(inp) ! ovde kako?? znamo vrednost za eps na zidu samim tim je Dirichlet

    t(ingc)  = (tin  -  t(inp)*fxp) * fxer ! Ako je isothermal wall

    ! Dirichlet for effective viscosity
    vis(ingc) = ( vis(ingc) - vis(inp)*fxp ) * fxer

  endif

  end do
  end do



! east
  do k=3,nkmm
  do j=3,njmm

  ! Index of current inner cell  
  inp=lk(k)+li(nimm)+j

  ! Index of adjecent ghost cell  
  ingc=inp+nj

  ! Index of a face dividing this cell with ghost cell
  iface = inp

  ! Index for boundary condition identifier array
  inbc=(j-1)*nk+k

  ! Boundary condition identifier lbhlp
  lbhlp=lbe(inbc)

  ! Interpolation factor for current inner cell value - fxp
  fxp = 1.0_dp-fx(inp)

  ! Interpolation factor for adjecent ghost cell and its reciprocal value
  fxe = 1.0_dp-fxp
  fxer = 1.0_dp / fxe

  if(lbhlp.eq.1) then
  ! inlet

    u(ingc)  = (uin  - u(inp)*fxp) * fxer
    v(ingc)  = (vin  - v(inp)*fxp) * fxer
    w(ingc)  = (win  - w(inp)*fxp) * fxer

    te(ingc) = (tein - te(inp)*fxp) * fxer
    ed(ingc) = (edin - ed(inp)*fxp) * fxer
    t(ingc) = (tin - t(inp)*fxp) * fxer

    vis(ingc) = 0.09_dp*te(ingc)**2/ed(ingc)

  endif

  if(lbhlp.eq.2) then
  ! outlet - extrapolate from two inside cells

    u(ingc) = u(inp)+(u(inp)-u(inp-nj))*(1.0_dp-fx(inp-nj))
    v(ingc) = v(inp)+(v(inp)-v(inp-nj))*(1.0_dp-fx(inp-nj))
    w(ingc) = w(inp)+(w(inp)-w(inp-nj))*(1.0_dp-fx(inp-nj))

    te(ingc) = te(inp)+(te(inp)-te(inp-nj))*(1.0_dp-fx(inp-nj))
    ed(ingc) = ed(inp)+(ed(inp)-ed(inp-nj))*(1.0_dp-fx(inp-nj))
    t(ingc) = t(inp)+(t(inp)-t(inp-nj))*(1.0_dp-fx(inp-nj))

    vis(ingc) = vis(inp)+(vis(inp)-vis(inp-nj))*(1.0_dp-fx(inp-nj))

  endif

  if(lbhlp.eq.3) then
  ! symmetry
  
    Unmag = u(inp)*ar1x(iface)+v(inp)*ar1y(iface)+w(inp)*ar1z(iface)

    u(ingc) = u(inp)-Unmag*ar1x(iface)
    v(ingc) = v(inp)-Unmag*ar1y(iface)
    w(ingc) = w(inp)-Unmag*ar1z(iface)

    te(ingc) = te(inp)
    ed(ingc) = ed(inp)
    t(ingc) = t(inp)
    vis(ingc) = vis(inp)

  endif

  if(lbhlp.eq.4) then
  ! wall 

    ! No-slip for velocity (zero at wall)
    u(ingc)  = ( - u(inp)*fxp) * fxer
    v(ingc)  = ( - v(inp)*fxp) * fxer
    w(ingc)  = ( - w(inp)*fxp) * fxer
    
    ! Zero for turbulence kinetic energy
    te(ingc) = ( - te(inp)*fxp) * fxer

    ed(ingc) = ed(inp) ! ovde kako?? znamo vrednost za eps na zidu samim tim je Dirichlet
    
    t(ingc) = (tin - t(inp)*fxp) * fxer ! Ako je isothermal wall

    ! Dirichlet for effective viscosity
    vis(ingc) = ( vis(ingc) - vis(inp)*fxp ) * fxer

  endif

  end do
  end do




! south
  do k=3,nkmm
  do i=3,nimm

  ! Index of current inner cell  
  inp=lk(k)+li(i)+3

  ! Index of adjecent ghost cell  
  ingc=inp-1

  ! Index of a face dividing this cell with ghost cell
  iface = inp-1

  ! Index for boundary condition identifier array
  inbc=(i-1)*nk+k

  ! Boundary condition identifier lbhlp
  lbhlp=lbs(inbc)

  ! Interpolation factor for current inner cell value - fxp
  fxp = 1.0_dp-fy(inp)

  ! Interpolation factor for adjecent ghost cell and its reciprocal value
  fxe = 1.0_dp-fxp
  fxer = 1.0_dp / fxe

  if(lbhlp.eq.1) then
  ! inlet

    u(ingc)  = (uin  - u(inp)*fxp) * fxer
    v(ingc)  = (vin  - v(inp)*fxp) * fxer
    w(ingc)  = (win  - w(inp)*fxp) * fxer

    te(ingc) = (tein - te(inp)*fxp) * fxer
    ed(ingc) = (edin - ed(inp)*fxp) * fxer
    t(ingc) = (tin - t(inp)*fxp) * fxer

    vis(ingc) = 0.09_dp*te(ingc)**2/ed(ingc)

  endif

  if(lbhlp.eq.2) then
  ! outlet - extrapolate from two inside cells

    u(ingc) = u(inp)-(u(inp+1)-u(inp))*fy(inp)
    v(ingc) = v(inp)-(v(inp+1)-v(inp))*fy(inp)
    w(ingc) = w(inp)-(w(inp+1)-w(inp))*fy(inp)

    te(ingc) = te(inp)-(te(inp+1)-te(inp))*fy(inp)
    ed(ingc) = ed(inp)-(ed(inp+1)-ed(inp))*fy(inp)
    t(ingc) = t(inp)-(t(inp+1)-t(inp))*fy(inp)

    vis(ingc) = vis(inp)-(vis(inp+1)-vis(inp))*fy(inp)

  endif

  if(lbhlp.eq.3) then
  ! symmetry
  
    Unmag = u(inp)*ar2x(iface)+v(inp)*ar2y(iface)+w(inp)*ar2z(iface)

    u(ingc) = u(inp)-Unmag*ar2x(iface)
    v(ingc) = v(inp)-Unmag*ar2y(iface)
    w(ingc) = w(inp)-Unmag*ar2z(iface)

    te(ingc) = te(inp)
    ed(ingc) = ed(inp)
    t(ingc) = t(inp)
    vis(ingc) = vis(inp)

  endif

  if(lbhlp.eq.4) then
  ! wall 

    ! No-slip for velocity (zero at wall)
    u(ingc)  = ( - u(inp)*fxp) * fxer
    v(ingc)  = ( - v(inp)*fxp) * fxer
    w(ingc)  = ( - w(inp)*fxp) * fxer
    
    ! Zero for turbulence kinetic energy
    te(ingc) = ( - te(inp)*fxp) * fxer

    ed(ingc) = ed(inp) ! ovde kako?? znamo vrednost za eps na zidu samim tim je Dirichlet

    t(ingc) = (tin - t(inp)*fxp) * fxer ! Ako ej isothermal wall

    ! Dirichlet for effective viscosity
    vis(ingc) = ( vis(ingc) - vis(inp)*fxp ) * fxer

  endif

  end do
  end do




! north
  do k=3,nkmm
  do i=3,nimm
  

  ! Index of current inner cell  
  inp=lk(k)+li(i)+njmm

  ! Index of adjecent ghost cell  
  ingc=inp+1

  ! Index of a face dividing this cell with ghost cell
  iface = inp

  ! Index for boundary condition identifier array
  inbc=(i-1)*nk+k

  ! Boundary condition identifier lbhlp
  lbhlp=lbn(inbc)

  ! Interpolation factor for current inner cell value - fxp
  fxp = 1.0_dp-fy(inp)

  ! Interpolation factor for adjecent ghost cell and its reciprocal value
  fxe = 1.0_dp-fxp
  fxer = 1.0_dp / fxe

  if(lbhlp.eq.1) then
  ! inlet

    u(ingc)  = (uin  - u(inp)*fxp) * fxer
    v(ingc)  = (vin  - v(inp)*fxp) * fxer
    w(ingc)  = (win  - w(inp)*fxp) * fxer

    te(ingc) = (tein - te(inp)*fxp) * fxer
    ed(ingc) = (edin - ed(inp)*fxp) * fxer
    t(ingc) = (tin - t(inp)*fxp) * fxer

    vis(ingc) = 0.09_dp*te(ingc)**2/ed(ingc)

  endif

  if(lbhlp.eq.2) then
  ! outlet - extrapolate from two inside cells

    u(ingc) = u(inp)+(u(inp)-u(inp-1))*(1.0_dp-fy(inp-1))
    v(ingc) = v(inp)+(v(inp)-v(inp-1))*(1.0_dp-fy(inp-1))
    w(ingc) = w(inp)+(w(inp)-w(inp-1))*(1.0_dp-fy(inp-1))

    te(ingc) = te(inp)+(te(inp)-te(inp-1))*(1.0_dp-fy(inp-1))
    ed(ingc) = ed(inp)+(ed(inp)-ed(inp-1))*(1.0_dp-fy(inp-1))
    t(ingc) = t(inp)+(t(inp)-t(inp-1))*(1.0_dp-fy(inp-1))

    vis(ingc) = vis(inp)+(vis(inp)-vis(inp-1))*(1.0_dp-fy(inp-1))
  endif

  if(lbhlp.eq.3) then
  ! symmetry
  
    Unmag = u(inp)*ar2x(iface)+v(inp)*ar2y(iface)+w(inp)*ar2z(iface)

    u(ingc) = u(inp)-Unmag*ar2x(iface)
    v(ingc) = v(inp)-Unmag*ar2y(iface)
    w(ingc) = w(inp)-Unmag*ar2z(iface)

    te(ingc) = te(inp)
    ed(ingc) = ed(inp)
    t(ingc) = t(inp)
    vis(ingc) = vis(inp)

  endif

  if(lbhlp.eq.4) then
  ! wall 

    ! No-slip for velocity (zero at wall)
    u(ingc)  = ( - u(inp)*fxp) * fxer
    v(ingc)  = ( - v(inp)*fxp) * fxer
    w(ingc)  = ( - w(inp)*fxp) * fxer
    
    ! Zero for turbulence kinetic energy
    te(ingc) = ( - te(inp)*fxp) * fxer

    ed(ingc) = ed(inp) ! ovde kako?? znamo vrednost za eps na zidu samim tim je Dirichlet

    t(ingc) = (tin - t(inp)*fxp) * fxer ! Ako je isothermal wall

    ! Dirichlet for effective viscosity
    vis(ingc) = ( vis(ingc) - vis(inp)*fxp ) * fxer
    
  endif



  end do
  end do




! botom
  do i=3,nimm
  do j=3,njmm

  ! Index of current inner cell  
  inp=lk(3)+li(i)+j

  ! Index of adjecent ghost cell  
  ingc=inp-nij

  ! Index of a face dividing this cell with ghost cell
  iface = inp-nij

  ! Index for boundary condition identifier array
  inbc=(i-1)*nj+j

  ! Boundary condition identifier lbhlp
  lbhlp=lbb(inbc)

  ! Interpolation factor for current inner cell value - fxp
  fxp = fz(ingc)

  ! Interpolation factor for adjecent ghost cell and its reciprocal value
  fxe = 1.0_dp-fxp
  fxer = 1.0_dp / fxe

  if(lbhlp.eq.1) then
  ! inlet

    u(ingc)  = (uin  - u(inp)*fxp) * fxer
    v(ingc)  = (vin  - v(inp)*fxp) * fxer
    w(ingc)  = (win  - w(inp)*fxp) * fxer

    te(ingc) = (tein - te(inp)*fxp) * fxer
    ed(ingc) = (edin - ed(inp)*fxp) * fxer
    t(ingc) = (tin - t(inp)*fxp) * fxer

    vis(ingc) = 0.09_dp*te(ingc)**2/ed(ingc)

  endif

  if(lbhlp.eq.2) then
  ! outlet - extrapolate from two inside cells

    u(ingc) = u(inp)-(u(inp+nij)-u(inp))*fz(inp)
    v(ingc) = v(inp)-(v(inp+nij)-v(inp))*fz(inp)
    w(ingc) = w(inp)-(w(inp+nij)-w(inp))*fz(inp)

    te(ingc) = te(inp)-(te(inp+nij)-te(inp))*fz(inp)
    ed(ingc) = ed(inp)-(ed(inp+nij)-ed(inp))*fz(inp)
    t(ingc) = t(inp)-(t(inp+nij)-t(inp))*fz(inp)

    vis(ingc) = vis(inp)-(vis(inp+nij)-vis(inp))*fz(inp)

  endif

  if(lbhlp.eq.3) then
  ! symmetry
  
    Unmag = u(inp)*ar3x(iface)+v(inp)*ar3y(iface)+w(inp)*ar3z(iface)

    u(ingc) = u(inp)-Unmag*ar3x(iface)
    v(ingc) = v(inp)-Unmag*ar3y(iface)
    w(ingc) = w(inp)-Unmag*ar3z(iface)

    te(ingc) = te(inp)
    ed(ingc) = ed(inp)
    t(ingc) = t(inp)
    vis(ingc) = vis(inp)

  endif

  if(lbhlp.eq.4) then
  ! wall 

    ! No-slip for velocity (zero at wall)
    u(ingc)  = ( - u(inp)*fxp) * fxer
    v(ingc)  = ( - v(inp)*fxp) * fxer
    w(ingc)  = ( - w(inp)*fxp) * fxer
    
    ! Zero for turbulence kinetic energy
    te(ingc) = ( - te(inp)*fxp) * fxer

    ed(ingc) = ed(inp) ! ovde kako?? znamo vrednost za eps na zidu samim tim je Dirichlet

    t(ingc) = (tin - t(inp)*fxp) * fxer ! Ako je isothermal wall

    ! Dirichlet for effective viscosity
    ! vis(ingc) = ( vis(ingc) - vis(inp)*fxp ) * fxer
    vis(ingc) = vis(inp)

  endif


  end do
  end do




! top
  do i=3,nimm
  do j=3,njmm

  ! Index of current inner cell  
  inp=lk(nkmm)+li(i)+j

  ! Index of adjecent ghost cell  
  ingc=inp+nij

  ! Index of a face dividing this cell with ghost cell
  iface = inp

  ! Index for boundary condition identifier array
  inbc=(i-1)*nj+j

  ! Boundary condition identifier lbhlp
  lbhlp=lbt(inbc)

  ! Interpolation factor for current inner cell value - fxp
  fxp = 1.0_dp-fz(inp)

  ! Interpolation factor for adjecent ghost cell and its reciprocal value
  fxe = 1.0_dp-fxp
  fxer = 1.0_dp / fxe

  if(lbhlp.eq.1) then
  ! inlet

    u(ingc)  = (uin  - u(inp)*fxp) * fxer
    v(ingc)  = (vin  - v(inp)*fxp) * fxer
    w(ingc)  = (win  - w(inp)*fxp) * fxer

    te(ingc) = (tein - te(inp)*fxp) * fxer
    ed(ingc) = (edin - ed(inp)*fxp) * fxer
    t(ingc) = (tin - t(inp)*fxp) * fxer

    vis(ingc) = 0.09_dp*te(ingc)**2/ed(ingc)

  endif

  if(lbhlp.eq.2) then
  ! outlet - extrapolate from two inside cells

    u(ingc) = u(inp)+(u(inp)-u(inp-nij))*(1.0_dp-fz(inp-nij))
    v(ingc) = v(inp)+(v(inp)-v(inp-nij))*(1.0_dp-fz(inp-nij))
    w(ingc) = w(inp)+(w(inp)-w(inp-nij))*(1.0_dp-fz(inp-nij))

    te(ingc) = te(inp)+(te(inp)-te(inp-nij))*(1.0_dp-fz(inp-nij))
    ed(ingc) = ed(inp)+(ed(inp)-ed(inp-nij))*(1.0_dp-fz(inp-nij))
    t(ingc) = t(inp)+(t(inp)-t(inp-nij))*(1.0_dp-fz(inp-nij))

    vis(ingc) = vis(inp)+(vis(inp)-vis(inp-nij))*(1.0_dp-fz(inp-nij))

  endif

  if(lbhlp.eq.3) then
  ! symmetry
  
    Unmag = u(inp)*ar3x(iface)+v(inp)*ar3y(iface)+w(inp)*ar3z(iface)

    u(ingc) = u(inp)-Unmag*ar3x(iface)
    v(ingc) = v(inp)-Unmag*ar3y(iface)
    w(ingc) = w(inp)-Unmag*ar3z(iface)

    te(ingc) = te(inp)
    ed(ingc) = ed(inp)
    t(ingc) = t(inp)
    vis(ingc) = vis(inp)

  endif

  if(lbhlp.eq.4) then
  ! wall 

    ! No-slip for velocity (zero at wall)
    u(ingc)  = ( - u(inp)*fxp) * fxer
    v(ingc)  = ( - v(inp)*fxp) * fxer
    w(ingc)  = ( - w(inp)*fxp) * fxer
    
    ! Zero for turbulence kinetic energy
    te(ingc) = ( - te(inp)*fxp) * fxer

    ed(ingc) = ed(inp) ! ovde kako?? znamo vrednost za eps na zidu samim tim je Dirichlet

    t(ingc) = (tin - t(inp)*fxp) * fxer ! Ako je isothermal wall

    ! Dirichlet for effective viscosity
    ! vis(ingc) = ( vis(ingc) - vis(inp)*fxp ) * fxer
    vis(ingc) = vis(inp)

  endif


  end do
  end do



! TODO sta radimo sa eps, sta se radi za symmetry??? da li se oduzima ta normalna komponenta
! Temperatura, koncentracija, varT,... sve sto se resava preko linearnog sistema.

end subroutine