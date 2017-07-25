!***********************************************************************
!
  subroutine calc_vis_les
!
!***********************************************************************
!
!     Update effective viscosity = molecular dynamic visc. + eddy visc.
!     We need fresh vel. gradients for some models - we got them in calcscm
!     when they are evaluated after vel. correction in calcp - these 
!     are the freshest they can be in this outer (SIMPLE) interation.
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use coef
  use variables
  use boundc

  implicit none
!
!***********************************************************************
!
  integer :: i, j, k, inp
  real(prec) :: visold
  real(prec) :: wldist
  real(prec) :: Uref,vis_smag,c_les,zplus
  real(prec), parameter :: C_delta = 0.158

  integer :: inbc, ingc, iface
  real(dp) :: fxp,fxe,ground_zero_x,ground_zero_y,ground_zero_z,dnw,are,nxf,nyf,nzf
  real(dp) :: vnp,xtp,ytp,ztp,vtp,ut2,wall_shear_stress,viscw



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! if(Smagor_sgs) then
!====================================================
!.....FOR Smagorinsky with Van Driest damping
!====================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  do k=3,nkmm
  do i=3,nimm
  do j=3,njmm

  inp=lk(k)+li(i)+j

  visold = vis(inp)

  C_Les = 0.1 ! 0.1-0.2; for channel flow even 0.065(Ferziger)   

  WlDist = wallDistance(inp)

  Uref = sqrt(U(Inp)*U(Inp)+V(Inp)*V(Inp)+W(Inp)*W(Inp)) 

  Utau = sqrt(Viscos*Uref/WlDist)

  Zplus = WlDist*Utau/Viscos

  !C_Les = C_Les*(1.-Exp(-Zplus/26.0))
  !Vis_Smag = (C_Les * Vol(Inp)**(1./3.))**2 * Strain(inp) * Den(Inp)
  !Vis(inp) = Vis_Smag + Viscos

  ! Verzija iz OpenFOAM-malo modifikovana, spojene ideje Schumann(1991) i Van Driest:
  C_Les = 0.17
  Vis_Smag = ( min(C_Les*Vol(Inp)**(1./3.),cappa/C_delta*Wldist) * (1.-Exp(-Zplus/26.0)) )**2 * Strain(inp) * Den(Inp)
  Vis(inp) = Vis_Smag + Viscos


  vis(inp)=urf(ivis)*vis(inp)+(1.-urf(ivis))*visold

  enddo
  enddo
  enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! endif


!----------------------------------------------------------------------------
! Wall boundaries - update Visw and Ypl
!----------------------------------------------------------------------------

!
! ovde moze da se parovi ijp,ijgc za sve granice sacuvaju, za svaku granicu zasepno i onda se ici po njima loop
!

  ! ! Bottom
  ! do i=3,nimm
  ! do j=3,njmm

  ! inbc=lk(3)+li(i)+j 
  ! ingc = inbc - nij

  ! iface  = ingc

  ! fxp = fz(ingc)
  ! fxe = 1.0_dp-fxp

  ! ! At what height is the wall
  ! ground_zero_x = xc(inbc)*fxp+xc(ingc)*fxe
  ! ground_zero_y = yc(inbc)*fxp+yc(ingc)*fxe
  ! ground_zero_z = zc(inbc)*fxp+zc(ingc)*fxe

  ! are = sqrt(ar3x(iface)**2+ar3y(iface)**2+ar3z(iface)**2)
  ! ! write(*,*) 'Afce ares components: ',ar3x(iface),ar3y(iface),ar3z(iface)

  ! ! Face normals
  ! nxf = ar3x(iface)/are
  ! nyf = ar3y(iface)/are
  ! nzf = ar3z(iface)/are
  ! ! write(*,*) 'Normals: ', nxf,nyf,nzf

  ! ! Magnitude of a cell center velocity projected on boundary face normal
  ! Vnp = U(inbc)*nxf+V(inbc)*nyf+W(inbc)*nzf
  ! ! write(*,*) 'Normal velocity magnitude Vnp = ',vnp

  ! ! Tangential velocity components 
  ! xtp = U(inbc)-Vnp*nxf
  ! ytp = V(inbc)-Vnp*nyf
  ! ztp = W(inbc)-Vnp*nzf
  ! ! write(*,*) 'Tangential velocity components = ',xtp,ytp,ztp

  ! ! Its magnitude
  ! Vtp = sqrt(xtp*xtp+ytp*ytp+ztp*ztp)
  ! ! write(*,*) 'Tangential velocity magnitude Vtp = ',vtp

  ! ! Tangent direction
  ! xtp = xtp/vtp
  ! ytp = ytp/vtp
  ! ztp = ztp/vtp
  ! ! write(*,*) 'Tangential direction components: ',xtp,ytp,ztp

  ! ! projektovanje razlike brzina na pravac tangente u cell centru ijp
  ! Ut2 = abs( U(inbc)*xtp + V(inbc)*ytp + W(inbc)*ztp) 
  ! ! write(*,*) 'ut2 = ', ut2

  ! ! Normal distance from wall of the cell center of wall adjecent cell
  ! dnw = abs( (xc(inbc)-ground_zero_x)*nx + (yc(inbc)-ground_zero_y)*ny + (zc(inbc)-ground_zero_z)*nz )
  ! ! write(*,*) 'Normal distance of the cell center - dnw = ',dnw

  ! wall_shear_stress = viscos*Ut2/dnw
  ! ! write(*,*) 'Tau = ',Tau

  ! ! u_tau iz mat. def. u_tau=sqrt(Tau_w/rho):                                                                   
  ! utau=sqrt(wall_shear_stress/densit)

  ! ! Ypl = rho* \Delta y * u_tau / \mu
  ! ypl = densit*utau*dnw/viscos
  ! ! write(*,*) 'ypl = ',yplus(ijk)

  !   ! ! Ima i ova varijanta u cisto turb. granicni sloj varijanti sa prvom celijom u log sloju
  !   ! ypl(i) = den(ijb)*cmu25*sqrt(te(ijp))*dnw(i)/viscos
  !   ! ! ...ovo je tehnicki receno ystar iliti y* a ne y+

  ! viscw = zero

  ! if(ypl > 11.61) then
  !   viscw = ypl*viscos*cappa/log(Elog*ypl)
  ! endif

  ! vis(ingc) = max(viscos,viscw)

  ! enddo
  ! enddo


  ! ! Top
  ! do i=3,nimm
  ! do j=3,njmm

  ! inbc=lk(nimm)+li(i)+j 
  ! ingc = inbc + nij

  ! iface  = inbc

  ! fxp = 1.0_dp-fz(inbc)
  ! fxe = 1.0_dp-fxp

  ! ! At what height is the wall
  ! ground_zero_x = xc(inbc)*fxp+xc(ingc)*fxe
  ! ground_zero_y = yc(inbc)*fxp+yc(ingc)*fxe
  ! ground_zero_z = zc(inbc)*fxp+zc(ingc)*fxe

  ! are = sqrt(ar3x(iface)**2+ar3y(iface)**2+ar3z(iface)**2)
  ! ! write(*,*) 'Afce ares components: ',ar3x(iface),ar3y(iface),ar3z(iface)

  ! ! Face normals
  ! nxf = ar3x(iface)/are
  ! nyf = ar3y(iface)/are
  ! nzf = ar3z(iface)/are
  ! ! write(*,*) 'Normals: ', nxf,nyf,nzf

  ! ! Magnitude of a cell center velocity projected on boundary face normal
  ! Vnp = U(inbc)*nxf+V(inbc)*nyf+W(inbc)*nzf
  ! ! write(*,*) 'Normal velocity magnitude Vnp = ',vnp

  ! ! Tangential velocity components 
  ! xtp = U(inbc)-Vnp*nxf
  ! ytp = V(inbc)-Vnp*nyf
  ! ztp = W(inbc)-Vnp*nzf
  ! ! write(*,*) 'Tangential velocity components = ',xtp,ytp,ztp

  ! ! Its magnitude
  ! Vtp = sqrt(xtp*xtp+ytp*ytp+ztp*ztp)
  ! ! write(*,*) 'Tangential velocity magnitude Vtp = ',vtp

  ! ! Tangent direction
  ! xtp = xtp/vtp
  ! ytp = ytp/vtp
  ! ztp = ztp/vtp
  ! ! write(*,*) 'Tangential direction components: ',xtp,ytp,ztp

  ! ! projektovanje razlike brzina na pravac tangente u cell centru ijp
  ! Ut2 = abs( U(inbc)*xtp + V(inbc)*ytp + W(inbc)*ztp) 
  ! ! write(*,*) 'ut2 = ', ut2

  ! ! Normal distance from wall of the cell center of wall adjecent cell
  ! dnw = abs( (xc(inbc)-ground_zero_x)*nx + (yc(inbc)-ground_zero_y)*ny + (zc(inbc)-ground_zero_z)*nz )
  ! ! write(*,*) 'Normal distance of the cell center - dnw = ',dnw

  ! wall_shear_stress = viscos*Ut2/dnw
  ! ! write(*,*) 'Tau = ',wall_shear_stress

  ! ! u_tau iz mat. def. u_tau=sqrt(Tau_w/rho):                                                                   
  ! utau=sqrt(wall_shear_stress/densit)

  ! ! Ypl = rho* \Delta y * u_tau / \mu
  ! ypl = densit*utau*dnw/viscos
  ! ! write(*,*) 'ypl = ',ypl

  !   ! ! Ima i ova varijanta u cisto turb. granicni sloj varijanti sa prvom celijom u log sloju
  !   ! ypl(i) = den(ijb)*cmu25*sqrt(te(ijp))*dnw(i)/viscos
  !   ! ! ...ovo je tehnicki receno ystar iliti y* a ne y+

  ! viscw = zero

  ! if(ypl > 11.61) then
  !   viscw = ypl*viscos*cappa/log(Elog*ypl)
  ! endif

  ! vis(ingc) = max(viscos,viscw)

  ! enddo
  ! enddo

  end subroutine
