!***********************************************************************
!
      subroutine calcuvw
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use geometry
      use coef
      use coefb
      use variables
      use buoy
      use time_mod
      use gradients

      implicit none
!
!***********************************************************************

!
!     local variables
!
      integer, parameter :: isol=1
      integer :: i, j, k, inp
      real(prec) :: urfrs, urfms, apotime, heat
      real(prec) :: sut, svt, swt

      logical :: ScndOrderWallBC_Model = .false.
      integer :: ijp, inbc, ingc, iface
      real(dp) :: fxp,fxe,ground_zero_x,ground_zero_y,ground_zero_z,dnw,dpb,are,nxf,nyf,nzf
      real(prec) :: utp,vtp,wtp,upb,vpb,wpb,vnp,viss,cf,vsol,fdui,fdvi,fdwi,fdne



!.....velocity gradients: 
      if (lstsq) then 
        call grad_lsq_qr(u,gradu,2)
        call grad_lsq_qr(v,gradv,2)
        call grad_lsq_qr(w,gradw,2)
      elseif (gauss) then
        call grad_gauss(u,gradu(1,:),gradu(2,:),gradu(3,:))
        call grad_gauss(v,gradv(1,:),gradv(2,:),gradv(3,:))
        call grad_gauss(w,gradw(1,:),gradw(2,:),gradw(3,:))
      endif

!.....calculate pressure gradients
      if (lstsq) then   
        call grad_lsq_qr(pp,gradp,2)
      elseif (gauss) then
        call grad_gauss(pp,gradp(1,:),gradp(2,:),gradp(3,:))
      endif

!.....calculate source terms integrated over volume
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j
!
!.....for u  sp => ap; for v  sp => spv; for w  sp => sp
!.....sum source terms
      ap(inp) = 0.0d0
      spv(inp) = 0.0d0
      sp(inp) = 0.0d0

!=======================================================================
!.....pressure source terms
!=======================================================================
!.....treating it as a volume source
      su(inp) = -gradp(1,inp)*vol(inp)
      sv(inp) = -gradp(2,inp)*vol(inp)
      sw(inp) = -gradp(3,inp)*vol(inp)

!.....constant mass flow forcing - used only on u velocity component
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(const_mflux) su(inp) = su(inp) + gradpcmf*vol(inp)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!=======================================================================
!.....buoyancy source terms
!=======================================================================
      if(lcal(ien).and.lbuoy) then
!-----------------------------------------------------------------------
!........[boussinesq-ova aproximacija: ]
!-----------------------------------------------------------------------
      heat=0.0d0
      if(boussinesq) then
        heat=beta*densit*(t(inp)-tref)*vol(inp)
      else !if(boussinesq.eq.0) 
        heat=(densit-den(inp))*vol(inp)
      endif
!-----------------------------------------------------------------------
      su(inp)=su(inp)-gravx*heat
      sv(inp)=sv(inp)-gravy*heat
      sw(inp)=sw(inp)-gravz*heat
      endif

!=======================================================================
!.....unsteady term
!=======================================================================
      if(bdf) then
!-----------------------------------------------------------------------
!    three level implicit time integration method:
!    in case that btime=0. --> implicit euler
!-----------------------------------------------------------------------
      apotime=den(inp)*vol(inp)/timestep
      sut=apotime*((1+btime)*uo(inp)-0.5*btime*uoo(inp))
      svt=apotime*((1+btime)*vo(inp)-0.5*btime*voo(inp))
      swt=apotime*((1+btime)*wo(inp)-0.5*btime*woo(inp))

      su(inp)=su(inp)+sut
      sv(inp)=sv(inp)+svt
      sw(inp)=sw(inp)+swt

      ap(inp)=ap(inp)+apotime*(1+0.5*btime)
      spv(inp)=spv(inp)+apotime*(1+0.5*btime)
      sp(inp)=sp(inp)+apotime*(1+0.5*btime)
!-----------------------------------------------------------------------
      endif

      end do
      end do
      end do

!
!.....calc reynols stresses explicitly
      call calcstress
!
!=======================================================================
!.....[additional terms: ]
!=======================================================================
      if(lturb.and.lasm) then
!-----the additional algebraic stress source terms ---------------------
      call additional_algebraic_stress_terms
!-----------------------------------------------------------------------
      end if  ![asm approach]


!.....fluxes

!.....east cell - face
      call fluxuvw(nj,1,nij, &
                   ar1x,ar1y,ar1z, &
                   fx,ae,aw,f1)

!.....north cell - face
      call fluxuvw(1,nij,nj, &
                   ar2x,ar2y,ar2z, &
                   fy,an,as,f2)

!.....top   cell - face
      call fluxuvw(nij,nj,1, &
                   ar3x,ar3y,ar3z, &
                   fz,at,ab,f3)  


  !----------------------------------------------------------------------------
  ! Wall boundaries
  !----------------------------------------------------------------------------

  ! !
  ! ! ovde moze da se imaju parovi ijp,ijgc za sve granice sacuvati, za svaku granicu zasepno i onda se ici po njima loop
  ! !

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

  ! ! Normal distance from wall of the cell center of wall adjecent cell
  ! dnw = abs( (xc(inbc)-ground_zero_x)*nx + (yc(inbc)-ground_zero_y)*ny + (zc(inbc)-ground_zero_z)*nz )

  ! viss = vis(ingc) ! Viscosity foor wall  calculated in calc_vis_les
  ! cf=viss*are/dnw ! cf v.s. vsol -> cf is calculated using normal distance!

  ! ! Dist from cc of owner cell to cf at wall boundary, cannot expect dpb to be normal to boundary face in general
  ! dpb = 0.5_dp * sqrt( (xc(inbc)-xc(ingc))**2 + (yc(inbc)-yc(ingc))**2 + (zc(inbc)-zc(ingc))**2 )

  ! ! Diffusion coef. 
  ! vsol = viss*are/dpb

  ! ! Velocity difference vector components
  ! upb = u(inbc)!-u(ijb) ! Change for moving wall u(ijb) \= 0
  ! vpb = v(inbc)!-v(ijb)
  ! wpb = w(inbc)!-w(ijb)

  ! ! Velocity difference vector projected to wall face normal.
  ! vnp = upb*nxf+vpb*nyf+wpb*nzf

  ! ! Velocity difference in tangential direction.
  ! utp = upb-vnp*nxf
  ! vtp = vpb-vnp*nyf
  ! wtp = wpb-vnp*nzf


  ! if (ScndOrderWallBC_Model) then

  !   ! Eksplicitna difuzija
  !   FdUi=viss*((gradu(1,ijp)+gradu(1,ijp))*nxf+(gradu(2,ijp)+gradv(1,ijp))*nyf+(gradu(3,ijp)+gradw(1,ijp))*nzf)
  !   FdVi=viss*((gradv(1,ijp)+gradu(2,ijp))*nxf+(gradv(2,ijp)+gradv(2,ijp))*nyf+(gradv(3,ijp)+gradw(2,ijp))*nzf)
  !   FdWi=viss*((gradw(1,ijp)+gradu(3,ijp))*nxf+(gradw(2,ijp)+gradv(3,ijp))*nyf+(gradw(3,ijp)+gradw(3,ijp))*nzf)
  !   ! Projektujes eksplicitnu difuziju na nomalu
  !   FdNe = FdUi*nxf + FdVi*nyf + FdWi*nzf
  !   ! oduzmes od eksplicitne difuzije onu komponentu duz normale
  !   FdUi = FdUi-FdNe*nxf
  !   FdVi = FdVi-FdNe*nyf
  !   FdWi = FdWi-FdNe*nzf

  !   ap(inbc)  = ap(inbc) + vsol
  !   spv(inbc) = spv(inbc) + vsol
  !   sp(inbc)  = sp(inbc)  + vsol

  !   su(inbc) = su(inbc)+Vsol*U(inbc)-(2*cf*Utp+FdUi*Are)
  !   sv(inbc) = sv(inbc)+Vsol*V(inbc)-(2*cf*Vtp+FdVi*Are)
  !   sw(inbc) = sw(inbc)+Vsol*W(inbc)-(2*cf*Wtp+FdWi*Are)

  ! else

  !   ap(inbc)  = ap(inbc) + vsol
  !   spv(inbc) = spv(inbc) + vsol
  !   sp(inbc)  = sp(inbc)  + vsol

  !   su(inbc) = su(inbc) + vsol*u(inbc) - cf*utp
  !   sv(inbc) = sv(inbc) + vsol*v(inbc) - cf*vtp
  !   sw(inbc) = sw(inbc) + vsol*w(inbc) - cf*wtp

  ! endif


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

  ! ! Normal distance from wall of the cell center of wall adjecent cell
  ! dnw = abs( (xc(inbc)-ground_zero_x)*nx + (yc(inbc)-ground_zero_y)*ny + (zc(inbc)-ground_zero_z)*nz )

  ! viss = vis(ingc) ! Viscosity foor wall  calculated in calc_vis_les
  ! cf=viss*are/dnw ! cf v.s. vsol -> cf is calculated using normal distance!

  ! ! Dist from cc of owner cell to cf at wall boundary, cannot expect dpb to be normal to boundary face in general
  ! dpb = 0.5_dp * sqrt( (xc(inbc)-xc(ingc))**2 + (yc(inbc)-yc(ingc))**2 + (zc(inbc)-zc(ingc))**2 )

  ! ! Diffusion coef. 
  ! vsol = viss*are/dpb

  ! ! Velocity difference vector components
  ! upb = u(inbc)!-u(ijb)  ! Change for moving wall u(ijb) \= 0
  ! vpb = v(inbc)!-v(ijb)
  ! wpb = w(inbc)!-w(ijb)

  ! ! Velocity difference vector projected to wall face normal.
  ! vnp = upb*nxf+vpb*nyf+wpb*nzf

  ! ! Velocity difference in tangential direction.
  ! utp = upb-vnp*nxf
  ! vtp = vpb-vnp*nyf
  ! wtp = wpb-vnp*nzf


  ! if (ScndOrderWallBC_Model) then

  !   ! Eksplicitna difuzija
  !   FdUi=viss*((gradu(1,ijp)+gradu(1,ijp))*nxf+(gradu(2,ijp)+gradv(1,ijp))*nyf+(gradu(3,ijp)+gradw(1,ijp))*nzf)
  !   FdVi=viss*((gradv(1,ijp)+gradu(2,ijp))*nxf+(gradv(2,ijp)+gradv(2,ijp))*nyf+(gradv(3,ijp)+gradw(2,ijp))*nzf)
  !   FdWi=viss*((gradw(1,ijp)+gradu(3,ijp))*nxf+(gradw(2,ijp)+gradv(3,ijp))*nyf+(gradw(3,ijp)+gradw(3,ijp))*nzf)
  !   ! Projektujes eksplicitnu difuziju na nomalu
  !   FdNe = FdUi*nxf + FdVi*nyf + FdWi*nzf
  !   ! oduzmes od eksplicitne difuzije onu komponentu duz normale
  !   FdUi = FdUi-FdNe*nxf
  !   FdVi = FdVi-FdNe*nyf
  !   FdWi = FdWi-FdNe*nzf

  !   ap(inbc)  = ap(inbc) + vsol
  !   spv(inbc) = spv(inbc) + vsol
  !   sp(inbc)  = sp(inbc)  + vsol

  !   su(inbc) = su(inbc)+Vsol*U(inbc)-(2*cf*Utp+FdUi*Are)
  !   sv(inbc) = sv(inbc)+Vsol*V(inbc)-(2*cf*Vtp+FdVi*Are)
  !   sw(inbc) = sw(inbc)+Vsol*W(inbc)-(2*cf*Wtp+FdWi*Are)

  ! else

  !   ap(inbc)  = ap(inbc) + vsol
  !   spv(inbc) = spv(inbc) + vsol
  !   sp(inbc)  = sp(inbc)  + vsol

  !   su(inbc) = su(inbc) + vsol*u(inbc) - cf*utp
  !   sv(inbc) = sv(inbc) + vsol*v(inbc) - cf*vtp
  !   sw(inbc) = sw(inbc) + vsol*w(inbc) - cf*wtp

  ! endif


  ! enddo
  ! enddo



!
!.....Assemble and solve system
!

      if(cn) then
!.....crank-nicolson stuff - only once:
      ae = 0.5d0*ae
      an = 0.5d0*an
      at = 0.5d0*at
      aw = 0.5d0*aw
      as = 0.5d0*as
      ab = 0.5d0*ab
      endif

      urfrs=urfr(iu)
      urfms=urfm(iu)

      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      if(cn) then
!.....crank-nicolson time stepping source terms
      apotime=den(inp)*vol(inp)/timestep
      su(inp)=su(inp)+(ae(inp)*uo(inp+nj)  + aw(inp)*uo(inp-nj)+    &
                       an(inp)*uo(inp+1)   + as(inp)*uo(inp-1)+     &
                       at(inp)*uo(inp+nij) + ab(inp)*uo(inp-nij))+  &
              (apotime-ae(inp)-aw(inp)                              &
                      -an(inp)-as(inp)                              &
                      -at(inp)-ab(inp))*uo(inp)               
      ap(inp)=ap(inp)+apotime
!.....end of crank-nicolson time stepping source terms
!.....end of crank-nicolson stuff
      endif

      ap(inp)=ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp)+ap(inp)
      ap(inp)=ap(inp)*urfrs
      su(inp)=su(inp)+urfms*ap(inp)*u(inp)
      apu(inp)=1./(ap(inp)+small) ! simple,piso
      !apu(inp)=1./(ap(inp)-(ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp))+small) ! simplec
      enddo
      enddo
      enddo
!
!.....solving f.d. equations
      if(isol.eq.1) call sipsol(u,iu)
      if(isol.eq.2) call cgstab_sip(u,iu) 
!
      urfrs=urfr(iv)
      urfms=urfm(iv)

      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      if(cn) then
!.....crank-nicolson time stepping source terms
      apotime=den(inp)*vol(inp)/timestep
      sv(inp)=sv(inp)+(ae(inp)*vo(inp+nj)  + aw(inp)*vo(inp-nj)+    &
                       an(inp)*vo(inp+1)   + as(inp)*vo(inp-1)+     &
                       at(inp)*vo(inp+nij) + ab(inp)*vo(inp-nij))+  &
              (apotime-ae(inp)-aw(inp)                              &
                      -an(inp)-as(inp)                              &
                      -at(inp)-ab(inp))*vo(inp)
      spv(inp)=spv(inp)+apotime
!.....end of crank-nicolson time stepping source terms
      endif

      ap(inp)=ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp)+spv(inp)
      ap(inp)=ap(inp)*urfrs
      su(inp)=sv(inp)+urfms*ap(inp)*v(inp)
      apv(inp)=1./(ap(inp)+small) ! simple,piso
      !apv(inp)=1./(ap(inp)-(ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp))+small) ! simplec
      enddo
      enddo
      enddo
!
!.....solving f.d. equations
!
      if(isol.eq.1) call sipsol(v,iv)
      if(isol.eq.2) call cgstab_sip(v,iv)
!
!.....problem modifications - boundary conditions for  w - comp.
      urfrs=urfr(iw)
      urfms=urfm(iw)

      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      if(cn) then
!.....crank-nicolson time stepping source terms
      apotime=den(inp)*vol(inp)/timestep
      sw(inp)=sw(inp)+(ae(inp)*wo(inp+nj)  + aw(inp)*wo(inp-nj)+    &
                       an(inp)*wo(inp+1)   + as(inp)*wo(inp-1)+     &
                       at(inp)*wo(inp+nij) + ab(inp)*wo(inp-nij))+  &
              (apotime-ae(inp)-aw(inp)                              &
                      -an(inp)-as(inp)                              &
                      -at(inp)-ab(inp))*wo(inp)
      sp(inp)=sp(inp)+apotime
!.....end of crank-nicolson time stepping source terms
      endif

      ap(inp)=ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp)+sp(inp)
      ap(inp)=ap(inp)*urfrs
      su(inp)=sw(inp)+urfms*ap(inp)*w(inp)
      apw(inp)=1./(ap(inp)+small) ! simple,piso
      !apw(inp)=1./(ap(inp)-(ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp))+small) ! simplec
      enddo
      enddo
      enddo
!
!.....solving f.d. equations
      if(isol.eq.1) call sipsol(w,iw)
      if(isol.eq.2) call cgstab_sip(w,iw)

      return
      end
