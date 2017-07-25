!***********************************************************************
!
      subroutine fluxuvw_explicit(idew,idns,idtb, &
                                  arvx,arvy,arvz, &
                                  fif,fcf)
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
      use fieldmanipulation

      implicit none
!
!***********************************************************************
! 
      integer, intent(in) :: idew,idns,idtb 
      real(prec), dimension(nxyza) :: arvx,arvy,arvz
      real(prec), dimension(nxyza) :: fif,fcf 
!
!     local variables
!
      integer :: i, j, k, ijk, inp,       &
                 ine,inw, ins,inb,inbs,inee
      integer :: dir

      !real(prec) :: dxu,dyu,dzu,dxu1,dyu1,dzu1
      !real(prec) :: uinu,vinu,winu,uinu1,vinu1,winu1
      !real(prec) :: dfidx_f,dfidy_f,dfidz_f
      !real(prec) :: dxet,dxzd,dyet,dyzd,dzet,dzzd

      real(prec) :: fxw,fxpw,fxee,fxew, &
                    dxc,dyc,dzc,dxe,dye,dze
      real(prec) :: dfidx_u,dfidx_u1, &
                    dfidy_u,dfidy_u1, &
                    dfidz_u,dfidz_u1
      real(prec) :: fxe,fxp,are,vole,game,&
                    ue,ve,we, & 
                    flcf,ce,cp, &
                    g11, g12, g21, g22
      real(prec) :: ae1,aw1,shigh1,shigh2,shigh3
      real(prec) :: r1,r2,r3,r4,r5,r6,                   &
                    psie1,psie2,psie3,psiw1,psiw2,psiw3, &
                    fuhigh,fvhigh,fwhigh
      real(prec) :: xf,yf,zf,xi,yi,zi,arx,ary,arz, &
                    duxi,duyi,duzi, &
                    dvxi,dvyi,dvzi, &
                    dwxi,dwyi,dwzi
      real(prec) :: xpn,ypn,zpn
      real(prec) :: nxx,nyy,nzz,ixi1,ixi2,ixi3,dpn,costheta,costn
      real(prec) :: fdue,fdve,fdwe
      real(prec) :: onethird, twothirds
      real(prec) :: d2x,d2y,d2z,d1x,d1y,d1z
      real(prec) :: duxii,dvxii,dwxii, &
                    duyii,dvyii,dwyii, &
                    duzii,dvzii,dwzii
!----------------------------------------------------------------------

      onethird = 1.0d0/3.0d0
      twothirds = 2.0d0/3.0d0

!.....Initialize values for high order conv. fluxes
      psie1=0.0d0
      psie2=0.0d0
      psie3=0.0d0
      psiw1=0.0d0
      psiw2=0.0d0
      psiw3=0.0d0 

      fuhigh=0.0d0
      fvhigh=0.0d0
      fwhigh=0.0d0

!.....Choose direction for gradients; if dir=1 we use dfidx or grad(inp,1), 2=>dfidy, 3=>dfidz
      if(idew.eq.nj)  dir = 1
      if(idew.eq.1)   dir = 2
      if(idew.eq.nij) dir = 3
!
!.....Calculate east,top,north  cell face
      do k=2,nkmm
      do i=2,nimm
      do j=2,njmm

      inp=lk(k)+li(i)+j

      ine=inp+idew
      inw=inp-idew
      ins=inp-idns
      inb=inp-idtb
      inbs=inb-idns

      inee=inp+2*idew
      inbs=inb-idns
!
!.....Interpolation factors in first,second and third direction
      fxe=fif(inp) 
      fxp=1.0d0-fxe

      fxpw = fif(inw)
      fxw = 1.0d0-fxpw
      fxee = fif(ine)
      fxew = 1.0d0-fxee 


!.....Distance from P to neighbor N
      xpn=xc(ine)-xc(inp)
      ypn=yc(ine)-yc(inp)
      zpn=zc(ine)-zc(inp)

      dpn=sqrt(xpn**2+ypn**2+zpn**2)     

!.....Components of the unit vector i_ksi
      ixi1=xpn/dpn
      ixi2=ypn/dpn
      ixi3=zpn/dpn


!.....Precomputed face areas
      arx=arvx(inp)
      ary=arvy(inp)
      arz=arvz(inp)

!.....Cell face area
      are=sqrt(arx**2+ary**2+arz**2)

!.....Unit vectors of the normal
      nxx=arx/are
      nyy=ary/are
      nzz=arz/are

!.....Angle between vectorsa n and i_xi - we need cosine
      costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

!.....Relaxation factor for higher-order cell face gradient
      ! Minimal correction: nrelax = +1 :
      !costn = costheta
      ! Orthogonal correction: nrelax =  0 : 
      costn = 1.0d0
      ! Over-relaxed approach: nrelax = -1 :
      !costn = 1./costheta
      ! In general, nrelax can be any signed integer from some 
      ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
      !costn = costheta**nrelax

      vole=xpn*arx+ypn*ary+zpn*arz

!.....Overrelaxed correction vector d2, where S=dpn+d2
      d1x = costn
      d1y = costn
      d1z = costn
      !
      d2x = xpn*costn
      d2y = ypn*costn
      d2z = zpn*costn

!.....Cell face viscosity
      game=vis(inp)*fxp+vis(ine)*fxe

!++++VELOCITIES AT CELL FACE CENTER and EXPLICIT DIFFUSION FLUXES+++++++

!.....Coordinates of cell-face center
      xf = 0.25d0*(x(inp)+x(ins)+x(inb)+x(inbs))
      yf = 0.25d0*(y(inp)+y(ins)+y(inb)+y(inbs))
      zf = 0.25d0*(z(inp)+z(ins)+z(inb)+z(inbs))

!.....Coordinates of point e'
      xi=xc(inp)*fxp+xc(ine)*fxe
      yi=yc(inp)*fxp+yc(ine)*fxe
      zi=zc(inp)*fxp+zc(ine)*fxe

!>>> U
!.....Interpolate gradients defined at CV centers to faces
      duxi = gradu(1,inp)*fxp+gradu(1,ine)*fxe
      duyi = gradu(2,inp)*fxp+gradu(2,ine)*fxe
      duzi = gradu(3,inp)*fxp+gradu(3,ine)*fxe

!.....du/dx_i interpolated at cell face:
      duxii = duxi*d1x + arx/vole*( u(ine)-u(inp)-duxi*d2x-duyi*d2y-duzi*d2z ) 
      duyii = duyi*d1y + ary/vole*( u(ine)-u(inp)-duxi*d2x-duyi*d2y-duzi*d2z ) 
      duzii = duzi*d1z + arz/vole*( u(ine)-u(inp)-duxi*d2x-duyi*d2y-duzi*d2z ) 

!>>> V
!.....Interpolate gradients defined at CV centers to faces
      dvxi = gradv(1,inp)*fxp+gradv(1,ine)*fxe
      dvyi = gradv(2,inp)*fxp+gradv(2,ine)*fxe
      dvzi = gradv(3,inp)*fxp+gradv(3,ine)*fxe

!.....dv/dx_i interpolated at cell face:
      dvxii = dvxi*d1x + arx/vole*( v(ine)-v(inp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 
      dvyii = dvyi*d1y + ary/vole*( v(ine)-v(inp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 
      dvzii = dvzi*d1z + arz/vole*( v(ine)-v(inp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 

!>>> W
!.....Interpolate gradients defined at CV centers to faces
      dwxi = gradw(1,inp)*fxp+gradw(1,ine)*fxe
      dwyi = gradw(2,inp)*fxp+gradw(2,ine)*fxe
      dwzi = gradw(3,inp)*fxp+gradw(3,ine)*fxe

!.....dw/dx_i interpolated at cell face:
      dwxii = dwxi*d1x + arx/vole*( w(ine)-w(inp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 
      dwyii = dwyi*d1y + ary/vole*( w(ine)-w(inp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 
      dwzii = dwzi*d1z + arz/vole*( w(ine)-w(inp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! EXPLICIT DIFFUSSION 
       fdue = game*( (dUxii+dUxii)*arx + (dUyii+dVxii)*ary + (dUzii+dWxii)*arz )
       fdve = game*( (dUyii+dVxii)*arx + (dVyii+dVyii)*ary + (dVzii+dWyii)*arz )
       fdwe = game*( (dUzii+dWxii)*arx + (dWyii+dVzii)*ary + (dWzii+dWzii)*arz )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++END: VELOCITIES AT CELL FACE CENTER and EXPLICIT DIFFUSION FLUXES+++++++

!.....CONVECTION FLUXES - UDS
!     If flow goes P=>E CP=FLCF, CE=0.
!     If flow goes E=>P CE=FLCF, CP=0.

      flcf=fcf(inp)

      ce=min(flcf,zero) 
      cp=max(flcf,zero)

!=====START CDS SCHEME===============================
      IF(LCDS) THEN
!
!.....Explicit convective fluxes for CDS
!

!        |________Ue'_________|_______________Ucorr___________________|
      ue=u(inp)*fxp+u(ine)*fxe!+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))

!        |________Ve'_________|_______________Vcorr___________________|
      ve=v(inp)*fxp+v(ine)*fxe!+(dvxi*(xf-xi)+dvyi*(yf-yi)+dvzi*(zf-zi))

!        |________We'_________|_______________Wcorr___________________|
      we=w(inp)*fxp+w(ine)*fxe!+(dwxi*(xf-xi)+dwyi*(yf-yi)+dwzi*(zf-zi))

!      UE = face_interpolated(U,gradU,inp,idew,idns,idtb,fxp,fxe)
!      VE = face_interpolated(V,gradV,inp,idew,idns,idtb,fxp,fxe)
!      WE = face_interpolated(W,gradW,inp,idew,idns,idtb,fxp,fxe)

      fuhigh=flcf*ue
      fvhigh=flcf*ve
      fwhigh=flcf*we


!=====END CDS SCHEME================================
      END IF

!=====START LUDS SCHEME=============================
      IF(LLUDS) THEN 

!.....EXPLICIT CONVECTIVE FLUXES FOR LUDS SCHEME
!     $Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_e[kg/s] * Phi_e[m/s]$
!     Phi_e is found by extrapolation from upwind nodes, see eq. (3.29) in Sasa's Thesis.

!.....DISTANCE FROM CONSIDERED NODE TO THE FACE
!     In the Waterson & Deconinck article there is a term 'DeltaXc/2', whole this term we denote DXC
!.....If flow goes from P to E
      dxc = xf - xc(inp)
      dyc = yf - yc(inp)
      dzc = zf - zc(inp)
!.....If flow goes from E to P
      dxe = xc(ine) - xf
      dye = yc(ine) - yf
      dze = zc(ine) - zf

!.....Interpolate gradients defined at CV centers to faces
      dfidx_u = gradu(dir,inw)*fxw+gradu(dir,inp)*fxpw
      dfidy_u = gradv(dir,inw)*fxw+gradv(dir,inp)*fxpw
      dfidz_u = gradw(dir,inw)*fxw+gradw(dir,inp)*fxpw
!
      dfidx_u1 = gradu(dir,inee)*fxee+gradu(dir,ine)*fxew
      dfidy_u1 = gradv(dir,inee)*fxee+gradv(dir,ine)*fxew
      dfidz_u1 = gradw(dir,inee)*fxee+gradw(dir,ine)*fxew

      fuhigh = ce*(u(ine) + dxe*dfidx_u1)+ &
               cp*(u(inp) + dxc*dfidx_u)
!       Mass Flux| Interpolation Of Velocity To Face |

      fvhigh = ce*(v(ine) + dye*dfidy_u1)+ &
               cp*(v(inp) + dyc*dfidy_u)
       
      fwhigh = ce*(w(ine) + dze*dfidz_u1)+ &
               cp*(w(inp) + dzc*dfidz_u)

!=====END LUDS SCHEME====================================
      END IF

!=====START QUICK SCHEME=================================
      IF(LQUDS) THEN

      shigh1=0.0d0
      shigh2=0.0d0
      shigh3=0.0d0
      aw1=0.0d0
      ae1=0.0d0

      ! IF(K.GT.2.AND.K.LT.NK-1.AND. &  ! Ukloni ovaj uslov za periodic mrezu!
      !    J.GT.2.AND.J.LT.NJ-1.AND. &
      !    I.GT.2.AND.I.LT.NI-1) THEN

!.....Find coefficients associated with quadratic function
      g11=((2-fif(inp-idew))*fif(inp)**2)/  &
           (1+fif(inp)-fif(inp-idew))
      g12=((1-fif(inp))*(1-fif(inp-idew))**2)/  &
           (1+fif(inp)-fif(inp-idew))
      g21=((1+fif(inp+idew))*(1-fif(inp))**2)/  &
           (1+fif(inp+idew)-fif(inp))
      g22=(fif(inp+idew)**2*fif(inp))/ &
           (1+fif(inp+idew)-fif(inp))

!.....EXPLICIT CONVECTIVE FLUXES FOR QUICK SCHEME
!     $Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_e[kg/s] * Phi_e[m/s]$
!     Phi_e is found by extrapolation from upwind nodes, see eq. (3.30) in Sasa's Thesis.

      fuhigh=ce*(g21*u(inp)+(1-g21+g22)*u(ine)-g22*u(inee))+ &
             cp*(g11*u(ine)+(1-g11+g12)*u(inp)-g12*u(inw))


      fvhigh=ce*(g21*v(inp)+(1-g21+g22)*v(ine)-g22*v(inee))+ &
             cp*(g11*v(ine)+(1-g11+g12)*v(inp)-g12*v(inw))
       
      fwhigh=ce*(g21*w(inp)+(1-g21+g22)*w(ine)-g22*w(inee))+ &
             cp*(g11*w(ine)+(1-g11+g12)*w(inp)-g12*w(inw))

      ! END IF

      END IF
!=============END QUICK SCHEME========================
!--------------------------------------------------------------------------------------------
!     BOUNDED HIGH-ORDER CONVECTIVE SCHEMES (Waterson & Deconinck JCP 224 (2007) pp. 182-207)
      if(flux_limiter) then

!+++++VERSION 1++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++Find r's. This is universal for all schemes. NON-ORTHOGONAL GRID! +++++++
!.....Distances of fictitious upwind point U from 'real' upwind cell center.
!.....If flow goes from P to E, distance from the new node U from closest cell-center
!#      DXU = (XC(INP) - (XC(INE)-XC(INP))) - XC(INW)
!#      DYU = (YC(INP) - (YC(INE)-YC(INP))) - YC(INW)
!#      DZU = (ZC(INP) - (ZC(INE)-ZC(INP))) - ZC(INW)
!.....If flow goes from P to E, distance from the new node U from closest cell-center
!#      DXU1 = (XC(INE) + (XC(INP)-XC(INW))) - XC(INEE)
!#      DYU1 = (YC(INE) + (YC(INP)-YC(INW))) - YC(INEE)
!#      DZU1 = (ZC(INE) + (ZC(INP)-ZC(INW))) - ZC(INEE)
!/////Correction for triple periodicity://///////////
!@      if(i.eq.nie.and.idew.eq.nj)  DXU=0.0D0
!@      if(j.eq.nje.and.idew.eq.1)   DYU=0.0D0
!@      if(k.eq.nke.and.idew.eq.nij) DZU=0.0D0
!@      if(i.eq.2.and.idew.eq.nj)  DXU1=0.0D0
!@      if(j.eq.2.and.idew.eq.1)   DYU1=0.0D0
!@      if(k.eq.2.and.idew.eq.nij) DZU1=0.0D0
!/////END: Correction for triple periodicity://///////////

!.....Variable value at fictitious upwind point
!#      Uinu=U(INW)+DXU*GRADU(1,INW)+DYU*GRADU(2,INW)+DZU*GRADU(3,INW)
!#      Vinu=V(INW)+DXU*GRADV(1,INW)+DYU*GRADV(2,INW)+DZU*GRADV(3,INW)
!#      Winu=W(INW)+DXU*GRADW(1,INW)+DYU*GRADW(2,INW)+DZU*GRADW(3,INW)
!
!#      Uinu1=U(INEE)+DXU1*GRADU(1,INEE)+DYU1*GRADU(2,INEE)+DZU1*GRADU(3,INEE)
!#      Vinu1=V(INEE)+DXU1*GRADV(1,INEE)+DYU1*GRADV(2,INEE)+DZU1*GRADV(3,INEE)
!#      Winu1=W(INEE)+DXU1*GRADW(1,INEE)+DYU1*GRADW(2,INEE)+DZU1*GRADW(3,INEE)

!.....Gradient ratio using value at fictitious upwind point:
!.....If flow goes from P to E
!#      r1 = (U(INE)-U(INP))/(U(INP)-Uinu)
!#      r2 = (V(INE)-V(INP))/(V(INP)-Vinu)
!#      r3 = (W(INE)-W(INP))/(W(INP)-Winu)
!.....If flow goes from E to P
!#      r4 = (U(INP)-U(INE))/(U(INE)-Uinu1)
!#      r5 = (V(INP)-V(INE))/(V(INE)-Vinu1)
!#      r6 = (W(INP)-W(INE))/(W(INE)-Winu1)
!+++++VERSION 1++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++VERSION 2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!.....DISTANCE FROM CONSIDERED NODE TO THE FACE
!     In the Waterson & Deconinck article there is a term 'DeltaXc/2', whole this term we denote DXC
!.....If flow goes from P to E
!      DXC = ( 0.25*(X(INP)+X(INS)+X(INB)+X(INBS)) ) - XC(INP)
!      DYC = ( 0.25*(Y(INP)+Y(INS)+Y(INB)+Y(INBS)) ) - YC(INP)
!      DZC = ( 0.25*(Z(INP)+Z(INS)+Z(INB)+Z(INBS)) ) - ZC(INP)
!.....If flow goes from E to P
!      DXE = XC(INE) - ( 0.25*(X(INP)+X(INS)+X(INB)+X(INBS)) )
!      DYE = YC(INE) - ( 0.25*(Y(INP)+Y(INS)+Y(INB)+Y(INBS)) )
!      DZE = ZC(INE) - ( 0.25*(Z(INP)+Z(INS)+Z(INB)+Z(INBS)) )

!.....Interpolate gradients defined at CV centers to faces
!      DFIDX_f = GRADU(dir,INP)*FXP+GRADU(dir,INE)*FXE 
!      DFIDY_f = GRADV(dir,INP)*FXP+GRADV(dir,INE)*FXE
!      DFIDZ_f = GRADW(dir,INP)*FXP+GRADW(dir,INE)*FXE
!
!      DFIDX_u = GRADU(dir,INW)*FXW+GRADU(dir,INP)*FXPW
!      DFIDY_u = GRADV(dir,INW)*FXW+GRADV(dir,INP)*FXPW
!      DFIDZ_u = GRADW(dir,INW)*FXW+GRADW(dir,INP)*FXPW
!
!      DFIDX_u1 = GRADU(dir,INEE)*FXEE+GRADU(dir,INE)*FXEW
!      DFIDY_u1 = GRADV(dir,INEE)*FXEE+GRADV(dir,INE)*FXEW
!      DFIDZ_u1 = GRADW(dir,INEE)*FXEE+GRADW(dir,INE)*FXEW

!.... Find r's. This is universal for all schemes.
!.....If flow goes from P to E
!      r1 = DFIDX_f/DFIDX_u  
!      r2 = DFIDY_f/DFIDY_u 
!      r3 = DFIDZ_f/DFIDZ_u  
!.....If flow goes from E to P
!      r4 = DFIDX_f/DFIDX_u1  
!      r5 = DFIDY_f/DFIDY_u1 
!      r6 = DFIDZ_f/DFIDZ_u1  
!+++++VERSION 2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++VERSION 3++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!.... Find r's. This is universal for all schemes.
!.....If flow goes from P to E
      r1 = (2*gradU(1,inp)*xpn + 2*gradU(2,inp)*ypn + 2*gradU(3,inp)*zpn)/(U(INE)-U(INP)) - 1.0d0  
      r2 = (2*gradV(1,inp)*xpn + 2*gradV(2,inp)*ypn + 2*gradV(3,inp)*zpn)/(V(INE)-V(INP)) - 1.0d0 
      r3 = (2*gradW(1,inp)*xpn + 2*gradW(2,inp)*ypn + 2*gradW(3,inp)*zpn)/(W(INE)-W(INP)) - 1.0d0 
!.....If flow goes from E to P
      r4 = (2*gradU(1,ine)*xpn + 2*gradU(2,ine)*ypn + 2*gradU(3,ine)*zpn)/(U(INP)-U(INE)) - 1.0d0 
      r5 = (2*gradV(1,ine)*xpn + 2*gradV(2,ine)*ypn + 2*gradV(3,ine)*zpn)/(V(INP)-V(INE)) - 1.0d0 
      r6 = (2*gradW(1,ine)*xpn + 2*gradW(2,ine)*ypn + 2*gradW(3,ine)*zpn)/(W(INP)-W(INE)) - 1.0d0  
!+++++VERSION 3++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!=====SMART SCHEME================================
      IF(LSMART) THEN
!.....PSI for SMART scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(2.*r1, 0.75*r1+0.25, 4.))
      PSIW2 = max(0., min(2.*r2, 0.75*r2+0.25, 4.))
      PSIW3 = max(0., min(2.*r3, 0.75*r3+0.25, 4.))
!.....If flow goes from E to P
      PSIE1 = max(0., min(2.*r4, 0.75*r4+0.25, 4.))
      PSIE2 = max(0., min(2.*r5, 0.75*r5+0.25, 4.))
      PSIE3 = max(0., min(2.*r6, 0.75*r6+0.25, 4.))
!=====END SMART SCHEME=============================
      END IF
!
!=====AVL-SMART SCHEME=============================
      IF(LAVL) THEN
!.....PSI for AVL-SMART scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(1.5*r1, 0.75*r1+0.25, 2.5))
      PSIW2 = max(0., min(1.5*r2, 0.75*r2+0.25, 2.5))
      PSIW3 = max(0., min(1.5*r3, 0.75*r3+0.25, 2.5))
!.....If flow goes from E to P
      PSIE1 = max(0., min(1.5*r4, 0.75*r4+0.25, 2.5))
      PSIE2 = max(0., min(1.5*r5, 0.75*r5+0.25, 2.5))
      PSIE3 = max(0., min(1.5*r6, 0.75*r6+0.25, 2.5))
!=====END AVL-SMART SCHEME==========================
      END IF
!
!=====MUSCL SCHEME=================================
      IF(LMUSCL) THEN
!.....PSI for MUSCL scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
      PSIW2 = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
      PSIW3 = max(0., min(2.*r3, 0.5*r3+0.5, 2.))
!.....If flow goes from E to P
      PSIE1 = max(0., min(2.*r4, 0.5*r4+0.5, 2.))
      PSIE2 = max(0., min(2.*r5, 0.5*r5+0.5, 2.))
      PSIE3 = max(0., min(2.*r6, 0.5*r6+0.5, 2.))

!.....PSI for Koren scheme:
!.....If flow goes from P to E
!      PSIW1 = max(0., min(2.*r1, twothirds*r1+onethird, 2.))
!      PSIW2 = max(0., min(2.*r2, twothirds*r2+onethird, 2.))
!      PSIW3 = max(0., min(2.*r3, twothirds*r3+onethird, 2.))
!.....If flow goes from E to P
!      PSIE1 = max(0., min(2.*r4, twothirds*r4+onethird, 2.))
!      PSIE2 = max(0., min(2.*r5, twothirds*r5+onethird, 2.))
!      PSIE3 = max(0., min(2.*r6, twothirds*r6+onethird, 2.))

!.....PSI for GPL-1/3-alpha-3/2 scheme:  !!!!NEW SCHEME>>>Koren i Ova schema su jako slicne
!.....If flow goes from P to E
!      PSIW1 = max(0., min(1.5*r1, twothirds*r1+onethird, 2.))
!      PSIW2 = max(0., min(1.5*r2, twothirds*r2+onethird, 2.))
!      PSIW3 = max(0., min(1.5*r3, twothirds*r3+onethird, 2.))
!.....If flow goes from E to P
!      PSIE1 = max(0., min(1.5*r4, twothirds*r4+onethird, 2.))
!      PSIE2 = max(0., min(1.5*r5, twothirds*r5+onethird, 2.))
!      PSIE3 = max(0., min(1.5*r6, twothirds*r6+onethird, 2.))

!.....PSI for SMARTER; CHARM NOTABLE; ISNAS
!.....If flow goes from P to E
!      PSIW1 = (r1+abs(r1))*(3*r1+1.)/(2*(r1+1.)**2)
!      PSIW2 = (r2+abs(r2))*(3*r2+1.)/(2*(r2+1.)**2)
!      PSIW3 = (r3+abs(r3))*(3*r3+1.)/(2*(r3+1.)**2)
!.....If flow goes from E to P
!      PSIE1 = (r4+abs(r4))*(3*r4+1.)/(2*(r4+1.)**2)
!      PSIE2 = (r5+abs(r5))*(3*r5+1.)/(2*(r5+1.)**2)
!      PSIE3 = (r6+abs(r6))*(3*r6+1.)/(2*(r6+1.)**2)

!.....PSI for OSPRE
!.....If flow goes from P to E
!      PSIW1 = 1.5*r1*(r1+1.)/(r1**2+r1+1.)
!      PSIW2 = 1.5*r2*(r2+1.)/(r2**2+r2+1.)
!      PSIW3 = 1.5*r3*(r3+1.)/(r3**2+r3+1.)
!.....If flow goes from E to P
!      PSIE1 = 1.5*r4*(r4+1.)/(r4**2+r4+1.)
!      PSIE2 = 1.5*r5*(r5+1.)/(r5**2+r5+1.)
!      PSIE3 = 1.5*r6*(r6+1.)/(r6**2+r6+1.)

!.....PSI for BSOU-BLUI-Chakravarthy-Osher scheme:
!.....If flow goes from P to E
!      PSIW1 = max(0., min(2.*r1,1.))
!      PSIW2 = max(0., min(2.*r2,1.))
!      PSIW3 = max(0., min(2.*r3,1.))
!.....If flow goes from E to P
!      PSIE1 = max(0., min(2.*r4,1.))
!      PSIE2 = max(0., min(2.*r5,1.))
!     PSIE3 = max(0., min(2.*r6,1.))
!=====END MUSCL SCHEME=============================
      END IF
!=====UMIST SCHEME=================================
      IF(LUMIST) THEN
!.....PSI for UMIST scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(2.*r1, 0.75*r1+0.25, 0.25*r1+0.75, 2.))
      PSIW2 = max(0., min(2.*r2, 0.75*r2+0.25, 0.25*r2+0.75, 2.))
      PSIW3 = max(0., min(2.*r3, 0.75*r3+0.25, 0.25*r3+0.75, 2.))
!.....If flow goes from E to P
      PSIE1 = max(0., min(2.*r4, 0.75*r4+0.25, 0.25*r4+0.75, 2.))
      PSIE2 = max(0., min(2.*r5, 0.75*r5+0.25, 0.25*r5+0.75, 2.))
      PSIE3 = max(0., min(2.*r6, 0.75*r6+0.25, 0.25*r6+0.75, 2.))
!=====END UMIST SCHEME=============================
      ENDIF
!=====GAMMA SCHEME================================
      IF(LGAMMA) THEN
!.....PSI for GAMMA scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(r1, 2.*r1/(r1+1.)))
      PSIW2 = max(0., min(r2, 2.*r2/(r2+1.)))
      PSIW3 = max(0., min(r3, 2.*r3/(r3+1.)))
!.....If flow goes from E to P
      PSIE1 = max(0., min(r4, 2.*r4/(r4+1.)))
      PSIE2 = max(0., min(r5, 2.*r5/(r5+1.)))
      PSIE3 = max(0., min(r6, 2.*r6/(r6+1.)))
!=====END GAMMA SCHEME=============================
      END IF

!.....EXPLICIT CONVECTIVE FLUXES FOR HIGH ORDER BOUNDED SCHEMES
!     $Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_e[kg/s] * Phi_e[m/s]$
!     Phi_e is found by extrapolation from upwind nodes, see eq. (3.29) in Sasa's Thesis.
!     Additional multiplication with PSI is application of flux limiters,
!     see eq. (10) in Waterson&Deconinck paper.


!.....Ver.1
!      FUHIGH=CE*(U(INE)*(1+FIF(INE)*PSIE1)-Uinu1*FIF(INE)*PSIE1)+ &
!       CP*(U(INP)*(1+(1-FIF(INW))*PSIW1)-Uinu*(1-FIF(INW))*PSIW1)

!      FVHIGH=CE*(V(INE)*(1+FIF(INE)*PSIE2)-Vinu1*FIF(INE)*PSIE2)+ &
!       CP*(V(INP)*(1+(1-FIF(INW))*PSIW2)-Vinu*(1-FIF(INW))*PSIW2)
       
!      FWHIGH=CE*(W(INE)*(1+FIF(INE)*PSIE3)-Winu1*FIF(INE)*PSIE3)+ &
!       CP*(W(INP)*(1+(1-FIF(INW))*PSIW3)-Winu*(1-FIF(INW))*PSIW3)

!.....Ver.2
!      FUHIGH = CE*(U(INE) + DXE*PSIE1*DFIDX_u1)+ &
!               CP*(U(INP) + DXC*PSIW1*DFIDX_u)
!       mass flux| bounded interpolation of velocity to face |

!      FVHIGH = CE*(V(INE) + DYE*PSIE2*DFIDY_u1)+ &
!               CP*(V(INP) + DYC*PSIW2*DFIDY_u)
       
!      FWHIGH = CE*(W(INE) + DZE*PSIE3*DFIDZ_u1)+ &
!               CP*(W(INP) + DZC*PSIW3*DFIDZ_u)

!.....Ver.3
!......Darwish-Moukalled TVD schemes for unstructured girds, IJHMT, 2003.
      fuhigh = ce*(u(ine) + fxe*psie1*(u(inp)-u(ine)))+ &
               cp*(u(inp) + fxp*psiw1*(u(ine)-u(inp)))
!       mass flux| bounded interpolation of velocity to face |

      fvhigh = ce*(v(ine) + fxe*psie2*(v(inp)-v(ine)))+ &
               cp*(v(inp) + fxp*psiw2*(v(ine)-v(inp)))
       
      fwhigh = ce*(w(ine) + fxe*psie3*(w(inp)-w(ine)))+ &
               cp*(w(inp) + fxp*psiw3*(w(ine)-w(inp)))
!......END: Darwish-Moukalled TVDschemes for unstructured girds, IJHMT, 2003.

!.....END OF BOUNDED HIGH-ORDER SCHEMES
      END IF 

!
!.....FINALE: EXPLICIT PART OF FLUXES!
!
!=======================================================================
      su(inp) = su(inp)-fuhigh+fdue
      sv(inp) = sv(inp)-fvhigh+fdve
      sw(inp) = sw(inp)-fwhigh+fdwe
!=======================================================================

!
!......[Vectorization Procedure: ]
!----------------------------------------------------
      bp(ine) = fuhigh-fdue
      bt(ine) = fvhigh-fdve
      bb(ine) = fwhigh-fdwe
!----------------------------------------------------

      ! 1/ap emulation for explicit scheme, needed for const mflux
      apu(inp) = apu(inp)-fuhigh+fdue

      end do
      end do
      end do


      do k=2,nkmm
      do i=2,nimm
      do j=2,njmm
      inp=lk(k)+li(i)+j+idew

      su(inp)=su(inp)+bp(inp)
      sv(inp)=sv(inp)+bt(inp)
      sw(inp)=sw(inp)+bb(inp)

      end do !j-loop
      end do !i-loop
      end do !k-loop


      do ijk=icst,icen
      bp(ijk)=0.0d0
      bt(ijk)=0.0d0
      bb(ijk)=0.0d0
      end do


      return
      end
