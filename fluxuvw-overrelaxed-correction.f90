!***********************************************************************
!
      subroutine fluxuvw(idew,idns,idtb, &
                         arvx,arvy,arvz, &
                         fif,acfe,acfw,fcf)
!
!***********************************************************************
!
!  fluid flow trough control volumes:
!   _______________________________________________________________
!  |           |            |            |            |            |
!  |           |            |            |            |            |
!  |     o ww ===>   o w   ===>   o p   ===>   o e   ===>   o  ee  |
!  |           |            |            |            |            |
!  |           |            |            |            |            |
!  |___________|____________|____________|____________|____________|
!
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
      real(prec), dimension(nxyza) :: fif,acfe,acfw,fcf 
!
!     local variables
!
      integer :: i, j, k, ijk, inp,       &
                 ine,inw,inn,ins,inb,int, &
                 inbs,inee
      integer :: dir
      integer :: ist,jst,kst,iend,jend,kend
      real(prec) :: fxw,fxpw,fxee,fxew, &
                    dxc,dyc,dzc,dxe,dye,dze
      !real(prec) :: dxu,dyu,dzu,dxu1,dyu1,dzu1,uinu,vinu,winu,uinu1,vinu1,winu1
      real(prec) :: dfidx_u,dfidx_u1, &
                    dfidy_u,dfidy_u1, &
                    dfidz_u,dfidz_u1
!      real(prec) :: dfidx_f,dfidy_f,dfidz_f
!      real(prec) :: dxet,dxzd,dyet,dyzd,dzet,dzzd
      real(prec) :: gam,fxe,fxp,are,vole,game,&
                    ue,ve,we, & 
                    de,flcf,ce,cp, &
                    g11, g12, g21, g22
      real(prec) :: ae1,aw1,shigh1,shigh2,shigh3
      real(prec) :: r1,r2,r3,r4,r5,r6,                   &
                    psie1,psie2,psie3,psiw1,psiw2,psiw3, &
                    fuuds,fvuds,fwuds,fuhigh,fvhigh,fwhigh
      real(prec) :: xf,yf,zf,xi,yi,zi,arx,ary,arz, &
                    duxi,duyi,duzi, &
                    dvxi,dvyi,dvzi, &
                    dwxi,dwyi,dwzi
      real(prec) :: xpn,ypn,zpn
      real(prec) :: nxx,nyy,nzz,ixi1,ixi2,ixi3,dpn,costheta,costn
      real(prec) :: fdue,fdve,fdwe,fdui,fdvi,fdwi
      real(prec) :: onethird, twothirds
      real(prec) :: d2x,d2y,d2z,d1x,d1y,d1z
      real(prec) :: duxii,dvxii,dwxii, &
                    duyii,dvyii,dwyii, &
                    duzii,dvzii,dwzii
!----------------------------------------------------------------------

      gam=gds(iu) 

      onethird = 1./3.0d0
      twothirds = 2./3.0d0

!.....initialize values for high order conv. fluxes
      psie1=0.0d0
      psie2=0.0d0
      psie3=0.0d0
      psiw1=0.0d0
      psiw2=0.0d0
      psiw3=0.0d0         
      fuhigh=0.0d0
      fvhigh=0.0d0
      fwhigh=0.0d0

!.....choose direction for gradients; if dir=1 we use dfidx or grad(inp,1), 2=>dfidy, 3=>dfidz
      if(idew.eq.nj)  dir = 1
      if(idew.eq.1)   dir = 2
      if(idew.eq.nij) dir = 3

      kst = 2
      kend = nkmm
      ist = 2
      iend = nimm
      jst = 2
      jend = njmm

      ! kod dir=3 hard kodirano je da su top i bottom zidovi pa ne radi fluks kroz te zidove nego ih racunamo kao BC
      if(dir == 1) then
        jst = 3
        kst = 3
      elseif(dir == 2) then
        ist = 3
        kst = 3
      elseif(dir == 3) then
        ist = 3
        jst = 3
        ! kst = 3
        ! kend = nkmm-1
      endif
!
!.....calculate east,top,north  cell face
      do k=kst,kend
      do i=ist,iend
      do j=jst,jend

      inp=lk(k)+li(i)+j

      ine=inp+idew
      inw=inp-idew
      inn=inp+idns
      ins=inp-idns
      inb=inp-idtb
      int=inp+idtb

      inee=inp+2*idew
      inbs=inb-idns
!
!.....interpolation factors
      fxe=fif(inp) 
      fxp=1.0d0-fxe

      fxpw = fif(inw)
      fxw = 1.0d0-fxpw
      fxee = fif(ine)
      fxew = 1.0d0-fxee 

      xpn=xc(ine)-xc(inp)
      ypn=yc(ine)-yc(inp)
      zpn=zc(ine)-zc(inp)


!.....distance from p to neighbor n
      dpn=sqrt(xpn**2+ypn**2+zpn**2)     

!.....components of the unit vector i_ksi
      ixi1=xpn/dpn
      ixi2=ypn/dpn
      ixi3=zpn/dpn

!.....precomputed face areas
      arx=arvx(inp)
      ary=arvy(inp)
      arz=arvz(inp)

!.....cell face area
      are=sqrt(arx**2+ary**2+arz**2)

!.....unit vectors of the normal
      nxx=arx/are
      nyy=ary/are
      nzz=arz/are
!.....angle between vectorsa n and i_xi - we need cosine
      costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

!.....relaxation factor for higher-order cell face gradient
      ! minimal correction: nrelax = +1 :
      !costn = costheta
      ! orthogonal correction: nrelax =  0 : 
      costn = 1.0d0
      ! over-relaxed approach: nrelax = -1 :
      !costn = 1./costheta
      ! in general, nrelax can be any signed integer from some 
      ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
      !costn = costheta**nrelax

!.....first .(second x third) = vol
      vole=xpn*arx+ypn*ary+zpn*arz
!.....overrelaxed correction vector d2, where s=dpn+d2
      d1x = costn
      d1y = costn
      d1z = costn
      !
      d2x = xpn*costn
      d2y = ypn*costn
      d2z = zpn*costn

!.....cell face viscosity
      game=(vis(inp)*fxp+vis(ine)*fxe)

!++++velocities at cell face center and explicit diffusion fluxes+++++++

!.....coordinates of cell-face center
      xf = 0.25*(x(inp)+x(ins)+x(inb)+x(inbs))
      yf = 0.25*(y(inp)+y(ins)+y(inb)+y(inbs))
      zf = 0.25*(z(inp)+z(ins)+z(inb)+z(inbs))

!.....coordinates of point e'
      xi=xc(inp)*fxp+xc(ine)*fxe
      yi=yc(inp)*fxp+yc(ine)*fxe
      zi=zc(inp)*fxp+zc(ine)*fxe

!.....interpolate gradients defined at cv centers to faces
      duxi = gradu(1,inp)*fxp+gradu(1,ine)*fxe
      duyi = gradu(2,inp)*fxp+gradu(2,ine)*fxe
      duzi = gradu(3,inp)*fxp+gradu(3,ine)*fxe

!.....du/dx_i interpolated at cell face:
      duxii = duxi*d1x + arx/vole*( u(ine)-u(inp)-duxi*d2x-duyi*d2y-duzi*d2z ) 
      duyii = duyi*d1y + ary/vole*( u(ine)-u(inp)-duxi*d2x-duyi*d2y-duzi*d2z ) 
      duzii = duzi*d1z + arz/vole*( u(ine)-u(inp)-duxi*d2x-duyi*d2y-duzi*d2z ) 


      dvxi = gradv(1,inp)*fxp+gradv(1,ine)*fxe
      dvyi = gradv(2,inp)*fxp+gradv(2,ine)*fxe
      dvzi = gradv(3,inp)*fxp+gradv(3,ine)*fxe

!.....dv/dx_i interpolated at cell face:
      dvxii = dvxi*d1x + arx/vole*( v(ine)-v(inp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 
      dvyii = dvyi*d1y + ary/vole*( v(ine)-v(inp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 
      dvzii = dvzi*d1z + arz/vole*( v(ine)-v(inp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 


      dwxi = gradw(1,inp)*fxp+gradw(1,ine)*fxe
      dwyi = gradw(2,inp)*fxp+gradw(2,ine)*fxe
      dwzi = gradw(3,inp)*fxp+gradw(3,ine)*fxe

!.....dw/dx_i interpolated at cell face:
      dwxii = dwxi*d1x + arx/vole*( w(ine)-w(inp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 
      dwyii = dwyi*d1y + ary/vole*( w(ine)-w(inp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 
      dwzii = dwzi*d1z + arz/vole*( w(ine)-w(inp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 


      ! if (lcds4.eq.1.and.(idew.eq.nj.or.idew.eq.1)) then
      ! duxii = 1./(24.0d0*xpn+small)*(27.0d0*u(ine)-27.0d0*u(inp)+u(inw)-u(inee))
      ! duyii = 1./(24.0d0*ypn+small)*(27.0d0*u(ine)-27.0d0*u(inp)+u(inw)-u(inee))
      ! duzii = 1./(24.0d0*zpn+small)*(27.0d0*u(ine)-27.0d0*u(inp)+u(inw)-u(inee))

      ! dvxii = 1./(24.0d0*xpn+small)*(27.0d0*v(ine)-27.0d0*v(inp)+v(inw)-v(inee))
      ! dvyii = 1./(24.0d0*ypn+small)*(27.0d0*v(ine)-27.0d0*v(inp)+v(inw)-v(inee))
      ! dvzii = 1./(24.0d0*zpn+small)*(27.0d0*v(ine)-27.0d0*v(inp)+v(inw)-v(inee))

      ! dwxii = 1./(24.0d0*xpn+small)*(27.0d0*w(ine)-27.0d0*w(inp)+w(inw)-w(inee))
      ! dwyii = 1./(24.0d0*ypn+small)*(27.0d0*w(ine)-27.0d0*w(inp)+w(inw)-w(inee))
      ! dwzii = 1./(24.0d0*zpn+small)*(27.0d0*w(ine)-27.0d0*w(inp)+w(inw)-w(inee))
      ! endif


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     we calculate explicit and implicit diffsion fde and fdi,
!     later se put their difference (fde-fdi) to rhs vector
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! explicit diffussion 
       fdue = game*( (duxii+duxii)*arx + (duyii+dvxii)*ary + (duzii+dwxii)*arz )
       fdve = game*( (duyii+dvxii)*arx + (dvyii+dvyii)*ary + (dvzii+dwyii)*arz )
       fdwe = game*( (duzii+dwxii)*arx + (dwyii+dvzii)*ary + (dwzii+dwzii)*arz )
       
       ! implicit diffussion 
       fdui = game*are/dpn*(duxi*xpn+duyi*ypn+duzi*zpn)
       fdvi = game*are/dpn*(dvxi*xpn+dvyi*ypn+dvzi*zpn)
       fdwi = game*are/dpn*(dwxi*xpn+dwyi*ypn+dwzi*zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++end: velocities at cell face center and explicit diffusion fluxes+++++++


!.....difusion coefficient
      de=game*are/dpn

!.....convection fluxes - uds
!     if flow goes p=>e cp=flcf, ce=0.
!     if flow goes e=>p ce=flcf, cp=0.
      flcf=fcf(inp)
      ce=min(flcf,zero) 
      cp=max(flcf,zero)
!
!.....coefficients ae(p) and aw(e) due to uds
!
      acfe(inp)=de-ce
      acfw(ine)=de+cp
!
!.....explicit convective fluxes for uds
!
      fuuds=cp*u(inp)+ce*u(ine)
      fvuds=cp*v(inp)+ce*v(ine)
      fwuds=cp*w(inp)+ce*w(ine)

!=====start cds scheme===============================
      if(lcds) then
!
!.....explicit convective fluxes for cds
!
!        |________ue'_________|_______________ucorr___________________|
      ue=u(inp)*fxp+u(ine)*fxe+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))
!        |________ve'_________|_______________vcorr___________________|
      ve=v(inp)*fxp+v(ine)*fxe+(dvxi*(xf-xi)+dvyi*(yf-yi)+dvzi*(zf-zi))
!        |________we'_________|_______________wcorr___________________|
      we=w(inp)*fxp+w(ine)*fxe+(dwxi*(xf-xi)+dwyi*(yf-yi)+dwzi*(zf-zi))

      fuhigh=flcf*ue
      fvhigh=flcf*ve
      fwhigh=flcf*we


!=====end cds scheme================================
      end if

!=====start cds4 scheme===============================
      if(lcds4.and.(idew.eq.nj.or.idew.eq.1)) then
!
!.....explicit convective fluxes for cds
!
      ! forth ordr.
      ue = 1./48.0d0*(27.0d0*u(inp)+27.0d0*u(ine)-3.0d0*u(inw)-3.0d0*u(inee))
      ve = 1./48.0d0*(27.0d0*v(inp)+27.0d0*v(ine)-3.0d0*v(inw)-3.0d0*v(inee))
      we = 1./48.0d0*(27.0d0*w(inp)+27.0d0*w(ine)-3.0d0*w(inw)-3.0d0*w(inee))

      fuhigh=flcf*ue
      fvhigh=flcf*ve
      fwhigh=flcf*we


!=====end cds scheme================================
      end if

!=====start luds scheme=============================
      if(lluds) then !!! iskulirao sam ga za neortogonalne mreze

!.....explicit convective fluxes for luds scheme
!     $flux_high_order_scheme[n] = mass_flow_rate_in_cell_face_e[kg/s] * phi_e[m/s]$
!     phi_e is found by extrapolation from upwind nodes, see eq. (3.29) in sasa's thesis.

!.....distance from considered node to the face
!     in the waterson & deconinck article there is a term 'deltaxc/2', whole this term we denote dxc
!.....if flow goes from p to e
      dxc = xf - xc(inp)
      dyc = yf - yc(inp)
      dzc = zf - zc(inp)
!.....if flow goes from e to p
      dxe = xc(ine) - xf
      dye = yc(ine) - yf
      dze = zc(ine) - zf

!.....interpolate gradients defined at cv centers to faces
      dfidx_u = gradu(dir,inw)*fxw+gradu(dir,inp)*fxpw
      dfidy_u = gradv(dir,inw)*fxw+gradv(dir,inp)*fxpw
      dfidz_u = gradw(dir,inw)*fxw+gradw(dir,inp)*fxpw
!
      dfidx_u1 = gradu(dir,inee)*fxee+gradu(dir,ine)*fxew
      dfidy_u1 = gradv(dir,inee)*fxee+gradv(dir,ine)*fxew
      dfidz_u1 = gradw(dir,inee)*fxee+gradw(dir,ine)*fxew

      fuhigh = ce*(u(ine) + dxe*dfidx_u1)+ &
               cp*(u(inp) + dxc*dfidx_u)
!       mass flux| interpolation of velocity to face |

      fvhigh = ce*(v(ine) + dye*dfidy_u1)+ &
               cp*(v(inp) + dyc*dfidy_u)
       
      fwhigh = ce*(w(ine) + dze*dfidz_u1)+ &
               cp*(w(inp) + dzc*dfidz_u)

!=====end luds scheme====================================
      end if

!=====start quick scheme=================================
      if(lquds) then

      shigh1=0.0d0
      shigh2=0.0d0
      shigh3=0.0d0
      aw1=0.0d0
      ae1=0.0d0

      if(k.gt.2.and.k.lt.nk-1.and. &
         j.gt.2.and.j.lt.nj-1.and. &
         i.gt.2.and.i.lt.ni-1) then

!.....find coefficients associated with quadratic function
      g11=((2-fif(inp-idew))*fif(inp)**2)/  &
           (1+fif(inp)-fif(inp-idew))
      g12=((1-fif(inp))*(1-fif(inp-idew))**2)/  &
           (1+fif(inp)-fif(inp-idew))
      g21=((1+fif(inp+idew))*(1-fif(inp))**2)/  &
           (1+fif(inp+idew)-fif(inp))
      g22=(fif(inp+idew)**2*fif(inp))/ &
           (1+fif(inp+idew)-fif(inp))

!.....explicit convective fluxes for quick scheme
!     $flux_high_order_scheme[n] = mass_flow_rate_in_cell_face_e[kg/s] * phi_e[m/s]$
!     phi_e is found by extrapolation from upwind nodes, see eq. (3.30) in sasa's thesis.

      fuhigh=ce*(g21*u(inp)+(1-g21+g22)*u(ine)-g22*u(inee))+ &
             cp*(g11*u(ine)+(1-g11+g12)*u(inp)-g12*u(inw))


      fvhigh=ce*(g21*v(inp)+(1-g21+g22)*v(ine)-g22*v(inee))+ &
             cp*(g11*v(ine)+(1-g11+g12)*v(inp)-g12*v(inw))
       
      fwhigh=ce*(g21*w(inp)+(1-g21+g22)*w(ine)-g22*w(inee))+ &
             cp*(g11*w(ine)+(1-g11+g12)*w(inp)-g12*w(inw))

      end if

      end if
!=============end quick scheme========================
!--------------------------------------------------------------------------------------------
!     bounded high-order convective schemes (waterson & deconinck jcp 224 (2007) pp. 182-207)
!     notes:
!     smart- converges slow, requires stronger underrelaxation (like urf = 0.2), otherwise
!     probably smallest error when converges
!     avl-smart - better convergence than smart, works fine with urf = 0.5
!     muscl - good convergence, urf = 0.6, even 0.8, the most robust one.
!     ...
!--------------------------------------------------------------------------------------------
      if(flux_limiter) then

! !+++++patch-version 1++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! !+++++find r's. this is universal for all schemes. non-orthogonal grid! +++++++
! !.....distances of fictitious upwind point u from 'real' upwind cell center.
! !.....if flow goes from p to e, distance from the new node u from closest cell-center
!       dxu = (xc(inp) - (xc(ine)-xc(inp))) - xc(inw)
!       dyu = (yc(inp) - (yc(ine)-yc(inp))) - yc(inw)
!       dzu = (zc(inp) - (zc(ine)-zc(inp))) - zc(inw)
! !.....if flow goes from p to e, distance from the new node u from closest cell-center
!       dxu1 = (xc(ine) + (xc(inp)-xc(inw))) - xc(inee)
!       dyu1 = (yc(ine) + (yc(inp)-yc(inw))) - yc(inee)
!       dzu1 = (zc(ine) + (zc(inp)-zc(inw))) - zc(inee)
! 
! !.....variable value at fictitious upwind point
!       uinu=u(inw)+dxu*gradu(1,inw)+dyu*gradu(2,inw)+dzu*gradu(3,inw)
!       vinu=v(inw)+dxu*gradv(1,inw)+dyu*gradv(2,inw)+dzu*gradv(3,inw)
!       winu=w(inw)+dxu*gradw(1,inw)+dyu*gradw(2,inw)+dzu*gradw(3,inw)
! !
!       uinu1=u(inee)+dxu1*gradu(1,inee)+dyu1*gradu(2,inee)+dzu1*gradu(3,inee)
!       vinu1=v(inee)+dxu1*gradv(1,inee)+dyu1*gradv(2,inee)+dzu1*gradv(3,inee)
!       winu1=w(inee)+dxu1*gradw(1,inee)+dyu1*gradw(2,inee)+dzu1*gradw(3,inee)
! 
! !.....gradient ratio using value at fictitious upwind point:
! !.....if flow goes from p to e
!       r1 = (u(ine)-u(inp))/(u(inp)-uinu)
!       r2 = (v(ine)-v(inp))/(v(inp)-vinu)
!       r3 = (w(ine)-w(inp))/(w(inp)-winu)
! !.....if flow goes from e to p
!       r4 = (u(inp)-u(ine))/(u(ine)-uinu1)
!       r5 = (v(inp)-v(ine))/(v(ine)-vinu1)
!       r6 = (w(inp)-w(ine))/(w(ine)-winu1)
! !+++++patch-version 1++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++patch-version 2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!.....distance from considered node to the face
!     in the waterson & deconinck article there is a term 'deltaxc/2', whole this term we denote dxc
!.....if flow goes from p to e
!      dxc = ( 0.25*(x(inp)+x(ins)+x(inb)+x(inbs)) ) - xc(inp)
!      dyc = ( 0.25*(y(inp)+y(ins)+y(inb)+y(inbs)) ) - yc(inp)
!      dzc = ( 0.25*(z(inp)+z(ins)+z(inb)+z(inbs)) ) - zc(inp)
!.....if flow goes from e to p
!      dxe = xc(ine) - ( 0.25*(x(inp)+x(ins)+x(inb)+x(inbs)) )
!      dye = yc(ine) - ( 0.25*(y(inp)+y(ins)+y(inb)+y(inbs)) )
!      dze = zc(ine) - ( 0.25*(z(inp)+z(ins)+z(inb)+z(inbs)) )

!.....interpolate gradients defined at cv centers to faces
!      dfidx_f = gradu(dir,inp)*fxp+gradu(dir,ine)*fxe 
!      dfidy_f = gradv(dir,inp)*fxp+gradv(dir,ine)*fxe
!      dfidz_f = gradw(dir,inp)*fxp+gradw(dir,ine)*fxe
!
!      dfidx_u = gradu(dir,inw)*fxw+gradu(dir,inp)*fxpw
!      dfidy_u = gradv(dir,inw)*fxw+gradv(dir,inp)*fxpw
!      dfidz_u = gradw(dir,inw)*fxw+gradw(dir,inp)*fxpw
!
!      dfidx_u1 = gradu(dir,inee)*fxee+gradu(dir,ine)*fxew
!      dfidy_u1 = gradv(dir,inee)*fxee+gradv(dir,ine)*fxew
!      dfidz_u1 = gradw(dir,inee)*fxee+gradw(dir,ine)*fxew

!.... find r's. this is universal for all schemes.
!.....if flow goes from p to e
!      r1 = dfidx_f/dfidx_u  
!      r2 = dfidy_f/dfidy_u 
!      r3 = dfidz_f/dfidz_u  
!.....if flow goes from e to p
!      r4 = dfidx_f/dfidx_u1  
!      r5 = dfidy_f/dfidy_u1 
!      r6 = dfidz_f/dfidz_u1  
!+++++patch-version 2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++patch-version 3++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!.... find r's. this is universal for all schemes.
!.....if flow goes from p to e
      r1 = (2*gradu(1,inp)*xpn + 2*gradu(2,inp)*ypn + 2*gradu(3,inp)*zpn)/(u(ine)-u(inp)) - 1.0d0  
      r2 = (2*gradv(1,inp)*xpn + 2*gradv(2,inp)*ypn + 2*gradv(3,inp)*zpn)/(v(ine)-v(inp)) - 1.0d0 
      r3 = (2*gradw(1,inp)*xpn + 2*gradw(2,inp)*ypn + 2*gradw(3,inp)*zpn)/(w(ine)-w(inp)) - 1.0d0 
!.....if flow goes from e to p
      r4 = (2*gradu(1,ine)*xpn + 2*gradu(2,ine)*ypn + 2*gradu(3,ine)*zpn)/(u(inp)-u(ine)) - 1.0d0 
      r5 = (2*gradv(1,ine)*xpn + 2*gradv(2,ine)*ypn + 2*gradv(3,ine)*zpn)/(v(inp)-v(ine)) - 1.0d0 
      r6 = (2*gradw(1,ine)*xpn + 2*gradw(2,ine)*ypn + 2*gradw(3,ine)*zpn)/(w(inp)-w(ine)) - 1.0d0  
!+++++patch-version 3++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!=====smart scheme================================
      if(lsmart) then
!.....psi for smart scheme:
!.....if flow goes from p to e
      psiw1 = max(0., min(2.*r1, 0.75*r1+0.25, 4.))
      psiw2 = max(0., min(2.*r2, 0.75*r2+0.25, 4.))
      psiw3 = max(0., min(2.*r3, 0.75*r3+0.25, 4.))
!.....if flow goes from e to p
      psie1 = max(0., min(2.*r4, 0.75*r4+0.25, 4.))
      psie2 = max(0., min(2.*r5, 0.75*r5+0.25, 4.))
      psie3 = max(0., min(2.*r6, 0.75*r6+0.25, 4.))
!=====end smart scheme=============================
      end if
!
!=====avl-smart scheme=============================
      if(lavl) then
!.....psi for avl-smart scheme:
!.....if flow goes from p to e
      psiw1 = max(0., min(1.5*r1, 0.75*r1+0.25, 2.5))
      psiw2 = max(0., min(1.5*r2, 0.75*r2+0.25, 2.5))
      psiw3 = max(0., min(1.5*r3, 0.75*r3+0.25, 2.5))
!.....if flow goes from e to p
      psie1 = max(0., min(1.5*r4, 0.75*r4+0.25, 2.5))
      psie2 = max(0., min(1.5*r5, 0.75*r5+0.25, 2.5))
      psie3 = max(0., min(1.5*r6, 0.75*r6+0.25, 2.5))
!=====end avl-smart scheme==========================
      end if
!
!=====muscl scheme=================================
      if(lmuscl) then
!.....psi for muscl scheme:
!.....if flow goes from p to e
      psiw1 = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
      psiw2 = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
      psiw3 = max(0., min(2.*r3, 0.5*r3+0.5, 2.))
!.....if flow goes from e to p
      psie1 = max(0., min(2.*r4, 0.5*r4+0.5, 2.))
      psie2 = max(0., min(2.*r5, 0.5*r5+0.5, 2.))
      psie3 = max(0., min(2.*r6, 0.5*r6+0.5, 2.))


!=====end muscl scheme=============================
      end if
!=====umist scheme=================================
      if(lumist) then
!.....psi for umist scheme:
!.....if flow goes from p to e
      psiw1 = max(0., min(2.*r1, 0.75*r1+0.25, 0.25*r1+0.75, 2.))
      psiw2 = max(0., min(2.*r2, 0.75*r2+0.25, 0.25*r2+0.75, 2.))
      psiw3 = max(0., min(2.*r3, 0.75*r3+0.25, 0.25*r3+0.75, 2.))
!.....if flow goes from e to p
      psie1 = max(0., min(2.*r4, 0.75*r4+0.25, 0.25*r4+0.75, 2.))
      psie2 = max(0., min(2.*r5, 0.75*r5+0.25, 0.25*r5+0.75, 2.))
      psie3 = max(0., min(2.*r6, 0.75*r6+0.25, 0.25*r6+0.75, 2.))
!=====end umist scheme=============================
      endif
!=====gamma scheme================================
      if(lgamma) then
!.....psi for gamma scheme:
!.....if flow goes from p to e
      psiw1 = max(0., min(r1, 2.*r1/(r1+1.)))
      psiw2 = max(0., min(r2, 2.*r2/(r2+1.)))
      psiw3 = max(0., min(r3, 2.*r3/(r3+1.)))
!.....if flow goes from e to p
      psie1 = max(0., min(r4, 2.*r4/(r4+1.)))
      psie2 = max(0., min(r5, 2.*r5/(r5+1.)))
      psie3 = max(0., min(r6, 2.*r6/(r6+1.)))
!=====end gamma scheme=============================
      end if

!.....explicit convective fluxes for high order bounded schemes
!     $flux_high_order_scheme[n] = mass_flow_rate_in_cell_face_e[kg/s] * phi_e[m/s]$
!     phi_e is found by extrapolation from upwind nodes, see eq. (3.29) in sasa's thesis.
!     additional multiplication with psi is application of flux limiters,
!     see eq. (10) in waterson&deconinck paper.


!.....ver.1
!      fuhigh=ce*(u(ine)*(1+fif(ine)*psie1)-uinu1*fif(ine)*psie1)+ &
!       cp*(u(inp)*(1+(1-fif(inw))*psiw1)-uinu*(1-fif(inw))*psiw1)

!      fvhigh=ce*(v(ine)*(1+fif(ine)*psie2)-vinu1*fif(ine)*psie2)+ &
!       cp*(v(inp)*(1+(1-fif(inw))*psiw2)-vinu*(1-fif(inw))*psiw2)
       
!      fwhigh=ce*(w(ine)*(1+fif(ine)*psie3)-winu1*fif(ine)*psie3)+ &
!       cp*(w(inp)*(1+(1-fif(inw))*psiw3)-winu*(1-fif(inw))*psiw3)

!.....ver.2
!      fuhigh = ce*(u(ine) + dxe*psie1*dfidx_u1)+ &
!               cp*(u(inp) + dxc*psiw1*dfidx_u)
!       mass flux| bounded interpolation of velocity to face |

!      fvhigh = ce*(v(ine) + dye*psie2*dfidy_u1)+ &
!               cp*(v(inp) + dyc*psiw2*dfidy_u)
       
!      fwhigh = ce*(w(ine) + dze*psie3*dfidz_u1)+ &
!               cp*(w(inp) + dzc*psiw3*dfidz_u)

!.....ver.3
!......darwish-moukalled tvd schemes for unstructured girds, ijhmt, 2003.
      fuhigh = ce*(u(ine) + fxe*psie1*(u(inp)-u(ine)))+ &
               cp*(u(inp) + fxp*psiw1*(u(ine)-u(inp)))
!       mass flux| bounded interpolation of velocity to face |

      fvhigh = ce*(v(ine) + fxe*psie2*(v(inp)-v(ine)))+ &
               cp*(v(inp) + fxp*psiw2*(v(ine)-v(inp)))
       
      fwhigh = ce*(w(ine) + fxe*psie3*(w(inp)-w(ine)))+ &
               cp*(w(inp) + fxp*psiw3*(w(ine)-w(inp)))
!......end: darwish-moukalled tvdschemes for unstructured girds, ijhmt, 2003.

!.....end of bounded high-order schemes
      end if 

!
!.....explicit part of diffusion fluxes and sources due to deffered correction
!     for all schemes!
!
      su(inp)=su(inp)+gam*(fuuds-fuhigh)+fdue-fdui
      sv(inp)=sv(inp)+gam*(fvuds-fvhigh)+fdve-fdvi
      sw(inp)=sw(inp)+gam*(fwuds-fwhigh)+fdwe-fdwi
!----------------------------------------------------
!......[vectorization procedure: ]
!----------------------------------------------------
      bp(ine)=-gam*(fuuds-fuhigh)-fdue+fdui
      bt(ine)=-gam*(fvuds-fvhigh)-fdve+fdvi
      bb(ine)=-gam*(fwuds-fwhigh)-fdwe+fdwi
!----------------------------------------------------
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



!.....psi for koren scheme:
!.....if flow goes from p to e
!      psiw1 = max(0., min(2.*r1, twothirds*r1+onethird, 2.))
!      psiw2 = max(0., min(2.*r2, twothirds*r2+onethird, 2.))
!      psiw3 = max(0., min(2.*r3, twothirds*r3+onethird, 2.))
!.....if flow goes from e to p
!      psie1 = max(0., min(2.*r4, twothirds*r4+onethird, 2.))
!      psie2 = max(0., min(2.*r5, twothirds*r5+onethird, 2.))
!      psie3 = max(0., min(2.*r6, twothirds*r6+onethird, 2.))

!.....psi for gpl-1/3-alpha-3/2 scheme:  !!!!new scheme>>>koren i ova schema su jako slicne
!.....if flow goes from p to e
!      psiw1 = max(0., min(1.5*r1, twothirds*r1+onethird, 2.))
!      psiw2 = max(0., min(1.5*r2, twothirds*r2+onethird, 2.))
!      psiw3 = max(0., min(1.5*r3, twothirds*r3+onethird, 2.))
!.....if flow goes from e to p
!      psie1 = max(0., min(1.5*r4, twothirds*r4+onethird, 2.))
!      psie2 = max(0., min(1.5*r5, twothirds*r5+onethird, 2.))
!      psie3 = max(0., min(1.5*r6, twothirds*r6+onethird, 2.))

!.....psi for smarter; charm notable; isnas
!.....if flow goes from p to e
!      psiw1 = (r1+abs(r1))*(3*r1+1.)/(2*(r1+1.)**2)
!      psiw2 = (r2+abs(r2))*(3*r2+1.)/(2*(r2+1.)**2)
!      psiw3 = (r3+abs(r3))*(3*r3+1.)/(2*(r3+1.)**2)
!.....if flow goes from e to p
!      psie1 = (r4+abs(r4))*(3*r4+1.)/(2*(r4+1.)**2)
!      psie2 = (r5+abs(r5))*(3*r5+1.)/(2*(r5+1.)**2)
!      psie3 = (r6+abs(r6))*(3*r6+1.)/(2*(r6+1.)**2)

!.....psi for ospre
!.....if flow goes from p to e
!      psiw1 = 1.5*r1*(r1+1.)/(r1**2+r1+1.)
!      psiw2 = 1.5*r2*(r2+1.)/(r2**2+r2+1.)
!      psiw3 = 1.5*r3*(r3+1.)/(r3**2+r3+1.)
!.....if flow goes from e to p
!      psie1 = 1.5*r4*(r4+1.)/(r4**2+r4+1.)
!      psie2 = 1.5*r5*(r5+1.)/(r5**2+r5+1.)
!      psie3 = 1.5*r6*(r6+1.)/(r6**2+r6+1.)

!.....psi for bsou-blui-chakravarthy-osher scheme:
!.....if flow goes from p to e
!      psiw1 = max(0., min(2.*r1,1.))
!      psiw2 = max(0., min(2.*r2,1.))
!      psiw3 = max(0., min(2.*r3,1.))
!.....if flow goes from e to p
!      psie1 = max(0., min(2.*r4,1.))
!      psie2 = max(0., min(2.*r5,1.))
!     psie3 = max(0., min(2.*r6,1.))