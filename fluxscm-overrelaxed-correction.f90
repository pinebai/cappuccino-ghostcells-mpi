!***********************************************************************
!
      subroutine fluxscm(idew,idns,idtb,fi,gradfi,ifi, &
                         arvx,arvy,arvz, &
                         fif,acfe,acfw,fcf)
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
      use omega_turb_models
      use fieldManipulation

      implicit none
!
!***********************************************************************
!
      integer, intent(in) :: idew, idns, idtb, ifi
      real(prec), dimension(nxyza) :: arvx,arvy,arvz
      real(prec), dimension(nxyza) :: fif
      real(prec), dimension(nxyza) :: acfe, acfw
      real(prec), dimension(nxyza) :: fcf, fi
      real(prec), dimension(3,nxyza) :: gradfi

!     local variables

      integer :: i, j, k, inp, ine, &
                 ins,inb,inbs
      real(prec) :: gam,prtr,fxe,fxp,are,vole, &
                    viste,game,de, &
                    ce,cp,fii,fm
      real(prec) :: xf,yf,zf,xi,yi,zi
      real(prec) :: arx,ary,arz
      real(prec) :: xpn,ypn,zpn
      real(prec) :: nxx,nyy,nzz,ixi1,ixi2,ixi3,dpn,costheta,costn
      real(prec) :: fdfie,fdfii,fcfie,fcfii,ffic,suadd
      real(prec) :: d2x,d2y,d2z,d1x,d1y,d1z
      real(prec) :: dfixi,dfiyi,dfizi
      real(prec) :: dfixii,dfiyii,dfizii
      real(prec) :: r1,r2,psie,psiw

!.....Blending coefficient for convection
      gam=gds(ifi)

!.....Usually it is constant:
      prtr=prtinv(ifi)
    
      dfixi = 0.0d0
      dfiyi = 0.0d0
      dfizi = 0.0d0
       
      ! Loop over faces
      do k=2,nkmm
      do i=2,nimm
      do j=2,njmm

      inp=lk(k)+li(i)+j

      ine=inp+idew
      ins=inp-idns
      inb=inp-idtb
      inbs=inb-idns

!.....Interpolation factor
      fxe=fif(inp)
      fxp=1.-fxe


!.....Distance between cell centers - components
      xpn=xc(ine)-xc(inp)
      ypn=yc(ine)-yc(inp)
      zpn=zc(ine)-zc(inp)

!.....Distance from P to neighbor N - magnitude
      dpn=sqrt(xpn**2+ypn**2+zpn**2)     

!.....Components of the unit vector i_ksi
      ixi1=xpn/dpn
      ixi2=ypn/dpn
      ixi3=zpn/dpn

!.....Precomputed face areas
      arx=arvx(inp)
      ary=arvy(inp)
      arz=arvz(inp)

!.....CELL FACE AREA
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
      ! reasonable interval [-nrelax,nrelax] (or even real number): 
      !costn = costheta**nrelax

!.....FIRST .(SECOND X THIRD) = VOL
      vole=xpn*arx+ypn*ary+zpn*arz
!.....Overrelaxed correction vector d2, where S=dpn+d2
      d1x = costn
      d1y = costn
      d1z = costn
      !
      d2x = xpn*costn
      d2y = ypn*costn
      d2z = zpn*costn

!.....For Menter SST model:
      if (sst.or.sas.or.earsm_wj.or.earsm_m) then
        if(ifi.eq.ite) prtr=prtinv_te(inp)*fxp+prtinv_te(ine)*fxe
        if(ifi.eq.ied) prtr=prtinv_ed(inp)*fxp+prtinv_ed(ine)*fxe
      end if

!.....Cell face diffussion coefficient
      viste=(vis(inp)-viscos)*fxp+(vis(ine)-viscos)*fxe
      game=(viste*prtr+viscos)

!.....Coordinates of cell-face center
      xf = 0.25*(x(inp)+x(ins)+x(inb)+x(inbs))
      yf = 0.25*(y(inp)+y(ins)+y(inb)+y(inbs))
      zf = 0.25*(z(inp)+z(ins)+z(inb)+z(inbs))

!.....Coordinates of point e'
      xi=xc(inp)*fxp+xc(ine)*fxe
      yi=yc(inp)*fxp+yc(ine)*fxe
      zi=zc(inp)*fxp+zc(ine)*fxe

!.....Interpolate gradients defined at CV centers to faces
      dfixi = gradfi(1,inp)*fxp+gradfi(1,ine)*fxe
      dfiyi = gradfi(2,inp)*fxp+gradfi(2,ine)*fxe
      dfizi = gradfi(3,inp)*fxp+gradfi(3,ine)*fxe


!.....The cell face interpolated gradient (d phi / dx_i)_j:
      dfixii = dfixi*d1x + arx/vole*( fi(ine)-fi(inp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
      dfiyii = dfiyi*d1y + ary/vole*( fi(ine)-fi(inp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
      dfizii = dfizi*d1z + arz/vole*( fi(ine)-fi(inp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Explicit diffusion
      fdfie = game*(dfixii*arx + dfiyii*ary + dfizii*arz)   
      ! Implicit diffussion 
      fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!.....Difusion coefficient
      de=game*are/dpn

!.....Convection fluxes - uds
      fm=fcf(inp)
      ce=min(fm,zero) 
      cp=max(fm,zero)

!.....System matrix coeffocients
      acfe(inp)=de-ce
      acfw(ine)=de+cp
!
!---------------------------------------------
!     CONVECTIVE FLUXES [ CENTRAL DIFFERENCING SCHEME (CDS) ]
!---------------------------------------------
      if (lcds) then
!.....Interpolate variable FI defined at CV centers to face using corrected CDS:
!         |__________Ue'_________|_______________Ucorr______________________|
      fii=fi(inp)*fxp+fi(ine)*fxe+(dfixi*(xf-xi)+dfiyi*(yf-yi)+dfizi*(zf-zi))
      !FII = face_interpolated(FI,gradFI,inp,idew,idns,idtb,fxp,fxe)
      fcfie=fm*fii
      else
!---------------------------------------------
!     Darwish-Moukalled TVD schemes for unstructured girds, IJHMT, 2003. 
!---------------------------------------------
!     Find r's - the gradient ratio. This is universal for all schemes.
!     If flow goes from P to E
     r1 = (2*gradFI(1,inp)*xpn + 2*gradFI(2,inp)*ypn + 2*gradFI(3,inp)*zpn)/(FI(INE)-FI(INP)) - 1.0d0  
!     If flow goes from E to P
     r2 = (2*gradFI(1,ine)*xpn + 2*gradFI(2,ine)*ypn + 2*gradFI(3,ine)*zpn)/(FI(INP)-FI(INE)) - 1.0d0 
!     Find Psi for [ MUSCL ] :
     psiw = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
     psie = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
!     High order flux at cell face
     fcfie =  ce*(fi(ine) + fxe*psie*(fi(inp)-fi(ine)))+ &
              cp*(fi(inp) + fxp*psiw*(fi(ine)-fi(inp)))
      endif

!.....First order upwind part
      fcfii=ce*fi(ine)+cp*fi(inp)
!.....Blend high order explicit and first order upwind 
      ffic = gam*(fcfie-fcfii)
!-------------------------------------------------------
!.....Explicit part of fluxes
!-------------------------------------------------------
      suadd=-ffic+fdfie-fdfii 
      ! sv(inp)=sv(inp)+suadd
      su(inp)=su(inp)+suadd
!-------------------------------------------------------
      bp(ine)=suadd

      end do !j-loop
      end do !i-loop
      end do !k-loop

      do k=2,nkmm
      do i=2,nimm
      do j=2,njmm
      inp=lk(k)+li(i)+j+idew

       su(inp)=su(inp)-bp(inp)
       ! sv(inp)=sv(inp)-bp(inp)

      end do !j-loop
      end do !i-loop
      end do !k-loop

      bp = 0.0d0
      
      return
      end
