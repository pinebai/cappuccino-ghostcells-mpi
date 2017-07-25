!***********************************************************************
!
      subroutine fluxmass(idew,idns,idtb, &
                          arvx,arvy,arvz, &
                          fif,acfe,acfw,fcf)
!
!***********************************************************************
!     Mass flux at cell face via a Rhie-Chow face interpolated velocity
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use geometry
      use coef
      use variables
      use gradients
      use fieldManipulation

      implicit none
!
!***********************************************************************
!
      integer, intent(in) :: idew, idns, idtb
      real(prec), dimension(nxyza) :: arvx,arvy,arvz,fif,acfe,acfw,fcf


!
!     Local variables
!
      integer :: i, j, k, inp, &
                 ine, ins, inb, inbs
      real(prec) :: fxe, fxp, &
                    arx, ary, arz, are, &
                    due, dve, dwe,cap,xpn,ypn,zpn,vole,dene
      real(prec) :: xf,yf,zf,xi,yi,zi
      real(prec) :: nxx,nyy,nzz,xpp,ypp,zpp,xep,yep,zep,dep
      real(prec) :: ui,vi,wi,ue,ve,we,dpe
      real(prec) :: dpdxi,dpdyi,dpdzi
      !real(prec) :: duxi,duyi,duzi
 

!.....TENTATIVE(!) VELOCITY GRADIENTS: 
      if (lstsq) then 
        call grad_lsq_qr(u,gradu,2)
        call grad_lsq_qr(v,gradv,2)
        call grad_lsq_qr(w,gradw,2)
      elseif (gauss) then
        call grad_gauss(u,gradu(1,:),gradu(2,:),gradu(3,:))
        call grad_gauss(v,gradv(1,:),gradv(2,:),gradv(3,:))
        call grad_gauss(w,gradw(1,:),gradw(2,:),gradw(3,:))
      endif

!
!.....Calculate east cell face
!
      do k=2,nkmm
      do i=2,nimm
      do j=2,njmm

      inp=lk(k)+li(i)+j
      ine=inp+idew
      ins=inp-idns
      inb=inp-idtb
      inbs=inb-idns

      fxe=fif(inp)
      fxp=1.-fxe


!.....GEOMETRICAL QUANTITIES...............................
!.....Coordinates of cell-face center
      xf = 0.25*(x(inp)+x(ins)+x(inb)+x(inbs))
      yf = 0.25*(y(inp)+y(ins)+y(inb)+y(inbs))
      zf = 0.25*(z(inp)+z(ins)+z(inb)+z(inbs))

!.....Coordinates of point e'
      xi=xc(inp)*fxp+xc(ine)*fxe
      yi=yc(inp)*fxp+yc(ine)*fxe
      zi=zc(inp)*fxp+zc(ine)*fxe

!.....Distance between P and E
      xpn=xc(ine)-xc(inp)
      ypn=yc(ine)-yc(inp)
      zpn=zc(ine)-zc(inp)

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

!.....Volume
      vole=arx*xpn+ary*ypn+arz*zpn
!.....END: GEOMETRICAL QUANTITIES...............................

!.....Cell face coefficients  1/ap(inp)
      due=-(apu(inp)*fxp+apu(ine)*fxe)
      dve=-(apv(inp)*fxp+apv(ine)*fxe)
      dwe=-(apw(inp)*fxp+apw(ine)*fxe)
!.....Density at the cell face
      dene=den(inp)*fxp+den(ine)*fxe
!
!.....Coefficients of pressure-correction equation
      cap=-dene*(due*arx**2+dve*ary**2+dwe*arz**2)
      acfe(inp)=cap
      acfw(ine)=cap
!
!.....CELL FACE PRESSURE GRADIENTS AND VELOCITIES
!
!////////////////////////////////////////////////////////
!     RHIE-CHOW velcity interolation at face
!
!   DPCV <- gradP(1/2/3,NXZYA)
!
!   Uf=UI+DPDXI-API*Sf*(Pn-Pp)
!         _
!   UI-> (U)f -> second order interpolation at face
!            _______________ 
!   DPDXI-> (dPdx*Vol*(1/ap))f -> second order interpolation at cell face
!          ______
!   API*Sf*(Pn-Pp) -> (1/ap)f*Sf*(Pn-Pp) cell face coefficient Ap x Area_f x (p1-p2)
!  
!   Finally:     
!         __     _______________     ______
!   Uf = (U)f + (dPdx*Vol*(1/ap))f - (1/ap)f*Sf*(Pn-Pp)
!
!   and:
!   Flmass = Densit*dot(Uf,Sf)
!/////////////////////////////////////////////////////////

!+++++Interpolate velocities to face center:+++++++++++++++++++++++++
!.....Interpolate gradients defined at CV centers to faces
!      duxi = gradu(1,inp)*fxp+gradu(1,ine)*fxe
!      duyi = gradu(2,inp)*fxp+gradu(2,ine)*fxe
!      duzi = gradu(3,inp)*fxp+gradu(3,ine)*fxe
!        |________Ue'_________|_______________Ucorr___________________|
!      ui=u(inp)*fxp+u(ine)*fxe+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))
      ui = face_interpolated(u,gradu,inp,idew,idns,idtb,fxp,fxe)

!      duxi = gradv(1,inp)*fxp+gradv(1,ine)*fxe
!      duyi = gradv(2,inp)*fxp+gradv(2,ine)*fxe
!      duzi = gradv(3,inp)*fxp+gradv(3,ine)*fxe
!        |________Ve'_________|_______________Vcorr___________________|
!      vi=v(inp)*fxp+v(ine)*fxe+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))
      vi = face_interpolated(v,gradv,inp,idew,idns,idtb,fxp,fxe)

!      duxi = gradw(1,inp)*fxp+gradw(1,ine)*fxe
!      duyi = gradw(2,inp)*fxp+gradw(2,ine)*fxe
!      duzi = gradw(3,inp)*fxp+gradw(3,ine)*fxe
!        |________We'_________|_______________Wcorr___________________|
!      wi=w(inp)*fxp+w(ine)*fxe+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)) 
      wi = face_interpolated(w,gradw,inp,idew,idns,idtb,fxp,fxe) 
  

!+++++Interpolate pressure gradients to cell face center++++++++++++++++++
      dpdxi=((apu(ine)*gradp(1,ine))*fxe+ &
             (apu(inp)*gradp(1,inp))*fxp)*vole
      dpdyi=((apv(ine)*gradp(2,ine))*fxe+ &
             (apv(inp)*gradp(2,inp))*fxp)*vole
      dpdzi=((apw(ine)*gradp(3,ine))*fxe+ &
             (apw(inp)*gradp(3,inp))*fxp)*vole

!+++++Pressure deriv. along normal+++++++++++++++++++++++++++++++++++++++++ 
      dep=xpn*nxx+ypn*nyy+zpn*nzz
!.....Values at points p' and e' due to non-orthogonality. 
      xpp=xf-(xf-xc(inp))*nxx; ypp=yf-(yf-yc(inp))*nyy; zpp=zf-(zf-zc(inp))*nzz
      xep=xf-(xf-xc(ine))*nxx; yep=yf-(yf-yc(ine))*nyy; zep=zf-(zf-zc(ine))*nzz
!.....Distances |P'P| and |E'E| projected ionto x,y,z-axis
      xpp=xpp-xc(inp); ypp=ypp-yc(inp); zpp=zpp-zc(inp)
      xep=xep-xc(ine); yep=yep-yc(ine); zep=zep-zc(ine)

      !# <-nonortho. correction; uncomment
      dpe = (p(ine)-p(inp))/(dep+small)  &                                   
       + ( gradp(1,ine)*xep+gradp(2,ine)*yep+gradp(3,ine)*zep - &            !#
           gradp(1,inp)*xpp+gradp(2,inp)*ypp+gradp(3,inp)*zpp  )/(dep+small) !#


!+++++Rhie-Chow Interpolation++++++++++++++++++++++++++++++++++++++++++++++ 
      ue = ui + due*dpe*vole + dpdxi
      ve = vi + dve*dpe*vole + dpdyi
      we = wi + dwe*dpe*vole + dpdzi

!
!.....Mass flux
      fcf(inp)=dene*(ue*arx+ve*ary+we*arz)

      end do
      end do
      end do

      return
      end
