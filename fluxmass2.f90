!***********************************************************************
!
      subroutine fluxmass(idew,idns,idtb, &
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
      use variables
      use gradients
      use fieldmanipulation

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
                    xpn,ypn,zpn,dene,smdpn,sfdpnr
      real(prec) :: xf,yf,zf,xi,yi,zi
      real(prec) :: nxx,nyy,nzz,xpp,ypp,zpp,xep,yep,zep
      real(prec) :: ui,vi,wi,ue,ve,we,dpe,dpex,dpey,dpez
      real(prec) :: dpxi,dpyi,dpzi
      real(prec) :: duxi,duyi,duzi

!
!.....calculate east cell face
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
      fxp=1.0d0-fxe

!.....geometrical quantities

!.....coordinates of cell-face center
      xf = 0.25d0*(x(inp)+x(ins)+x(inb)+x(inbs))
      yf = 0.25d0*(y(inp)+y(ins)+y(inb)+y(inbs))
      zf = 0.25d0*(z(inp)+z(ins)+z(inb)+z(inbs))

!.....coordinates of point e'
      xi=xc(inp)*fxp+xc(ine)*fxe
      yi=yc(inp)*fxp+yc(ine)*fxe
      zi=zc(inp)*fxp+zc(ine)*fxe

!.....distance between p and e
      xpn=xc(ine)-xc(inp)
      ypn=yc(ine)-yc(inp)
      zpn=zc(ine)-zc(inp)

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

!.....volume
      sfdpnr=1.d0/(arx*xpn*nxx+ary*ypn*nyy+arz*zpn*nzz)

!.....end: geometrical quantities

!.....density at the cell face
      dene=den(inp)*fxp+den(ine)*fxe
!
!.....coefficients of pressure-correction equation
      smdpn=(arx*arx+ary*ary+arz*arz)*sfdpnr
      acfe(inp)=dene*(fxp*vol(inp)*apu(inp)+fxe*vol(ine)*apu(ine))*smdpn
      acfw(ine)=acfe(inp) 

!
!.....cell face pressure gradients and velocities
!
!////////////////////////////////////////////////////////
!     Rhie-Chow velcity interolation at face
!
!   dpcv <- gradp(1/2/3,nxzya)
!
!   uf=ui+dpdxi-api*sf*(pn-pp)
!         _
!   ui-> (u)f -> second order interpolation at face
!            _______________ 
!   dpdxi-> (dpdx*vol*(1/ap))f -> second order interpolation at cell face
!          ______
!   api*sf*(pn-pp) -> (1/ap)f*sf*(pn-pp) cell face coefficient ap x area_f x (p1-p2)
!  
!   finally:     
!         __     _______________     ______
!   uf = (u)f + (dpdx*vol*(1/ap))f - (1/ap)f*sf*(pn-pp)
!
!   and:
!   flmass = densit*dot(uf,sf)
!/////////////////////////////////////////////////////////

!+++++interpolate velocities to face center:+++++++++++++++++++++++++
!.....interpolate gradients defined at cv centers to faces
      duxi = gradu(1,inp)*fxp+gradu(1,ine)*fxe
      duyi = gradu(2,inp)*fxp+gradu(2,ine)*fxe
      duzi = gradu(3,inp)*fxp+gradu(3,ine)*fxe
!        |________ue'_________|_______________ucorr___________________|
      ui=u(inp)*fxp+u(ine)*fxe+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))
!      ui = face_interpolated(u,gradu,inp,idew,idns,idtb,fxp,fxe)

      duxi = gradv(1,inp)*fxp+gradv(1,ine)*fxe
      duyi = gradv(2,inp)*fxp+gradv(2,ine)*fxe
      duzi = gradv(3,inp)*fxp+gradv(3,ine)*fxe
!        |________ve'_________|_______________vcorr___________________|
      vi=v(inp)*fxp+v(ine)*fxe+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))
!      vi = face_interpolated(v,gradv,inp,idew,idns,idtb,fxp,fxe)

      duxi = gradw(1,inp)*fxp+gradw(1,ine)*fxe
      duyi = gradw(2,inp)*fxp+gradw(2,ine)*fxe
      duzi = gradw(3,inp)*fxp+gradw(3,ine)*fxe
!        |________we'_________|_______________wcorr___________________|
      wi=w(inp)*fxp+w(ine)*fxe+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)) 
!      wi = face_interpolated(w,gradw,inp,idew,idns,idtb,fxp,fxe) 


!+++++interpolate pressure gradients to cell face center++++++++++++++++++
      dpxi=(fxe*vol(ine)*apu(ine)+fxp*vol(inp)*apu(inp))*(fxe*gradp(1,ine)+fxp*gradp(1,inp))*xpn*nxx
      dpyi=(fxe*vol(ine)*apv(ine)+fxp*vol(inp)*apv(inp))*(fxe*gradp(2,ine)+fxp*gradp(2,inp))*ypn*nyy
      dpzi=(fxe*vol(ine)*apw(ine)+fxp*vol(inp)*apw(inp))*(fxe*gradp(3,ine)+fxp*gradp(3,inp))*zpn*nzz


!+++++pressure deriv. along normal+++++++++++++++++++++++++++++++++++++++++ 
!.....values at points p' and e' due to non-orthogonality. 
      xpp=xf-(xf-xc(inp))*nxx
      ypp=yf-(yf-yc(inp))*nyy
      zpp=zf-(zf-zc(inp))*nzz

      xep=xf-(xf-xc(ine))*nxx
      yep=yf-(yf-yc(ine))*nyy
      zep=zf-(zf-zc(ine))*nzz

!.....distances |p'p| and |e'e| projected ionto x,y,z-axis
      xpp=xpp-xc(inp)
      ypp=ypp-yc(inp)
      zpp=zpp-zc(inp)

      xep=xep-xc(ine)
      yep=yep-yc(ine)
      zep=zep-zc(ine)

      dpe = (p(ine)-p(inp)) + &
      ( gradp(1,ine)*xep+gradp(2,ine)*yep+gradp(3,ine)*zep - & !<<--correction
        gradp(1,inp)*xpp+gradp(2,inp)*ypp+gradp(3,inp)*zpp  )  !<<|

      dpex = (fxe*vol(ine)*apu(ine)+fxp*vol(inp)*apu(inp))*(arx*sfdpnr)*dpe
      dpey = (fxe*vol(ine)*apv(ine)+fxp*vol(inp)*apv(inp))*(ary*sfdpnr)*dpe
      dpez = (fxe*vol(ine)*apw(ine)+fxp*vol(inp)*apw(inp))*(arz*sfdpnr)*dpe


!+++++rhie-chow interpolation++++++++++++++++++++++++++++++++++++++++++++++ 
      ue = ui - dpex + dpxi
      ve = vi - dpey + dpyi
      we = wi - dpez + dpzi


!
!.....mass flux via rhie-chow interpolation of velocity
      fcf(inp)=dene*(ue*arx+ve*ary+we*arz)

      end do
      end do
      end do

      return
      end
