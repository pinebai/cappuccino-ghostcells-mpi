!***********************************************************************
!
      subroutine piso_assemble_pressure_eq(idew,idns,idtb, &
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

      implicit none
!
!***********************************************************************
!
      integer, intent(in) :: idew, idns, idtb
      real(prec), dimension(nxyza) :: arvx,arvy,arvz,fif,acfe,acfw,fcf
!
!     local variables
!
      integer :: i, j, k, inp, &
                 ine, ins, inb, inbs
      real(prec) :: fxe, fxp
      real(prec) :: arx, ary, arz, are, smdpn
      real(prec) :: xpn,ypn,zpn,dene,nxx,nyy,nzz,sfdpnr
      real(prec) :: xf,yf,zf,xi,yi,zi
      real(prec) :: ui,vi,wi
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
      fxp=1.-fxe


!.....geometrical quantities...............................
!.....coordinates of cell-face center
      xf = 0.25*(x(inp)+x(ins)+x(inb)+x(inbs))
      yf = 0.25*(y(inp)+y(ins)+y(inb)+y(inbs))
      zf = 0.25*(z(inp)+z(ins)+z(inb)+z(inbs))

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

      sfdpnr=1./(arx*xpn*nxx+ary*ypn*nyy+arz*zpn*nzz)
!.....end: geometrical quantities...............................

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
!+++++interpolate velocities to face center:+++++++++++++++++++++++++
!.....interpolate gradients defined at cv centers to faces
      duxi = gradu(1,inp)*fxp+gradu(1,ine)*fxe
      duyi = gradu(2,inp)*fxp+gradu(2,ine)*fxe
      duzi = gradu(3,inp)*fxp+gradu(3,ine)*fxe
!        |________ue'_________|_______________ucorr_________________|
      ui=u(inp)*fxp+u(ine)*fxe+duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)
!
      duxi = gradv(1,inp)*fxp+gradv(1,ine)*fxe
      duyi = gradv(2,inp)*fxp+gradv(2,ine)*fxe
      duzi = gradv(3,inp)*fxp+gradv(3,ine)*fxe
!        |________ve'_________|_______________vcorr_________________|
      vi=v(inp)*fxp+v(ine)*fxe+duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)
!
!
      duxi = gradw(1,inp)*fxp+gradw(1,ine)*fxe
      duyi = gradw(2,inp)*fxp+gradw(2,ine)*fxe
      duzi = gradw(3,inp)*fxp+gradw(3,ine)*fxe
!        |________we'_________|_______________wcorr_________________|
      wi=w(inp)*fxp+w(ine)*fxe+duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)
!+++++end: interpolate velocities to face center:+++++++++++++++++++++++++

!
!.....mass flux
!// calculate the fluxes by dotting the interpolated velocity (to cell faces) with face normals
!     phi = (fvc::interpolate(u) & mesh.sf()) 
      fcf(inp)=dene*(ui*arx+vi*ary+wi*arz)

      end do
      end do
      end do

      return
      end
