!***********************************************************************
!
      subroutine fluxmass_plain(idew,idns,idtb, &
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
      use variables
      use gradients

      implicit none
!
!***********************************************************************
!
      integer, intent(in) :: idew, idns, idtb
      real(prec), dimension(nxyza) :: arvx,arvy,arvz,fif,fcf
!
!     Local variables
!

      integer :: i, j, k, inp, & 
                 ine, ins, inb, inbs
      real(prec) :: fxe, fxp
      real(prec) :: arx, ary, arz
      real(prec) :: xpn,ypn,zpn
      real(prec) :: dene
      real(prec) :: xf,yf,zf,xi,yi,zi
      real(prec) :: ui,vi,wi
      real(prec) :: duxi,duyi,duzi
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
      fxp=1.0d0-fxe

!.....Geometrical quantities

!.....Coordinates of cell-face center
      xf = 0.25d0*(x(inp)+x(ins)+x(inb)+x(inbs))
      yf = 0.25d0*(y(inp)+y(ins)+y(inb)+y(inbs))
      zf = 0.25d0*(z(inp)+z(ins)+z(inb)+z(inbs))

!.....Coordinates of point e'
      xi=xc(inp)*fxp+xc(ine)*fxe
      yi=yc(inp)*fxp+yc(ine)*fxe
      zi=zc(inp)*fxp+zc(ine)*fxe

!.....Distance between p and e
      xpn=xc(ine)-xc(inp)
      ypn=yc(ine)-yc(inp)
      zpn=zc(ine)-zc(inp)

!.....Precomputed face areas
      arx=arvx(inp)
      ary=arvy(inp)
      arz=arvz(inp)

!.....end: geometrical quantities

!.....Density at the cell face
      dene=den(inp)*fxp+den(ine)*fxe
!
!.....cell face pressure gradients and velocities
!
!+++++interpolate velocities to face center:+++++++++++++++++++++++++
!.....interpolate gradients defined at cv centers to faces
      duxi = gradu(1,inp)*fxp+gradu(1,ine)*fxe
      duyi = gradu(2,inp)*fxp+gradu(2,ine)*fxe
      duzi = gradu(3,inp)*fxp+gradu(3,ine)*fxe
!        |________ue'_________|_______________ucorr_________________|
      ui=u(inp)*fxp+u(ine)*fxe!+duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)
!
      duxi = gradv(1,inp)*fxp+gradv(1,ine)*fxe
      duyi = gradv(2,inp)*fxp+gradv(2,ine)*fxe
      duzi = gradv(3,inp)*fxp+gradv(3,ine)*fxe
!        |________ve'_________|_______________vcorr_________________|
      vi=v(inp)*fxp+v(ine)*fxe!+duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)
!
!
      duxi = gradw(1,inp)*fxp+gradw(1,ine)*fxe
      duyi = gradw(2,inp)*fxp+gradw(2,ine)*fxe
      duzi = gradw(3,inp)*fxp+gradw(3,ine)*fxe
!        |________we'_________|_______________wcorr_________________|
      wi=w(inp)*fxp+w(ine)*fxe!+duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)

!
!.....Mass flux

!     // calculate the fluxes by dotting the interpolated velocity (to cell faces) with face normals
!     phi = (fvc::interpolate(u) & mesh.sf()) 

      fcf(inp)=dene*(ui*arx+vi*ary+wi*arz)

      end do
      end do
      end do

      return
      end
