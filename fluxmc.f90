!***********************************************************************
!
      subroutine fluxmc(inp,idew,idns,idtb, &
                        fif,fmcor)
!
!***********************************************************************
!
!     this routine calculates mass flux correction in the
!     second pressure-correction step which accounts for the
!     effects of non-orthogonality
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use coef
      use geometry
      use variables
      use gradients

      implicit none
!
!***********************************************************************
!
      integer, intent(in) :: inp,idew,idns,idtb
      real(prec), dimension(nxyza), intent(in) :: fif
      real(prec), intent(inout) :: fmcor
!
!     local variables
!
      integer :: ine,ins,inb,inbs
      real(prec) :: fxe,fxp
      real(prec) :: arx,ary,arz, &
                    xpn,ypn,zpn
      real(prec) :: rapr
      real(prec) :: are
      real(prec) :: nxx,nyy,nzz
      real(prec) :: xf,yf,zf
      real(prec) :: xep,yep,zep,xpp,ypp,zpp
      real(prec) :: dppnnr

      ine=inp+idew ! neighbour cv

      ins=inp-idns
      inb=inp-idtb
      inbs=inb-idns

      fxe=fif(inp)
      fxp=1.-fxe

!
!.....geometrical quantities

!.....precomputed face areas
      if(idew.eq.nj) then
        arx=ar1x(inp)
        ary=ar1y(inp)
        arz=ar1z(inp)
      elseif(idew.eq.1) then
        arx=ar2x(inp)
        ary=ar2y(inp)
        arz=ar2z(inp)
      else
        arx=ar3x(inp)
        ary=ar3y(inp)
        arz=ar3z(inp)
      endif

!.....distance vector components
      xpn=xc(ine)-xc(inp)
      ypn=yc(ine)-yc(inp)
      zpn=zc(ine)-zc(inp)

!.....surface vector magnitude squared
      are=sqrt(arx**2+ary**2+arz**2)

!.....unit normal vector
      nxx = arx/are
      nyy = ary/are
      nzz = arz/are

!.....distance from p' to n'-reciprocal value
      dppnnr = 1./((xpn*nxx)+(ypn*nyy)+(zpn*nzz))

!.....coordinates of cell face center 'f'
      xf = 0.25*(x(inp)+x(ins)+x(inb)+x(inbs))
      yf = 0.25*(y(inp)+y(ins)+y(inb)+y(inbs))
      zf = 0.25*(z(inp)+z(ins)+z(inb)+z(inbs))

!.....values at points p' and e' due to non-orthogonality. 
      xpp=xf-(xf-xc(inp))*nxx; ypp=yf-(yf-yc(inp))*nyy; zpp=zf-(zf-zc(inp))*nzz
      xep=xf-(xf-xc(ine))*nxx; yep=yf-(yf-yc(ine))*nyy; zep=zf-(zf-zc(ine))*nzz
!.....distances |p'p| and |e'e| projected ionto x,y,z-axis
      xpp=xpp-xc(inp); ypp=ypp-yc(inp); zpp=zpp-zc(inp)
      xep=xep-xc(ine); yep=yep-yc(ine); zep=zep-zc(ine)

!.....apu==1./ap x density - interpolated            
      rapr = (apu(inp)*den(inp)*vol(inp)*fxp+apu(ine)*den(ine)*vol(ine)*fxe)


!.....mass flux correction for the second p'-equation (source term)
      fmcor = rapr*are*((gradp(1,ine)*xep-gradp(1,inp)*xpp)   & 
                       +(gradp(2,ine)*yep-gradp(2,inp)*ypp)   &
                       +(gradp(3,ine)*zep-gradp(3,inp)*zpp))  &
                       *dppnnr 

      return
      end
