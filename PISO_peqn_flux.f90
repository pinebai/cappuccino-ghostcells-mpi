!***********************************************************************
!
  subroutine piso_correct_flmass_peqn_flux(fcf,arvx,arvy,arvz,idew,idns,idtb,acfe,fif)
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
  use gradients

  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: idew, idns, idtb
  real(prec), dimension(nxyza) :: arvx,arvy,arvz,fif,acfe,fcf

!
!     Local variables
!

  integer :: i, j, k, inp, &
             ine, ins, inb,  inbs!, inn!, int
  ! integer:: idew, idns, idtb
  real(prec) :: fxe, fxp
  real(prec) :: arx, ary, arz, are
  real(prec) :: xpn,ypn,zpn,dene,nxx,nyy,nzz
  real(prec) :: xf,yf,zf,xi,yi,zi

  real(prec) :: ixi1,ixi2,ixi3,dpn,costheta,costn
  real(prec) :: fdfie,fdfii
  real(prec) :: d2x,d2y,d2z,d1x,d1y,d1z
  real(prec) :: dfixi,dfiyi,dfizi
  real(prec) :: dfixii,dfiyii,dfizii
  real(prec) :: suadd,vole,game


  !  INNER FACES


  ! Lopp inner faces on east, nort or top
  do k=2,nkm
  do i=2,nim
  do j=2,njm

  inp=lk(k)+li(i)+j
  ine=inp+idew
  ins=inp-idns
  inb=inp-idtb   
  inbs=inb-idns

  fxe=fif(inp)
  fxp=1.-fxe

  ! geometrical quantities:

  ! coordinates of cell-face center
  xf = 0.25*(x(inp)+x(ins)+x(inb)+x(inbs))
  yf = 0.25*(y(inp)+y(ins)+y(inb)+y(inbs))
  zf = 0.25*(z(inp)+z(ins)+z(inb)+z(inbs))

  ! coordinates of point e'
  xi=xc(inp)*fxp+xc(ine)*fxe
  yi=yc(inp)*fxp+yc(ine)*fxe
  zi=zc(inp)*fxp+zc(ine)*fxe

  ! distance between p and e
  xpn=xc(ine)-xc(inp)
  ypn=yc(ine)-yc(inp)
  zpn=zc(ine)-zc(inp)

!..

  ! distance from p to neighbor n
  dpn=sqrt(xpn**2+ypn**2+zpn**2)  

  ! components of the unit vector i_ksi
  ixi1=xpn/dpn
  ixi2=ypn/dpn
  ixi3=zpn/dpn

  ! precomputed face areas
  arx=arvx(inp)
  ary=arvy(inp)
  arz=arvz(inp)

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! ___ __
  ! dpn.sf
  vole=xpn*arx+ypn*ary+zpn*arz

  ! unit vectors of the normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are

  ! angle between vectorsa n and i_xi - we need cosine
  costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

  ! relaxation factor for higher-order cell face gradient
  ! minimal correction: nrelax = +1 :
  !costn = costheta
  ! orthogonal correction: nrelax =  0 : 
  costn = 1.0d0
  ! over-relaxed approach: nrelax = -1 :
  !costn = 1./costheta
  ! in general, nrelax can be any signed integer from some 
  ! reasonable interval [-nrelax,nrelax] (or even real number): 
  !costn = costheta**nrelax

  ! overrelaxed correction vector d2, where s=dpn+d2
  d1x = costn
  d1y = costn
  d1z = costn
  !
  d2x = xpn*costn
  d2y = ypn*costn
  d2z = zpn*costn

  ! diffusion coef
  dene=den(inp)*fxp+den(ine)*fxe
  game = dene*(fxp*vol(inp)*apu(inp)+fxe*vol(ine)*apu(ine))

  ! interpolate gradients defined at cv centers to faces
  dfixi = gradp(1,inp)*fxp+gradp(1,ine)*fxe
  dfiyi = gradp(2,inp)*fxp+gradp(2,ine)*fxe
  dfizi = gradp(3,inp)*fxp+gradp(3,ine)*fxe


  ! the cell face interpolated gradient (d phi / dx_i)_j:
  dfixii = dfixi*d1x + arx/vole*( p(ine)-p(inp) - (dfixi*d2x + dfiyi*d2y + dfizi*d2z) ) 
  dfiyii = dfiyi*d1y + ary/vole*( p(ine)-p(inp) - (dfixi*d2x + dfiyi*d2y + dfizi*d2z) ) 
  dfizii = dfizi*d1z + arz/vole*( p(ine)-p(inp) - (dfixi*d2x + dfiyi*d2y + dfizi*d2z) ) 

  !-------------------------------------------------------------------------------------
  ! explicit difussion
  fdfie = game*(dfixii*arx + dfiyii*ary + dfizii*arz)   

  ! implicit difussion 
  fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)
  !-------------------------------------------------------------------------------------

  ! explicit part of diffusion fluxes
  suadd=fdfie-fdfii 

                                                                                    
  ! Correct mass fluxes at inner cv-faces only (only inner flux)                      
  fcf(inp)=fcf(inp)-acfe(inp)*(p(ine)-p(inp))-suadd  

  end do
  end do
  end do


  return
  end


