!***********************************************************************
!
  subroutine piso_assemble_pressure_eq(idew,idns,idtb, &
                                       arvx,arvy,arvz, &
                                       fif,acfe,acfw,fcf,ipc,bp)
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use coef
  !use coefb
  use variables
  use gradients

  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: idew, idns, idtb,ipc
  real(prec), dimension(nxyza) :: arvx,arvy,arvz,fif,acfe,acfw,fcf,bp
!
!     Local variables
!

  integer :: i, j, k, inp, &
             ine, ins, inb, inbs
  real(prec) :: fxe, fxp
  real(prec) :: arx, ary, arz, are!, smdpn
  real(prec) :: xpn,ypn,zpn,dene,nxx,nyy,nzz!,sfdpnr
  real(prec) :: xf,yf,zf,xi,yi,zi
  real(prec) :: ui,vi,wi
  real(prec) :: duxi,duyi,duzi

  real(prec) :: ixi1,ixi2,ixi3,dpn,costheta,costn
  real(prec) :: fdfie,fdfii
  real(prec) :: d2x,d2y,d2z,d1x,d1y,d1z
  real(prec) :: dfixi,dfiyi,dfizi
  real(prec) :: dfixii,dfiyii,dfizii
  real(prec) :: suadd,vole,game

  ! zero the bp array
  bp = 0.0_dp

  ! Lopp inner faces on east, nort or top
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

  ! geometrical quantities:

  ! coordinates of cell-face center
  xf = 0.25d0*(x(inp)+x(ins)+x(inb)+x(inbs))
  yf = 0.25d0*(y(inp)+y(ins)+y(inb)+y(inbs))
  zf = 0.25d0*(z(inp)+z(ins)+z(inb)+z(inbs))

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
  costn = costheta
  ! orthogonal correction: nrelax =  0 : 
  ! costn = 1.0d0
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

!..


! !.....precomputed face areas
!       arx=arvx(inp)
!       ary=arvy(inp)
!       arz=arvz(inp)

! !.....cell face area
!       are=sqrt(arx**2+ary**2+arz**2)

! !.....unit vectors of the normal
!       nxx=arx/are
!       nyy=ary/are
!       nzz=arz/are

!       sfdpnr=1./(arx*xpn*nxx+ary*ypn*nyy+arz*zpn*nzz)
! !.....end: geometrical quantities...............................

! !.....density at the cell face
!       dene=den(inp)*fxp+den(ine)*fxe
! !
! !.....coefficients of pressure-correction equation
!       smdpn=(arx*arx+ary*ary+arz*arz)*sfdpnr
!       acfe(inp)=dene*(fxp*vol(inp)*apu(inp)+fxe*vol(ine)*apu(ine))*smdpn
!       acfw(ine)=acfe(inp) 

!..
  ! diffusion coef
  dene=den(inp)*fxp+den(ine)*fxe
  game = dene*(fxp*vol(inp)*apu(inp)+fxe*vol(ine)*apu(ine))
  ! system matrix coefficients
  acfe(inp)=game*are/dpn
  acfw(ine)=acfe(inp) 

  !..
  if(ipc.gt.0) then ! on second and higher non-orthogonal correction for pressure
  !..
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
  suadd = fdfie-fdfii 

  ! add source to this cell
  su(inp) = su(inp)+suadd
  ! su(ine) = su(ine)-suadd

  ! store source
  bp(inp) = suadd

  !..
  !else ! mass flux
  !..

  !.....cell face pressure gradients and velocities

  !+++++interpolate velocities to face center:+++++++++++++++++++++++++
  !.....interpolate gradients defined at cv centers to faces
  duxi = gradu(1,inp)*fxp+gradu(1,ine)*fxe
  duyi = gradu(2,inp)*fxp+gradu(2,ine)*fxe
  duzi = gradu(3,inp)*fxp+gradu(3,ine)*fxe
  !  |________ue'_________|_______________ucorr_________________|
  ui=u(inp)*fxp+u(ine)*fxe+duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)
  !
  duxi = gradv(1,inp)*fxp+gradv(1,ine)*fxe
  duyi = gradv(2,inp)*fxp+gradv(2,ine)*fxe
  duzi = gradv(3,inp)*fxp+gradv(3,ine)*fxe
  !  |________ve'_________|_______________vcorr_________________|
  vi=v(inp)*fxp+v(ine)*fxe+duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)
  !
  !
  duxi = gradw(1,inp)*fxp+gradw(1,ine)*fxe
  duyi = gradw(2,inp)*fxp+gradw(2,ine)*fxe
  duzi = gradw(3,inp)*fxp+gradw(3,ine)*fxe
  !  |________we'_________|_______________wcorr_________________|
  wi=w(inp)*fxp+w(ine)*fxe+duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)
  !+++++end: interpolate velocities to face center:+++++++++++++++++++++++++

  !
  !.....mass flux
  !// calculate the fluxes by dotting the interpolated velocity (to cell faces) with face normals
  !     phi = (fvc::interpolate(u) & mesh.sf()) 

  fcf(inp)=dene*(ui*arx+vi*ary+wi*arz)

  !..
  endif ! ipc condition
  !..

  end do
  end do
  end do

  ! print*,ipc,':',sum(abs(bp))

  if(ipc.gt.0) then

    ! now add source to neighbour cells ine
    do k=3,nkmm
    do i=3,nimm
    do j=3,njmm
    inp=lk(k)+li(i)+j+idew

     su(inp)=su(inp)-bp(inp)

    end do !j-loop
    end do !i-loop
    end do !k-loop

    ! bp = 0.0d0

  endif

  return
  end

! Trenutna situacija je (dok je ipc.gt.0) da u prvoj piso korekciji (icorr=1 u piso_multiple_correction.f90)
! nema neortogonalne delove izvornog clana jer su gradijet pritiska i pritiska
! inicijalizovani nulom, dok u drugoj, trecoj, piso korekciji ima.
