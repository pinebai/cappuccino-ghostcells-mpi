!***********************************************************************
!
  subroutine outbc
!
!***********************************************************************
!
! Scale mass flow in domain to satisfy global conservation!
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use coef
  use variables
  use bc
  use boundc
  use buoy
  use time_mod
  use obstacle

 implicit none
!
!***********************************************************************
!

  integer :: i, j, k, inp, inbc, inh
  real(prec) :: fac, flow

!.....check if cavity
  if(.not.lout) return
  flow=0.
!----------------------
!.....west
!----------------------
  if(low) then
    do k=2,nkm
    do j=2,njm
    inp=lk(k)+li(2)+j
    inbc=(j-1)*nk+k
    lbhlp=lbw(inbc+jks)
    if(lbhlp.eq.2) then
    inh=inp+nj
    f1(inp)=den(inp)*(ar1x(inp)*u(inh) &
                     +ar1y(inp)*v(inh) &
                     +ar1z(inp)*w(inh))
    flow=flow-f1(inp)
    endif
    enddo
    enddo
  endif
!------------------------
!.....east
!------------------------
  if(loe) then
    do k=2,nkm
    do j=2,njm
    inp=lk(k)+li(nim)+j
    inbc=(j-1)*nk+k
    lbhlp=lbe(inbc+jks)
    if(lbhlp.eq.2) then
    inh=inp-nj
    f1(inh)=den(inp)*(ar1x(inp)*u(inh) &
                     +ar1y(inp)*v(inh) &
                     +ar1z(inp)*w(inh))
    flow=flow+f1(inh)
    endif
    enddo
    enddo
  endif
!--------------------
!.....south
!--------------------
  if(los) then
    do k=2,nkm
    do i=2,nim
    inp=lk(k)+li(i)+2
    inbc=(i-1)*nk+k
    lbhlp=lbs(inbc+iks)
    if(lbhlp.eq.2) then
    inh=inp+1
    f2(inp)=den(inp)*(ar1x(inp)*u(inh) &
                     +ar1y(inp)*v(inh) &
                     +ar1z(inp)*w(inh))
    flow=flow-f2(inp)
    endif
    enddo
    enddo
  endif
!----------------------
!.....north
!----------------------
  if(lon) then
    do k=2,nkm
    do i=2,nim
    inp=lk(k)+li(i)+njm
    inbc=(i-1)*nk+k
    lbhlp=lbn(inbc+iks)
    if(lbhlp.eq.2) then
    inh=inp-1
    f2(inh)=den(inp)*(ar1x(inp)*u(inh) &
                     +ar1y(inp)*v(inh) &
                     +ar1z(inp)*w(inh))
    flow=flow+f2(inh)
    endif
    enddo
    enddo
  endif
!---------------------
!.....bottom
!---------------------
  if(lob) then
    do i=2,nim
    do j=2,njm
    inp=lk(2)+li(i)+j
    inbc=(i-1)*nj+j
    lbhlp=lbb(inbc+ijs)
    if(lbhlp.eq.2) then
    inh=inp+nij
    f3(inp)=den(inp)*(ar1x(inp)*u(inh) &
                     +ar1y(inp)*v(inh) &
                     +ar1z(inp)*w(inh))
    flow=flow-f3(inp)
    endif
    enddo
    enddo
  endif
!-----------------------------
!.....top
!-----------------------------
  if(lot) then
    do i=2,nim
    do j=2,njm
    inp=lk(nkm)+li(i)+j
    inbc=(i-1)*nj+j
    lbhlp=lbt(inbc+ijs)
    if(lbhlp.eq.2) then
    inh=inp-nij
    f3(inh)=den(inp)*(ar1x(inp)*u(inh) &
                     +ar1y(inp)*v(inh) &
                     +ar1z(inp)*w(inh))
    flow=flow+f3(inh)
    endif
    enddo
    enddo
  endif

! zero gradient b.c. employed at outlet
! calculate correction factor
  fac=flowin/(flow+small)
! write(*,'(2x,a5,1p10e15.6)') ' fac: ',fac,' flowin: ',flowin,' flow: ',flow

!
!.....west
  if(low) then
    do k=2,nkm
    do j=2,njm
    inp=lk(k)+li(2)+j
    inbc=(j-1)*nk+k
    lbhlp=lbw(inbc+jks)
    if(lbhlp.eq.2) then
    inh=inp+nj
    u(inp)=u(inh)*fac
    v(inp)=v(inh)*fac
    w(inp)=w(inh)*fac
    f1(inp)=f1(inp)*fac

    te(inp)=te(inh)
    ed(inp)=ed(inh)
    vis(inp)=vis(inh)
    t(inp)=t(inh)
    vart(inp)=vart(inh)
    con(inp)=con(inh)
    utt(inp)=utt(inh)
    vtt(inp)=vtt(inh)
    wtt(inp)=wtt(inh)
    endif
    enddo
    enddo
  endif
!.....east
  if(loe) then
    do k=2,nkm
    do j=2,njm
    inp=lk(k)+li(nim)+j
    inbc=(j-1)*nk+k
    lbhlp=lbe(inbc+jks)
    if(lbhlp.eq.2) then
    inh=inp-nj
    u(inp)=u(inh)*fac
    v(inp)=v(inh)*fac
    w(inp)=w(inh)*fac
    f1(inh)=f1(inh)*fac

    te(inp)=te(inh)
    ed(inp)=ed(inh)
    vis(inp)=vis(inh)
    t(inp)=t(inh)
    vart(inp)=vart(inh)
    con(inp)=con(inh)
    utt(inp)=utt(inh)
    vtt(inp)=vtt(inh)
    wtt(inp)=wtt(inh)
    endif
    enddo
    enddo
  endif
!.....south
  if(los) then
    do k=2,nkm
    do i=2,nim
    inp=lk(k)+li(i)+2
    inbc=(i-1)*nk+k
    lbhlp=lbs(inbc+iks)
    if(lbhlp.eq.2) then
    inh=inp+1
    u(inp)=u(inh)*fac
    v(inp)=v(inh)*fac
    w(inp)=w(inh)*fac
    f2(inp)=f2(inp)*fac
  
    te(inp)=te(inh)
    ed(inp)=ed(inh)
    t(inp)=t(inh)
    vis(inp)=vis(inh)
    vart(inp)=vart(inh)
    con(inp)=con(inh)
    utt(inp)=utt(inh)
    vtt(inp)=vtt(inh)
    wtt(inp)=wtt(inh)
    endif
    enddo
    enddo
  endif
!.....north
  if(lon) then
    do k=2,nkm
    do i=2,nim
    inp=lk(k)+li(i)+njm
    inbc=(i-1)*nk+k
    lbhlp=lbn(inbc+iks)
    if(lbhlp.eq.2) then
    inh=inp-1
    u(inp)=u(inh)*fac
    v(inp)=v(inh)*fac
    w(inp)=w(inh)*fac
    f2(inh)=f2(inh)*fac

    te(inp)=te(inh)
    ed(inp)=ed(inh)
    t(inp)=t(inh)
    vis(inp)=vis(inh)
    vart(inp)=vart(inh)
    con(inp)=con(inh)
    utt(inp)=utt(inh)
    vtt(inp)=vtt(inh)
    wtt(inp)=wtt(inh)
    endif
    enddo
    enddo
  endif
!.....bottom
  if(lob) then
    do i=2,nim
    do j=2,njm
    inp=lk(2)+li(i)+j
    inbc=(i-1)*nj+j
    lbhlp=lbb(inbc+ijs)
    if(lbhlp.eq.2) then
    inh=inp+nij
    u(inp)=u(inh)*fac
    v(inp)=v(inh)*fac
    w(inp)=w(inh)*fac
    f3(inp)=f3(inp)*fac

    te(inp)=te(inh)
    ed(inp)=ed(inh)
    t(inp)=t(inh)
    vis(inp)=vis(inh)
    vart(inp)=vart(inh)
    con(inp)=con(inh)
    utt(inp)=utt(inh)
    vtt(inp)=vtt(inh)
    wtt(inp)=wtt(inh)
    endif
    enddo
    enddo
  endif
!.....top
  if(lot) then
    do i=2,nim
    do j=2,njm
    inp=lk(nkm)+li(i)+j
    inbc=(i-1)*nj+j
    lbhlp=lbt(inbc+ijs)
    if(lbhlp.eq.2) then
    inh=inp-nij
    u(inp)=u(inh)*fac
    v(inp)=v(inh)*fac
    w(inp)=w(inh)*fac
    f3(inh)=f3(inh)*fac

    te(inp)=te(inh)
    ed(inp)=ed(inh)
    vis(inp)=vis(inh)
    t(inp)=t(inh)
    vart(inp)=vart(inh)
    con(inp)=con(inh)
    utt(inp)=utt(inh)
    vtt(inp)=vtt(inh)
    wtt(inp)=wtt(inh)
    endif
    enddo
    enddo
  endif

  return
  end
