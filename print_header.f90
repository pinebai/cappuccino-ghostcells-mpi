!***********************************************************************
!
subroutine print_header
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use boundc
  use title_mod
  use buoy
  use time_mod

  implicit none 
!
!***********************************************************************
!
  write(66,'(//,45x,a)') '######################################################################'
  write(66,'(45x,a70)')                             title
  write(66,'(45x,a)')    '######################################################################'
  write(66,'(/,50x,a,1pe10.4)')   'flow inlet       : flow = ',flowin
  write(66,'(50x,a,1pe10.4)')   'fluid density    :  den = ',densit
  write(66,'(50x,a,1pe10.4)')   'dynamic viscosity:  vis = ',viscos
  write(66,'(50x,a,0pf4.2)')    'alfa  parameter  : alfa = ',alfa
  write(66,'(50x,a,1pe10.4)')   'conv. criterion  :  sor = ',sormax
  write(66,'(50x,a,1pe10.4)')   '                  small = ',small
  write(66,'(50x,a,1pe10.4)')   '                  great = ',great

  if(lcal(iu))   write(66,'(/,50x,a,f5.2,a,f5.2,a,i3)') 'urf( u )=',urf(iu), '  gds( u )=',gds(iu), '  nsw=', nsw(iu)
  if(lcal(iv))   write(66,'(50x,a,f5.2,a,f5.2,a,i3)')   'urf( v )=',urf(iv), '  gds( v )=',gds(iu), '  nsw=', nsw(iv)
  if(lcal(iw))   write(66,'(50x,a,f5.2,a,f5.2,a,i3)')   'urf( w )=',urf(iw), '  gds( w )=',gds(iu), '  nsw=', nsw(iw)
  if(lcal(ip))   write(66,'(50x,a,f5.2,15x,a,i3)')      'urf( p )=',urf(ip),                       '   nsw=', nsw(ip)
  if(lcal(ite))  write(66,'(50x,a,f5.2,a,f5.2,a,i3)')   'urf( te)=',urf(ite),'  gds( te)=',gds(ite),'  nsw=', nsw(ite)
  if(lcal(ied))  write(66,'(50x,a,f5.2,a,f5.2,a,i3)')   'urf( ed)=',urf(ied),'  gds( ed)=',gds(ied),'  nsw=', nsw(ied)
  if(lcal(ivis)) write(66,'(50x,a,f5.2)')               'urf(vis)=',urf(ivis)
  if(lcal(ien))  write(66,'(50x,a,f5.2,a,f5.2,a,i3)')   'urf( t )=',urf(ien),'  gds( t )=',gds(ien),'  nsw=',nsw(ien)
  if(lcal(ivart)) then 
    write(66,'(50x,a,f5.2,a,f5.2,a,i3)') 'urf(var)=',urf(ivart),'  gds(var)=',gds(ivart),'  nsw=',nsw(ivart)
  end if
  if(lcal(icon)) then 
    write(66,'(50x,a,f5.2,a,f5.2,a,i3)') 'urf(con)=',urf(icon),'  gds(con)=',gds(icon),'  nsw=',nsw(icon)
  end if
  write(66,*)
  write(66,'(45x,a)') '================================================================='
  write(66,'(50x,a,e10.4)') 'time step= ',timestep
  write(66,'(45x,a)') '-----------------------------------------------------------------'
  if(.not.lbuoy) write(66,'(50x,a)') 'buoyancy activated:  - no -    '
  if(lbuoy) write(66,'(50x,a)') 'buoyancy activated:  - yes-    '
  if(boussinesq) then
    write(66,'(50x,a)') 'boussinesq aproximaton: - yes- '
  else !if(boussinesq.eq.0) 
    write(66,'(50x,a)') 'boussinesq aproximaton: - no - '
  endif
  write(66,'(45x,a)') '-----------------------------------------------------------------'
  write(66,'(50x,a,e10.4,2x,e10.4,2x,e10.4)') 'gravity: (x-y-z): ',gravx,gravy,gravz 
  write(66,'(45x,a)') '-----------------------------------------------------------------'
  write(66,'(50x,3(a,i3),a,i8)')'no. cells = ',ni-2,' x ',nj-2,' x ',nk-2,' = ',(ni-2)*(nj-2)*(nk-2)
  write(66,'(45x,a)') '================================================================='
  write(66,'(a3)') '   '
  write(66,'(a3)') '   '

  end subroutine
