subroutine read_input
!********************************************************************************
!  Open & read input file
!********************************************************************************
  use types
  use parameters
  use indexes
  use variables, only: conumfixvalue, gradPcmf
  use boundc
  use title_mod
  use buoy
  use time_mod
  use printing
  use inlet

  implicit none

  include 'mpif.h'

  integer :: i


! Root processor opens files
if (myid .eq. 0) then

  open(unit=5,file=input_file)
  rewind 5

  read(5,'(a70)') title 
  read(5,*) lread,lwrite,ltest
  read(5,*) louts,loute,lpri,lprj,lprk
  ! selection of dependent variables
  read(5,*) (lcal(i),i=1,nphi)
  ! Iteration and conv. limits, output control, fluid properties
  read(5,*) imon,jmon,kmon,ipr,jpr,kpr,mpoints
  read(5,*) slarge,sormax
  read(5,*) densit,viscos,dis
  read(5,*) pranl,tref,beta
  read(5,*) lbuoy,gravx,gravy,gravz,boussinesq
  ! Turbulence
  read(5,*) roughwall,erough,zzero
  read(5,*) phit,sksi,eta,rcost,facnap,facflx
  read(5,*) ltransient,btime
  read(5,*) levm,lasm,lles,ldes
  read(5,*) lsgdh,lggdh,lafm
  read(5,*) iturbmodel
  ! Initialisation of variables
  read(5,*) uin,vin,win,tein,edin,tin,vartin,conin
  read(5,*) iconvective_scheme
  read(5,*) (gds(i),i=1,nphi)
  read(5,*) (urf(i),i=1,nphi)
  read(5,*) (sor(i),i=1,nphi)
  read(5,*) (nsw(i),i=1,nphi)
  read(5,*) numstep,timestep,nzapis,maxit
  read(5,*) gauss, lstsq
  read(5,*) npcor, nigrad
  read(5,*) bdf,cn
  read(5,*) simple,piso,pimple,ncorr
  read(5,*) periodic_boundary,const_mflux,gradpcmf
  read(5,*) conumfix, conumfixvalue

  close (5)

  ! Create an input file reading log:
  write(66,'(a)') ' input file log: '
  write(66,'(a)') '--------------------------------------------------------------------------------'
  write(66,'(a70)') title
  write(66,'(3(l1,1x),5x,a)') lread,lwrite,ltest,'read3,writ3,ltest'
  write(66,'(5(l1,1x),5x,a)') louts,loute,lpri,lprj,lprk,'louts,loute,lpri,lprj,lprk'
  write(66,'(10(l1,1x),5x,a)') (lcal(i),i=1,nphi),'(ical(i),i=1,nphi),ien=7!!!,ivis=8,ivart=9;icon=10!'
  write(66,'(7(i3,1x),5x,a)') imon,jmon,kmon,ipr,jpr,kpr,mpoints,'imon,jmon,kmon,ipr,jpr,kpr,mpoints'
  write(66,'(2(es11.4,1x),5x,a)') slarge,sormax,'slarge,sormax'
  write(66,'(3(es11.4,1x),a)') densit,viscos,dis,'densit,viscos,dis'
  write(66,'(3(es11.4,1x),a)') pranl,tref,beta,'pranl,tref,beta'
  write(66,'(l1,1x,3f5.2,1x,l1,1x,a)') lbuoy,gravx,gravy,gravz,boussinesq,'lbuoy,gravx,gravy,gravz,boussinesq'
  write(66,'(l1,1x,f5.2,1x,es11.4,1x,a)') roughwall,erough,zzero,'roughwall,erough,zzero'
  write(66,'(6(f4.2,1x),a)') phit,sksi,eta,rcost,facnap,facflx,'phit,sksi,eta,rcost,facnap,facflx'
  write(66,'(l1,1x,f4.2,1x,a)') ltransient,btime,'ltransient,btime'
  write(66,'(4(l1,1x),a)') levm,lasm,lles,ldes,'levm,lasm,lles,ldes'
  write(66,'(3(l1,1x),a)') lsgdh,lggdh,lafm,'lsgdh,lggdh,lafm'
  write(66,'((i2,1x),a)') iturbmodel,'iturbmodel'
  write(66,'(8(es11.4,1x),a)') uin,vin,win,tein,edin,tin,vartin,conin,'uin,vin,win,tein,edin,tin,vartin,conin'
  write(66,'((i2,1x),a)') iconvective_scheme,'iconvective_scheme'
  write(66,'(10(f4.2,1x),a)') (gds(i),i=1,nphi),'(gds(i),i=1,nphi), muscl velocity, cds other'
  write(66,'(10(f4.2,1x),a)') (urf(i),i=1,nphi),'(urf(i),i=1,nphi)'
  write(66,'(10(es9.2,1x),a)') (sor(i),i=1,nphi),'(sor(i),i=1,nphi)'
  write(66,'(10(i3,1x),a)') (nsw(i),i=1,nphi),'(nsw(i),i=1,nphi)'
  write(66,'(i6,1x,es9.2,1x,i5,1x,i4,1x,a)') numstep,timestep,nzapis,maxit,'numstep,timestep,nzapis,maxit'
  write(66,'(l1,1x,l1,1x,a)') gauss, lstsq,'gauss, lstsq'
  write(66,'(i1,1x,i1,1x,a)') npcor, nigrad,'npcor, nigrad'
  write(66,'(2(l1,1x),1x,a)') bdf,cn,'bdf,cn'
  write(66,'(3(l1,1x),i1,1x,a)') simple,piso,pimple,ncorr,'simple,piso,pimple,ncorr'
  write(66,'(2(l1,1x),es11.4,5x,a)') periodic_boundary,const_mflux,gradpcmf,'periodic_boundary, const_mflux, gradpcmf'
  write(66,'(l1,es11.4,5x,a)') conumfix, conumfixvalue,'conumfix, conumfixvalue'
  write(66,'(a)') '--------------------------------------------------------------------------------'
  write(66,'(a)') ' '

endif

! Broadcast input data to other processes
  call MPI_BCAST(title,70,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lread,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lwrite,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(ltest,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(louts,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(loute,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lpri,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lprj,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lprk,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lcal,NPHI,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(IMON,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(JMON,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(KMON,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(IPR,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(JPR,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(KPR,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(MPOINTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  ! Treba naci kom procesu pripadaju monitoring tacke - pogledaj getpidlm.f kod Sase.

  call MPI_BCAST(slarge,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(sormax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(densit,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(viscos,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(dis,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(pranl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(tref,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(beta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lbuoy,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gravx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gravy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gravz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(boussinesq,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(roughWall,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(erough,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(zzero,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(phit,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(sksi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(eta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(rcost,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(facnap,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(facflx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(ltransient,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(btime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(levm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lasm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lles,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(ldes,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lsgdh,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lggdh,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lafm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(iturbmodel,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  ! Set turbulence model:
  stdkeps = .false.
  durbin = .false.
  rng = .false.
  realizable = .false.
  alphamodel = .false.
  LowRe_LB = .false.
  wilcox = .false.
  sst = .false.
  sas = .false.
  earsm_wj = .false.
  earsm_m = .false.
  lowre = .false.

  if ( iturbmodel == 1 ) then
    stdkeps = .true.
  elseif ( iturbmodel == 2 ) then
    durbin = .true.
  elseif ( iturbmodel == 3 ) then
    rng = .true.
  elseif ( iturbmodel == 4 ) then
    realizable = .true.
  elseif ( iturbmodel == 5 ) then
    alphamodel = .true.
  elseif ( iturbmodel == 6 ) then
    LowRe_LB = .true.
  elseif ( iturbmodel == 7 ) then
    wilcox = .true.
  elseif ( iturbmodel == 8 ) then
    sst = .true.
  elseif ( iturbmodel == 9 ) then
    sas = .true.
  elseif ( iturbmodel == 10 ) then
    earsm_wj = .true.
  elseif ( iturbmodel == 11 ) then
    earsm_m = .true.
  elseif ( iturbmodel == 12 ) then
    sst = .true.
    lowre = .true.
  ! elseif ( iturbmodel == 13 ) then
  !    = .true.
  ! elseif ( iturbmodel == 14 ) then
  !    = .true.
  ! elseif ( iturbmodel == 15 ) then
  !    = .true.
  ! elseif ( iturbmodel == 16 ) then
  !    = .true.
  endif


  call MPI_BCAST(uin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(vin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(win,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(tein,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(edin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(tin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(vartin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(conin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)


  call MPI_BCAST(iconvective_scheme,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  ! Convective scheme choice
  lcds = .false.
  lluds = .false.
  lquds = .false.
  lsmart = .false.
  lavl = .false.
  lmuscl = .false.
  lumist = .false.
  lgamma = .false.
  lcds4 = .false.

  if ( iconvective_scheme == 1 ) then
    lcds = .true.
  elseif ( iconvective_scheme == 2 ) then
    lluds = .true.
  elseif ( iconvective_scheme == 3 ) then
    lquds = .true.
  elseif ( iconvective_scheme == 4 ) then
    lsmart = .true.
  elseif ( iconvective_scheme == 5 ) then
    lavl = .true.
  elseif ( iconvective_scheme == 6 ) then
    lmuscl= .true.
  elseif ( iconvective_scheme == 7 ) then
    lumist = .true.
  elseif ( iconvective_scheme == 8 ) then
    lgamma = .true.
  elseif ( iconvective_scheme == 9 ) then
    lcds4 = .true.
  ! elseif ( iconvective_scheme == 10 ) then
  !    = .true.
  ! elseif ( iconvective_scheme == 11 ) then
  !    = .true.
  ! elseif ( iconvective_scheme == 12 ) then
  !    = .true.
  endif


  call MPI_BCAST(GDS(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(URF(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(SOR(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(NSW(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(numstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(timestep,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(nzapis,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(maxit,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)


  call MPI_BCAST(gauss,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lstsq,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(npcor,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(nigrad,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(bdf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(cn,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(simple,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(piso,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(pimple,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(ncorr,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(periodic_boundary,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(const_mflux,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gradpcmf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(conumfix,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(conumfixvalue,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

if (myid == 0) then
  write(66,*)' '
  write(66,*)'  ->input data O.K.'
  write(66,*)' '
endif

end subroutine read_input