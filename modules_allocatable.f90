module types
!
! define precision for floating-point real
!
   ! single precision
!   integer, parameter :: dp = selected_real_kind(6,37)
   ! double precision    
   integer, parameter :: dp = kind(1.0d0) !selected_real_kind(15,307) 
   ! quadruple precision
!   integer, parameter :: dp = selected_real_kind(33,4931)
   integer, parameter :: prec = dp
!
end module types

module parameters

  use types

  integer, parameter :: nphi=10  ! number of varibales to solve-fields
  integer, parameter :: ngit=1   ! for multigrid-how many grid levels
  integer, parameter :: iobst=0
  integer, parameter :: nobst=1
  integer, parameter :: nxo=1
  integer, parameter :: nyo=1
  integer, parameter :: nzo=1

  integer :: nicv, njcv, nkcv  ! important - no. of cvs - program reads this from the input file

  integer :: nxyo,nxzo,nyzo, &
             nx,ny,nz, &
             nxy,nxz,nyz,nxyz, &
             nxa,nya,nza, &
             nxya,nxza, nyza, &
             nxyza, &
             nxyzm

  ! Identifiers
  integer, parameter :: iu=1,iv=2,iw=3,ip=4,ite=5,ied=6,ien=7,ivis=8,ivart=9,icon=10

  real(prec), parameter :: small = 1e-30
  real(prec), parameter :: great = 1e+30

  ! MPI related variables
  integer :: myid
  integer :: nproc
  integer :: ierr
  integer :: this

end module parameters

module indexes
  !
  ! determine location on our structural grid
  !
  use types
  use parameters
   

    integer :: ni,nj,nk,nim,njm,nkm,nij,nik,njk,nijk, &
               ijs,iks,jks,ijks,icst,icen,&
               nig,njg,nkg,nimg,njmg,nkmg,nijg,nikg,njkg,nijkg,&
               ijsg,iksg,jksg,ijksg,icstg,iceng,&
               kgrid,iter,iterf,&
               imon,jmon,kmon,mpoints, &
               ijkmon,iim,jjm,kkm,ipr,jpr,kpr,ijkpr,idir
    integer :: nimm,njmm,nkmm

   integer, dimension(:), allocatable :: lig, li
   integer, dimension(:), allocatable :: lkg, lk

   ! parameters related to multigrid procedure:
   integer, dimension(ngit) :: nigit,njgit,nkgit,&
                               ijkgit,&
                               isbij,isbik,isbjk,&
                               lsg,lsr,lsi,mit

   ! those with nphi are related to each field that we calculate u,v,w,p,t,kin,dis...:
   integer, dimension(nphi) :: nsw
   logical, dimension(nphi) :: lcal
   real(prec), dimension(nphi) ::  sor,resor,snorin,prtinv,urf,urfr,urfm,gds

 
   !  stuff read from the input file
   real(prec) :: flowin,xmonin,c1,c2,c3,prm1,prt1,cmu,&
                 phit,sksi,eta,rcost,zero,dis,&
                 alfa,densit,viscos,sormax,slarge,&
                 bns1,bns2,&
                 c1asm,c2asm,c3asm,facnap,facflx

   real(prec), dimension(:), allocatable ::  bnuselt1,bnuselt2 ! dimension(nx*ny)

   integer :: iconvective_scheme
   logical :: lcds,lluds,lquds,lsmart,lavl,lmuscl,lumist,lgamma,lcds4
    
   integer :: iturbmodel ! Integer identifying turbulence model

   ! LOGICALS (yes/no) , MOSTLY READ FROM SIMULATION-INPUT FILE:
   logical :: lturb,lread,lwrite,ltest,louts,loute, &  ! turbulent simulation, read restart file, write restart file, print residual of the linear solver,.,..      
              lsor,lsol,ltransient, &             ! LTRANSIENT is TRUE for transient (non-stationary) simulations              
              levm,lasm,lles,ldes,lsgdh,lggdh,lafm, &  ! eddy-viscosity, algebraic stress model or LES, simple gradient or generalized gradient hypothesis, algerbaic flux model
              stdkeps,durbin,             &       ! which turbulence model: Standard k-epsilon, k-eps with Durbin limiter,
              rng, realizable, alphamodel, &      ! RNG-k-epsilon, Shih.et.al. Realizable k-epsilon, Kenjeres-Hanjalic ALPHA Seamless hybrid,
              LowRe_LB, &                         ! Low-Re Lam-Bremhorst k-epsilon,
              SST,Wilcox,  &                      ! Menter SST, Wilcox k-omega,
              earsm_wj,earsm_m,sas,lowre, &  ! EARSM by Wallin and Johansson, Menter 2009 verion of EARSM-WJ, Scale Adaptive Simulation-SAS hybrid, Low-Re version of any k-omega model.
              lstsq,gauss, &  ! gradient discretization approach
              bdf,cn,      &  ! control for the time-stepping algorithm
              simple,piso,pimple, &  ! control for the velocity-pressure coupling algorithm
              const_mflux     ! control for constant flow rate 
   ! PISO control parameter: no. of Piso corrections.
   integer :: ncorr
   ! PISO iteration no.: icorr=1..ncorr
   integer :: icorr
   ! No. of pressure-corrections; non-orthogonality correctors
   integer :: npcor
   ! Iteration no.: ipcorr=1..npcor
   integer :: ipcorr
   ! No. of iters. for iterative cell-centered gradient calculation
   integer :: nigrad

   ! Switches for periodicity; for channel its 'periodic_boundary' set to True.
   ! These enable last cell layer in one direction to exchange fluxes with first layer.
   logical :: periodic_boundary
   logical :: CoNumFix
   logical :: roughWall
   logical :: flux_limiter 

end module indexes

module geometry
!%%%%%%%%%%%%%%
   use types
    
   ! Geometry parameters defined cellwise
   real(prec), dimension(:), allocatable :: x, y, z
   real(prec), dimension(:), allocatable :: fx,fy,fz
   real(prec), dimension(:), allocatable :: xc,yc,zc,vol
   real(prec), dimension(:), allocatable :: ar1x,ar1y,ar1z
   real(prec), dimension(:), allocatable :: ar2x,ar2y,ar2z
   real(prec), dimension(:), allocatable :: ar3x,ar3y,ar3z
   real(prec), dimension(:), allocatable :: wallDistance
end module geometry

module coef
!%%%%%%%%%%%
   use types
     
  ! Coefficients resulting form fvm discretization:
  real(prec), dimension(:), allocatable :: ae,aw,an,as,at,ab,ap
  real(prec), dimension(:), allocatable :: sp,res,spv
  real(prec), dimension(:), allocatable :: su, sv, sw
  real(prec), dimension(:), allocatable :: apu, apv, apw
end module coef

module coefb
!%%%%%%%%%%%%%
   use types
   !    
   ! arrays used in sip solver
   !  
    real(prec), dimension(:), allocatable ::  bb,bs,bw,bp,bn,be,bt

end module coefb

module hcoef
!%%%%%%%%%%%%%
   use types
   !    
   ! arrays used in piso algorithm
   !  
      real(prec), dimension(:), allocatable ::  hb,hs,hw,hn,he,ht

end module hcoef

module omega_turb_models
!%%%%%%%%%%%%%%
! needed for menter sst model and wj-earsm-k-omega
!%%%%%%%%%%%%%%
   use types  
   real(prec) :: sigmk1,sigmk2,sigmom1,sigmom2, & ! constants prescribed in modinp routine
                 betai1,betai2,a1,alpha1,alpha2
   real(prec) :: alpha,betta,bettast
   real(prec), dimension(:), allocatable :: domega,alphasst,bettasst, &     ! cross diffusion coed and sst coefs
                                            qsas, &                         ! sas model additional production term
                                            prtinv_te,prtinv_ed, &          ! 1/sigma_k; 1/sigma_epsilon
                                            cmueff                          ! size(nxyza) ! cmueff the effective cmu for earsm
   real(prec), dimension(:,:), allocatable :: bij                           ! size(5,nxyza) reynolds stress anisotropy, used in earsm
   real(prec), dimension(:), allocatable :: lvk                             ! the on karman length scale for sas model.
end module omega_turb_models

module bc
!%%%%%%%%%%
   use types
  !
  !  for seting up boundary conditions 
  !

  !  mechanical field boundary conditions
  integer :: lbhlp
  integer, dimension(:), allocatable :: lbw,  lbe
  integer, dimension(:), allocatable :: lbs,  lbn
  integer, dimension(:), allocatable :: lbb,  lbt

  !  obstacle related
  logical :: low, loe, los, lon, lob, lot, lout

  ! energy equation related       
  integer :: lthermbc
  integer, dimension(:), allocatable :: lbwt,lbet
  integer, dimension(:), allocatable :: lbst,lbnt
  integer, dimension(:), allocatable :: lbbt,lbtt
  !  temperatures and heat fluxes
  integer, dimension(:,:), allocatable :: twest, teast,&
                                          qflxwest, qflxeast
  integer , dimension(:,:), allocatable :: tsouth, tnorth,&
                               qflxsouth, qflxnorth
  integer, dimension(:,:), allocatable :: tbottom, ttop,&
                               qflxbottom, qflxtop
  ! concentration related
  integer :: lconbc 
  integer, dimension(:), allocatable :: lbwc, lbec
  integer, dimension(:), allocatable :: lbsc, lbnc
  integer, dimension(:), allocatable :: lbbc, lbtc 
  integer, dimension(:,:), allocatable :: conwest,coneast
  integer , dimension(:,:), allocatable :: consouth,connorth
  integer, dimension(:,:), allocatable :: conbottom,contop
end module bc

module boundc
!%%%%%%%%%%
   use types
  !
  ! coefficient needed near the boundary 
  !

  integer ::      lw,lww
  real(prec) ::   suu,svu,swu,sup,svp,swp,cmu25,cmu75, &
                  cappa,elog,erough,zzero,ypl,tau,gent,sued,pranl,prant, &
                  pfun,hcoef,dx1,dx2,dy1,dy2,dz1,dz2,deln, &
                  t_hflx,psi,cu
end module boundc

module buoy
!%%%%%%%
   use types 
  !  
  ! Needed for buoyancy computations
  !
  logical :: lbuoy
  logical :: boussinesq
  real(prec) :: tref
  real(prec) :: beta
  real(prec) :: gravx, gravy, gravz
  real(prec) :: tybw,tybe,tybs,tybn,tybb,tybt
  real(prec) :: temperaturatop
  real(prec) :: facvis
end module buoy

module variables
!%%%%%%%%%%%%%%
   use types
    
   ! these are cellwise defined variables, that is - the fields
   real(prec), dimension(:), allocatable :: u ,v ,w
   real(prec), dimension(:), allocatable :: f1 ,f2 ,f3
   real(prec), dimension(:), allocatable :: p ,pp ,te ,ed ,t
   real(prec), dimension(:), allocatable :: vis ,den ,gen ,vart
   real(prec), dimension(:), allocatable :: edd ,utt ,vtt ,wtt
   real(prec), dimension(:), allocatable :: ret ,con
   real(prec), dimension(:), allocatable :: yplus ,upl,kplus 
   real(prec), dimension(:), allocatable :: uu ,uv ,uw, &
                                                vv ,vw, &
                                                    ww

   ! friction velocity
   real(prec) :: utau 
   ! Magnitude of the bulk velocity,and pressure grad that will drive the constant mass-flux flow (cmf)
   real(prec) :: magUbar, gradPcmf
   ! Continuity errrs
   real(prec) :: LocalContErr,sumLocalContErr, globalContErr, cumulativeContErr
   ! Fixed value for Courant number - set in modinp for now - may be read in input
   real(prec) :: CoNumFixValue

   !  Related to seamless-alpha turbulence model:
   real(prec), dimension(:), allocatable :: alph, al_les, al_rans, diff          !  Related to seamless-alpha turbulence model
   real(prec), dimension(:), allocatable :: timelimit                            ! Durbin Time-scale limiter
   real(prec), dimension(:), allocatable :: strain                               ! Strain magnitude
   real(prec), dimension(:), allocatable :: Vorticity                            ! Vorticity magnitude

end module variables

module nusselt
!%%%%%%%%%%%%
  use types
  use parameters
   
  real(prec), dimension(:), allocatable :: bnusmeanw,bnusmeane, &
                                           bnusmeanb, bnusmeant

  real(prec), dimension(:), allocatable :: bnuslocalw,bnuslocale, &
                                           bnuslocalb,bnuslocalt  
end module nusselt

module title_mod
!%%%%%%%%%%%%
  use types
  use parameters
   
  character(len=70) :: title
  character(len=4), dimension(10) ::  chvar = (/'  U ', '  V ', '  W ', '  P ', ' TE ', ' ED ', ' IEN', ' VIS', 'VART', ' CON' /)
  character(len=7), dimension(10) ::  chvarSolver = (/'U      ', 'V      ', 'W      ', 'p      ', 'k      ', 'epsilon','Energy ',&
                                                  &   'Visc   ', 'VarTemp', 'Conc   ' /)
  character(len=100):: input_file,inlet_file,grid_file,monitor_file,restart_file,out_folder_path
end module title_mod

module time_mod
! Variables for transient simulation.
   use types
   use parameters
   !
   ! timesteping control
   !
   integer:: numstep,itime, &
             nzapis,maxit
   real(prec) :: timestep,time,btime  
   real(prec) :: CoNum,meanCoNum ! courant number.         

   ! values from n-1 timestep
   real(prec), dimension(:), allocatable :: uo,vo,wo, &
                                            to,teo,edo, &
                                            varto,cono
   ! values from n-2 time step
   real(prec), dimension(:), allocatable :: uoo,voo,woo, &
                                            too,teoo,edoo, &
                                            vartoo,conoo
end module time_mod

module statistics
! Variables for collecting statistics
! these are ensamble averaged values over time steps
! used in t_rans-urans-hybrid rans/les approach.
   use types

   integer :: n_sample,istat,ifreq

   real(prec), dimension(:), allocatable :: u_aver,v_aver,w_aver, &
                                            te_aver,t_aver, &
                                            uu_aver,vv_aver,ww_aver, &
                                            uv_aver,uw_aver,vw_aver, &
                                            ut_aver,vt_aver,wt_aver, &
                                            tt_aver
end module statistics

module wall
!%%%%%%%%%%%%%%%%%
   use types

   real(prec), dimension(:), allocatable :: gentw,gente, &
                                            suedw,suede
   real(prec), dimension(:), allocatable :: gents,gentn, &
                                            sueds,suedn
   real(prec), dimension(:), allocatable :: gentb,gentt, &
                                            suedb,suedt 
end module wall


module printing
!%%%%%%%%%%%
  logical :: lpri,lprj,lprk
end module printing

module inlet
!%%%%%%%%%%%%
   use types
   real(prec) :: uin,vin,win,tein,edin,tin,vartin,conin   
   real(prec), dimension(:), allocatable :: u_inl,v_inl,w_inl,p_inl, &
                                            t_inl,te_inl,ed_inl

end module inlet

module gradients
   use types

   real(prec),dimension(:,:), allocatable :: gradu,gradv,gradw, &  ! size(nxyza,3)
                                             gradte,graded, &
                                             gradp,gradt, &
                                             gradvart,gradcon
  real(prec),dimension(:,:), allocatable ::  d   !  d(6,nxyza) - when using bn, or dm version of the subroutine
  real(prec),dimension(:,:,:), allocatable ::  dqr             ! when using qr version of the subroutine size(3,6,nxyza)!
end module gradients

module obstacle
  !
  ! variables for defining obsacles.
  !
  use types
   

  integer, dimension(:), allocatable :: ios,ioe,jos,joe,kos,koe

  real(prec),dimension(:,:), allocatable :: dnow, dnoe,  &
                                            genow, genoe,  &
                                            suedow, suedoe

  real(prec),dimension(:,:), allocatable :: dnos, dnon,  &
                                            genos, genon,  &
                                            suedos, suedon

  real(prec),dimension(:,:), allocatable :: dnob ,dnot , &
                                            genot ,genob, &
                                            suedot ,suedob 

  !  used in whole program to store molecular dynamic viscosity 
  real(prec),dimension(:), allocatable :: visob

  real(prec),dimension(:,:), allocatable :: tobwest, &
                                            tobeast, &
                                            cobwest, &
                                            cobeast

  real(prec),dimension(:,:), allocatable :: tobsouth, &
                                            tobnorth, &
                                            cobsouth, &
                                            cobnorth

  real(prec),dimension(:,:), allocatable :: tobbottom, &
                                                   tobtop, &
                                                   cobbottom, &
                                                   cobtop

  integer, dimension(:,:), allocatable :: typobw, typobe, &
                                          typobtw, typobte, &
                                          typobcw, typobce

  integer, dimension(:,:), allocatable :: typobs, typobn, &
                                          typobts, typobtn, &
                                          typobcs, typobcn                    

  integer, dimension(:,:), allocatable :: typobb, typobt, &
                                          typobtb, typobtt, & 
                                          typobcb, typobct

end module obstacle