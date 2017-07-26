!***********************************************************************
!
  SUBROUTINE INIT
!
!***********************************************************************
! Precomute and store cell face areas
! Set Coefficient values for Turbulence models
! Set Initial timestepping control values
! Various initialisations
!     Field Initialisation
! Read Restart File And Set Field Values
! Initial Gradient Calculation
! Calculate distance to the nearest wall.
!
!***********************************************************************
  USE TYPES
  USE PARAMETERS
  USE INDEXES
  USE GEOMETRY
  USE COEF
  USE VARIABLES
  USE BC
  USE BOUNDC
  USE TITLE_MOD
  USE BUOY
  USE TIME_MOD
  USE PRINTING
  USE OBSTACLE
  USE INLET
  USE GRADIENTS
  USE OMEGA_Turb_Models
  use fieldManipulation

  IMPLICIT NONE
!
!***********************************************************************
!
  integer :: I, J, K, INP,intc,inbc
  integer :: inn,inc,ins,inb,inbs
  real(prec) :: xf,yf,zf
  real(prec) :: ednom, tenom
  real(prec) :: dxet,dxzd,dyet,dyzd,dzet,dzzd
  real(prec) :: perturb
  ! real(prec) :: zz
  ! real(prec) :: UST 

  ! Power law profile for hills:
  !REAL(PREC), PARAMETER :: blthick = 0.2 
  !REAL(PREC), PARAMETER :: expp = 0.135  
  !REAL(PREC), PARAMETER :: Umag = 5.5 

  real(prec), parameter :: pi = 4.*atan(1.) 


!
! Precomute and store cell face areas
!

  ! East faces  
  do k=2,nkm
  do i=1,nim
  do j=2,njm

    inp=lk(k)+li(i)+j
    ins=inp-1
    inb=inp-nij
    inbs=inb-1

    dxet=0.5*(x(inp)-x(ins)+x(inb)-x(inbs))
    dyet=0.5*(y(inp)-y(ins)+y(inb)-y(inbs))
    dzet=0.5*(z(inp)-z(ins)+z(inb)-z(inbs))

    dxzd=0.5*(x(inp)-x(inb)+x(ins)-x(inbs))
    dyzd=0.5*(y(inp)-y(inb)+y(ins)-y(inbs))
    dzzd=0.5*(z(inp)-z(inb)+z(ins)-z(inbs))

    ar1x(inp)=dyet*dzzd-dyzd*dzet
    ar1y(inp)=dxzd*dzet-dxet*dzzd
    ar1z(inp)=dxet*dyzd-dyet*dxzd


    ! interpolation factor lambda = pj'/ pn 
    inn=inp+nj

    xf = 0.25d0*(x(inp)+x(ins)+x(inb)+x(inbs))
    yf = 0.25d0*(y(inp)+y(ins)+y(inb)+y(inbs))
    zf = 0.25d0*(z(inp)+z(ins)+z(inb)+z(inbs))

    !Treba naci j' tacku na kojoj se sece ravan face-a i linija spajanja celijskih centara => find_intersection_point
    call find_intersection_point( &
  !            plane defined by face corner, bottom and south:
           x(inp),y(inp),z(inp),&
           x(ins),y(ins),z(ins), &
           x(inb),y(inb),z(inb), &
  !            line defined by cell center and neighbour center:
           xc(inp),yc(inp),zc(inp), &
           xc(inn),yc(inn),zc(inn), &
  !            intersection point:
           xf,yf,zf &
           )

    fx(inp) =   sqrt( (xf-xc(inp))**2 + (yf-yc(inp))**2 + (zf-zc(inp))**2 )  &
        / ( sqrt( (xc(inn)-xc(inp))**2 + (yc(inn)-yc(inp))**2 + (zc(inn)-zc(inp))**2 ) + small )


  enddo
  enddo
  enddo

  ! North faces
  do k=2,nkm
  do i=2,nim
  do j=1,njm

    inp=lk(k)+li(i)+j-nj
    ins=lk(k)+li(i)+j
    inb=inp-nij
    inbs=ins-nij

    dxet=0.5*(x(inp)-x(ins)+x(inb)-x(inbs))
    dyet=0.5*(y(inp)-y(ins)+y(inb)-y(inbs))
    dzet=0.5*(z(inp)-z(ins)+z(inb)-z(inbs))

    dxzd=0.5*(x(inp)-x(inb)+x(ins)-x(inbs))
    dyzd=0.5*(y(inp)-y(inb)+y(ins)-y(inbs))
    dzzd=0.5*(z(inp)-z(inb)+z(ins)-z(inbs))

    ar2x(inp)=dyet*dzzd-dyzd*dzet
    ar2y(inp)=dxzd*dzet-dxet*dzzd
    ar2z(inp)=dxet*dyzd-dyet*dxzd


    ! interpolation factor lambda = pj'/ pn 
    inc=lk(k)+li(i)+j
    inn=inp+1

    xf = 0.25d0*(x(inp)+x(ins)+x(inb)+x(inbs))
    yf = 0.25d0*(y(inp)+y(ins)+y(inb)+y(inbs))
    zf = 0.25d0*(z(inp)+z(ins)+z(inb)+z(inbs))

    !Treba naci j' tacku na kojoj se sece ravan face-a i linija spajanja celijskih centara => find_intersection_point
    call find_intersection_point( &
  !            plane defined by face corner, bottom and south:
           x(inp),y(inp),z(inp),&
           x(ins),y(ins),z(ins), &
           x(inb),y(inb),z(inb), &
  !            line defined by cell center and neighbour center:
           xc(inc),yc(inc),zc(inc), &
           xc(inn),yc(inn),zc(inn), &
  !            intersection point:
           xf,yf,zf &
           )

    fy(inc) =   sqrt( (xf-xc(inc))**2 + (yf-yc(inc))**2 + (zf-zc(inc))**2 )  &
        / ( sqrt( (xc(inn)-xc(inc))**2 + (yc(inn)-yc(inc))**2 + (zc(inn)-zc(inc))**2 ) + small )


  enddo
  enddo
  enddo

  ! Top faces
  do k=1,nkm
  do i=2,nim
  do j=2,njm

    inp=lk(k)+li(i)+j-nj
    ins=inp-1
    inb=lk(k)+li(i)+j
    inbs=inb-1

    dxet=0.5*(x(inp)-x(ins)+x(inb)-x(inbs))
    dyet=0.5*(y(inp)-y(ins)+y(inb)-y(inbs))
    dzet=0.5*(z(inp)-z(ins)+z(inb)-z(inbs))

    dxzd=0.5*(x(inp)-x(inb)+x(ins)-x(inbs))
    dyzd=0.5*(y(inp)-y(inb)+y(ins)-y(inbs))
    dzzd=0.5*(z(inp)-z(inb)+z(ins)-z(inbs))

    ar3x(inp)=dyet*dzzd-dyzd*dzet
    ar3y(inp)=dxzd*dzet-dxet*dzzd
    ar3z(inp)=dxet*dyzd-dyet*dxzd


  !   ! interpolation factor lambda = pj'/ pn 
  !   inc=lk(k)+li(i)+j
  !   inn=inp+nij

  !   xf = 0.25d0*(x(inp)+x(ins)+x(inb)+x(inbs))
  !   yf = 0.25d0*(y(inp)+y(ins)+y(inb)+y(inbs))
  !   zf = 0.25d0*(z(inp)+z(ins)+z(inb)+z(inbs))

  !   !Treba naci j' tacku na kojoj se sece ravan face-a i linija spajanja celijskih centara => find_intersection_point
  !   call find_intersection_point( &
  ! !            plane defined by face corner, bottom and south:
  !          x(inp),y(inp),z(inp),&
  !          x(ins),y(ins),z(ins), &
  !          x(inb),y(inb),z(inb), &
  ! !            line defined by cell center and neighbour center:
  !          xc(inc),yc(inc),zc(inc), &
  !          xc(inn),yc(inn),zc(inn), &
  ! !            intersection point:
  !          xf,yf,zf &
  !          )
    
  !   fz(inp) =   sqrt( (xf-xc(inc))**2 + (yf-yc(inc))**2 + (zf-zc(inc))**2 )  &
  !       / ( sqrt( (xc(inn)-xc(inc))**2 + (yc(inn)-yc(inc))**2 + (zc(inn)-zc(inc))**2 ) + small )
      

  enddo
  enddo
  enddo
!......................................................................./

!
! Calculation of interpolation factors - new method
!.......................................................................

  ! call intfac(nj,1,nij, fx)
  ! call intfac(1,nj,nij, fy)
  ! call intfac(nij,1,nj, fz)




!
! Set Coefficient values for Turbulence models
!

  ! LTURB set to True if TKE field is calculated -> true for all RANS models
  LTURB=LCAL(ITE).OR.LCAL(IED)

  ! Flux limiter?
  if (lsmart.or.lavl.or.lmuscl.or.lumist.or.lgamma) then
    flux_limiter = .true.
  else
    flux_limiter = .false.
  endif

!=====Factor for SIP solver==================================
  ALFA = 0.92_dp

!=====Law of the wall parameters=============================
  CAPPA = 0.41_dp 
  ELOG = 8.432_dp

!=====Coefficients for Sasa's algebraic flux model
  C1asm = 1.8_dp
  C2asm = 0.6_dp 
  C3asm = 0.6_dp

!=====STANDARD k-epsilon=============================
  IF (STDKEPS) THEN
  CMU = 0.09_dp   
  C1 = 1.44_dp
  C2 = 1.92_dp
  C3 = 1.44_dp
!=====STANDARD k-epsilon=============================
  ENDIF

!=====RENORMALIZATION GROUP (RNG) k-epsilon=============================
  IF (RNG) THEN
  CMU = 0.0845   ! za RNG model!!!
  C1 = 1.42  ! za RNG model!!!
  C2 = 1.68  ! za RNG model!!!
!=====RENORMALIZATION GROUP (RNG) k-epsilon=============================
  ENDIF

!=====REALIZABLE k-epsilon==============================================
  IF (REALIZABLE) THEN
  C1 = 1.44  ! za Realizable k-eps model !!!
  C2 = 1.9   ! za Realizable k-eps model !!!
!=====END:REALIZABLE k-epsilon==========================================
  END IF

!====================================================
!    Define turbulence model constants here.
!    k-omega model: Sigma_k=2.0; Sigma_omega=2.0
!    REFERENCES:
!    Wilcox1998:
!    * Wilcox, D. C., "Reassessment of the Scale-Determining Equation for Advanced Turbulence Models," AIAA Journal, Vol. 26, No. 11, 1988, pp. 1299-1310.
!    * Wilcox, D. C., Turbulence Modeling for CFD, 1st edition, DCW Industries, Inc., La Canada CA, 1993. 
!    * Wilcox, D. C., Turbulence Modeling for CFD, 2nd edition, DCW Industries, Inc., La Canada CA, 1998. <- This is where Wilcox1998 formulation comes from.
!    Wilcox2006:
!    * Wilcox, D. C., "Formulation of the k-omega Turbulence Model Revisited," AIAA Journal, Vol. 46, No. 11, 2008, pp. 2823-2838.
!    * Wilcox, D. C., Turbulence Modeling for CFD, 3rd edition, DCW Industries, Inc., La Canada CA, 2006. 
!====================================================
  IF (Wilcox) THEN
!.....Wilcox1998 (These are in Fluent 13):
  ALPHA=13./25.
  BETTA=0.072
  BETTAST=0.09
!.....Wilcox2006:
!%  ALPHA=13./25.
!%  BETTA=0.0708
!%  BETTAST=0.09
!=====================================================
  ENDIF
!====================================================
!     Define SST, ad SAS-SST model constants.
!     REFERENCES:
!     * ANSYS FLUENT THheory Guide p.71
!     * Menter, F. R., "Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications," AIAA Journal, Vol. 32, No. 8, August 1994, pp. 1598-1605. 
!     * Menter, F. R., Kuntz, M., and Langtry, R., "Ten Years of Industrial Experience with the SST Turbulence Model," Turbulence, Heat and Mass Transfer 4, ed: K. Hanjalic, Y. Nagano, and M. Tummers, Begell House, Inc., 2003, pp. 625 - 632. 
!=====================================================
  IF (SST.OR.SAS) THEN
  SIGMK1=1./1.176D0
  SIGMK2=1.0D0
  SIGMOM1=1./2.0D0
  SIGMOM2=1./1.168D0
  BETAI1=0.075
  BETAI2=0.0828
  A1=0.31D0
!.....SST-1994 coefficients
!%  ALPHA1=(BETAI1/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMAOM1)
!%  ALPHA2=(BETAI2/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMAOM2)
!.....SST-2003 coefficients. The first is higher than the original constant
!     definition by approximately 0.43%, and the second is lower by less than 0.08%. 
  ALPHA1=5./9.
  ALPHA2=0.44
  BETTAST=0.09
  END IF

  IF (EARSM_WJ) THEN
  SIGMK1=1.1d0
  SIGMK2=1.1d0
  SIGMOM1=0.53d0
  SIGMOM2=1.0
  BETAI1=0.0747
  BETAI2=0.0828
  A1=0.31D0
  ALPHA1=0.518d0
  ALPHA2=0.44
  BETTAST=0.09
  END IF

  IF (EARSM_M) THEN
  SIGMK1=0.5d0
  SIGMK2=1.1d0
  SIGMOM1=0.5d0
  SIGMOM2=0.856d0
  BETAI1=0.075
  BETAI2=0.0828
  A1=0.31D0
  BETTAST=0.09
  ALPHA1=(BETAI1/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMOM1)
  ALPHA2=(BETAI2/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMOM2)
  END IF

  CMU25=DSQRT(DSQRT(CMU))
  CMU75=CMU25**3

!==========================
!.....PRANDTL NUMBER FOR FLUID
!.....PRANL FOR WATER, PRANL=7.0
!.....PRANL FOR AIR  , PRANL=0.7
!==========================
  PRM1=1./PRANL
  PRANT=0.86
  PRT1=1./PRANT
!
!.....RECIPROCAL VALUES OF PRANDTL NUMBERS
  DO I=1,NPHI
  PRTINV(I)=1.0d0
  END DO

!=====TURBULENT MODEL SIGMA's===========================================
  IF(STDKEPS.OR.DURBIN) THEN
!.....[Standard k-epsilon Coefficient: ]
  PRTINV(IED)=1./1.3d0
!     [Beljaars (Askervein hill paper), 1987 Coefficients: ]
  !PRTINV(IED)=1./1.85
  ENDIF

!=====RENORMALIZATION GROUP (RNG) k-epsilon=============================
  IF (RNG) THEN
  PRTINV(ITE)=1./0.7194d0
  PRTINV(IED)=1./0.7194d0
!=====END:RENORMALIZATION GROUP (RNG) k-epsilon=========================
  ENDIF

!=====REALIZABLE k-epsilon==============================================
  IF (REALIZABLE) THEN
  PRTINV(ITE) = 1.0d0
  PRTINV(IED) = 1./1.20d0
!=====END:REALIZABLE k-epsilon==========================================
  ENDIF

!=====Wilcox k-omega==============================================
!.....Wilcox1998:
  IF (Wilcox) THEN
  PRTINV(IED)=0.5
  PRTINV(ITE)=0.5
!.....Wilcox2006:
!%  PRTINV(IED)=0.5
!%     PRTINV(ITE)=0.6
!=====END:Wilcox k-omega==========================================
  ENDIF

!.....Prandtl-Schmidt number for energy equation:
  PRTINV(IEN)=1./PRANL
!  IF(LTURB) PRTINV(IEN)=1./PRANT
!@  PRRAT=PRANL/PRANT
!@  PFUN=9.24*(PRRAT**0.75-1.0)*(1.0+0.28*EXP(-0.007*PRRAT))
!@  WRITE(66,*)'pfun: ',pfun

!.....RECIPROCAL VALUES OF UNDERRRELAXATION FACTORS
  DO I=1,NPHI
  URFR(I)=1./URF(I)
  URFM(I)=1.-URF(I)
  ENDDO
  TENOM=UIN**2+VIN**2+WIN**2
  IF(TENOM.LT.SMALL) TENOM=1.
  EDNOM=TENOM*DSQRT(TENOM)/(0.1*DIS) ! FIXME add for omega models also



! 6)  Set Initial timestepping control values
  ZERO=0.0D0
  ITIME=0
  TIME=-TIMESTEP

  !Set to zero cumulative error in continuity
  cumulativeContErr = 0.0_dp



! 7)  Various initialisations

!     7.0)  Field parameter Initialisation

! # Bulk velocity - estimated from Reynolds number, and written in 'input' file
  magUbar = UIN !!!<----veoma bitno za const_mflux flow!

!     7.1)  Field Initialisation - Turbulence kinetic energy and dissipation at inlet

  inbc = li(1)+1
  intc = lk(nk)+li(1)+1
  IF(Wilcox.or.SST.or.SAS.or.EARSM_WJ.or.EARSM_M) THEN
    TEIN=1e-6*UIN*UIN
!.......0.056894 = 0.08948-0.028 -> visina kanala!
    EDIN=5.*UIN/(zc(intc)-zc(inbc))
  ELSEIF(stdkeps.or.Durbin.or.RNG.or.Realizable.or.LowRe_LB) THEN 
    TEIN=1.5*(0.05*UIN)**2 ! for Ti=0.05  
    EDIN=cmu75*tein**1.5/(zc(intc)-zc(inbc)) 
  ENDIF

!     7.2)  Field Initialisation - reading inlet file

!.....INITIALIZATION OF FIELD VARIABLES FROM INLET FILE:
!.....READ INLET VALUES FIRST:
  ! DO K=2,NKM
  ! DO J=2,NJM-1
  !   READ(7,*)
  ! ENDDO
  ! READ(7,*)       U_INL(K),V_INL(K),W_INL(K),P_INL(K),TE_INL(K),ED_INL(K),T_INL(K)
  ! WRITE(66,'(7(es11.4))') U_INL(K),V_INL(K),W_INL(K),P_INL(K),TE_INL(K),ED_INL(K),T_INL(K)
  ! ENDDO
  ! REWIND 7

!     7.3)  Field Initialisation

!-----Field initialisation loop over inner & ghostcells--------------------------------
  DO K=2,NKM
  DO I=2,NIM
  DO J=2,NJM

  INP=LK(K)+LI(I)+J
  INBC=LK(2)+LI(I)+J

!.....INITIALIZATION OF FIELD VARIABLES FROM INLET FILE:
  ! U(INP)=U_INL(K)
  ! V(INP)=V_INL(K)
  ! W(INP)=W_INL(K)
  ! TE(INP)=TE_INL(K)
  ! ED(INP)=ED_INL(K)

  !ed(inp)=ed(inp)/(cmu*te(inp)+small) ! k-eps --> k-omega

!.......Channel flow:
  ! Random number based fluctuation of mean profile    
  CALL init_random_seed()
  CALL RANDOM_NUMBER(perturb)
  
  perturb = 0.9+perturb/5. ! Max perturbation is 10% of mean profile

  U(INP) = perturb*magUbar*(1.25*(1.-zc(inp)**4))
  V(INP) = small
  W(INP) = small

!.....INITIALIZATION OF FIELD VARIABLES FROM INPUT FILE:
  !U(INP)=UIN
  !V(INP)=VIN
  !W(INP)=WIN
  !TE(INP)=TEIN
  !ED(INP)=EDIN!/(CMU*TE(INP))

!.....USER INITIALIZATION:

!.....Power law profile:
  !U(INP) = Umag*(zc(inp)/blthick)**expp
  !V(INP) = small
  !W(INP) = small
  !TE(INP)=-3.71354692995576e-5*ZC(INP)**7+0.0008508042*ZC(INP)**6-0.0074501233*ZC(INP)**5 &
  !    +0.0287202493*ZC(INP)**4-0.0279210481*ZC(INP)**3-0.1060775951*ZC(INP)**2 &
  !    +0.1756108394*ZC(INP)+0.2104104465
  !ED(INP)=CMU75*TE(INP)**1.5/(0.07*blthick) !/(CMU*TE(INP))

!.....Bolund hill prescribed inlet:
!  U(IJK)=LOG((ZC(IJK)*120-ZC(INBC)*120)/0.0003) ! log profil case 270
!  V(IJK)=0.
!  W(IJK)=0.
!  TE(IJK)=0.928 ! vrednost data za blind comparison
!  ED(IJK)=(0.4)**3/(CAPPA*(ZC(IJK)-ZC(INBC))*120)  !CMU75*TE(IJK)**1.5/(0.4*0.1) 

!  UO(INP) = small
!  VO(INP) = small
!  WO(INP) = small

!.....log profile FOR Dobric case;
!  inbc = LI(I)+J 
!  if ( (zc(inp)-zc(inbc)) .lt. 500. ) then
!  !if(k.eq.2) print*, i, j, (zc(inp)-zc(inbc))
!    ust=0.345843
!    U(INP) = min(-(ust/cappa)*log((zc(inp)-zc(inbc))/zzero),-small)
!    TE(INP) = (ust)**2/sqrt(cappa)*(1.-((zc(inp)-zc(inbc))/500.))**2
!  else
!    U(INP) = U(INP)
!    TE(INP) = SMALL
!  endif
!  ed(inp) = max(TE(INP)**1.5/10.,small) !<-- epsilon eq.
!  !ed(inp) = ed(inp)/0.02973/TE(INP) !<-- omega eq.
!  !ed(inp)=ed(inp)/(cmu*te(inp)+small)! k-eps --> k-omega


  ! Molecular viscosity
  VIS(INP)=VISCOS
  VISOB(INP)=VISCOS
  ! Temperature
  T(INP)=TIN
  ! Temperature variance
  VART(INP)=VARTIN
  ! Concentration
  CON(INP)=CONIN

  ! Effective viscosity
  IF(LTURB) THEN
    IF(STDKEPS.or.DURBIN.or.RNG.or.REALIZABLE.or.LowRE_LB) VIS(INP)=VIS(INP)+DENSIT*TE(INP)**2*CMU/(ED(INP)+SMALL)   
    IF(Wilcox.or.SST.or.SAS.or.EARSM_WJ.or.EARSM_M) VIS(INP)=VIS(INP)+DENSIT*TE(INP)/(ED(INP)+SMALL)
  ENDIF

  ! Reynolds stress tensor components
  UU(INP) = 0.0d0
  VV(INP) = 0.0d0
  WW(INP) = 0.0d0
  UV(INP) = 0.0d0
  UW(INP) = 0.0d0
  VW(INP) = 0.0d0

  ! Turbulent heat fluxes
  UTT(INP) = 0.0d0
  VTT(INP) = 0.0d0
  WTT(INP) = 0.0d0

  ! Reynolds stress anisotropy
  IF(EARSM_WJ.OR.EARSM_M) BIJ(:,INP)=0.0d0

  END DO
  END DO
  END DO
  
!-----Field initialisation loop over inner cells--------------------------------
    
!-----Field initialisation loop over ALL cells----------------------------------
  do inp=1,nxyza
    den(inp)=densit
  enddo
!-----Field initialisation loop over ALL cells----------------------------------





! 8)  Read Restart File And Set Field Values

!
!.....Read field values - results of previous run
  if(lread) then
  call readfiles
    do inp=icst,icen
      pp(inp)=p(inp)
    end do
  end if



! 9)  Initial Gradient Calculation
  gradU=0.0D0
  gradV=0.0D0
  gradW=0.0D0
  gradP=0.0D0
  gradTE=0.0D0
  gradED=0.0D0

!.....Prepare System Matrix For Least-Squares Gradient Calculation (Allways)
  IF (LSTSQ) THEN
    ! It is done by setting this --v value to one.
    CALL GRAD_LSQ_QR(U,gradU,1)

!.....Initial Calculation Of Gradients
!   For field update set this ---v value to two allways.
    CALL GRAD_LSQ_QR(U,gradU,2)
    CALL GRAD_LSQ_QR(V,gradV,2)
    CALL GRAD_LSQ_QR(W,gradW,2)

  ELSE ! IF (GAUSS)

!.....Initial Calculation Of Gradients - Using Gauss_Grad    
    CALL GRAD_GAUSS(U,gradU(1,:),gradU(2,:),gradU(3,:))
    CALL GRAD_GAUSS(V,gradV(1,:),gradV(2,:),gradV(3,:))
    CALL GRAD_GAUSS(W,gradW(1,:),gradW(2,:),gradW(3,:))

  ENDIF

! 10) Calculate distance to the nearest wall.

  ! Source term
!  do k=2,nkm; do i=2,nim; do j=2,njm
!  inp=lk(k)+li(i)+j
!      ! Wall distance Poisson eq. source :
!      su(inp) = Vol(inp) 
!  enddo; enddo; enddo  

  ! Initialize solution and set fixedValue boundaries
  pp=0. 
  do k=1,nk
   do i=1,ni
    do j=1,nj
      inp=lk(k)+li(i)+j
      intc = lk(nk)+li(i)+j
      inbc = lk(1)+li(i)+j
      ! 1. If channel - upper and lower boundaries are Wall:
      !pp(inp) = min(zc(intc)-zc(inp),zc(inp)-zc(inbc))
      ! 2. If Boundary layer - only lower boundary is Wall:
      pp(inp) = zc(inp)-zc(inbc)
      enddo
    enddo
  enddo   

  wallDistance = pp
  pp=p

  ! Initial internal field
  !do k=2,nkm; do i=2,nim; do j=2,njm
  !inp=lk(k)+li(i)+j
  !    pp(inp) = 0.
  !enddo; enddo; enddo 
 
!  sv = 1.0d0        ! Unit coefficient array
!  call fvm_laplacian(sv,pp) ! Laplacian operator and BCs

!  call cgstab_sip(pp,ip) 

  ! Gradient of solution field stored in pp (gradient stored in gradP) :
!  call grad_lsq_qr(pp,gradp,2,d,nxyza,nx,nz, &
!           nkm,nim,njm,lk,li,xc,yc,zc,nj,1,nij)

  ! Wall distance computation from Poisson eq. solution stored in pp:
!  wallDistance = -sqrt(  gradp(1,:)*gradp(1,:)+gradp(2,:)*gradp(2,:)+gradp(3,:)*gradp(3,:)  ) + &
!          sqrt(  gradp(1,:)*gradp(1,:)+gradp(2,:)*gradp(2,:)+gradp(3,:)*gradp(3,:) + 2.*pp  )

 
!  su = 0.0d0; sv = 0.0d0 ! Clear arrays

!  call plot_3D_field_vtk (88, trim(out_folder_path)//'/wallDistance_scalar_field', 'scalar', &
!             'vtk', 'WDIS_field', 'wall-distance ',          &
!          NI, NJ, NK, 1, 1, 1, Xc, Yc, Zc, wallDistance, 0.0, 0.0)

!   Open(Unit=87,File=Trim(Out_Folder_Path)//'/wallDistance.plt') 
!   Rewind 87
!   Write(87,*) 'Title     = " "'
!   Write(87,*) 'Variables = "X"'
!   Write(87,*) '"Y"'
!   Write(87,*) '"Z"'
!   Write(87,*) '"Wdist"'
!   Write(87,*) 'Zone T=" "'
!   Write(87,*) 'I=',Ni, ' ,J=',Nj, ' ,K=',Nk,', F=Point'
!   Do k=1,nk; do j=1,nj; do i=1,ni
!   Inp=Lk(K)+Li(I)+J
!   Write(87,*) Xc(Inp),Yc(Inp),Zc(Inp),wallDistance(Inp)
!   Enddo; Enddo; Enddo 
!   Close(87)    

  RETURN
  END
