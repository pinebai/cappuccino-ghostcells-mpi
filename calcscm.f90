!***********************************************************************
!
      SUBROUTINE CALCSCM(FI,GRADFI,IFI)
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE GEOMETRY
      USE COEF
      USE COEFB
      USE VARIABLES
      USE BUOY
      USE TIME_MOD
      USE OBSTACLE
      USE GRADIENTS
      USE OMEGA_Turb_Models
      USE fieldManipulation ! Volume Weighted Average function

      IMPLICIT NONE
!
!***********************************************************************
!
      INTEGER, INTENT(IN) :: IFI
      REAL(PREC), DIMENSION(NXYZA) :: FI
      REAL(PREC), DIMENSION(3,NXYZA) :: GRADFI

!
!     LOCAL VARIABLES
!
      INTEGER ::    I, J, K, INP
      REAL(PREC) :: GAM, PRTR, APOTIME, CONST, URFRS, URFMS, &
                    UTP, VTP, WTP, UTN, VTN, WTN, &
                    DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DWDX, DWDY, DWDZ, &
                    GENP, GENN, SUT, &
                    UTTBUOY, VTTBUOY, WTTBUOY, &
                    ETARLZB,ETARNG,RETA
      REAL(PREC) :: Ssq, sqrt2r
      REAL(PREC) :: f_1, f_2, f_mu,Re_y, wdis! Za Low-Re k-epsilon modele
      REAL(PREC) :: tmp,alphast,DOMEGAP,DOMEGAN,VIST ! Za k-omega modele
      real(prec) :: lengthrke, lengthles, lengthdes, DimMax, DesDissipTke ! DES model

      sqrt2r = 1./DSQRT(2.0D0)

!.....Variable specific coefficients:
      IDIR=IFI
      GAM=GDS(IFI)
      PRTR=PRTINV(IFI)

!.....Calculate gradient: 
      if (lstsq) then
        call grad_lsq_qr(fi,gradfi,2)
      elseif (gauss) then
        call grad_gauss(fi,gradfi(1,:),gradfi(2,:),gradfi(3,:))
      endif

!
!.....CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
      IF(IFI.EQ.ITE) THEN

!.....CALCULATE PRODUCTION ... GEN(IJ)
      CALL FIND_STRAIN_RATE
!
      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM

      INP=LK(K)+LI(I)+J

      IF(LEVM) THEN
!=========================================================
!     STANDARD PRODUCTION
!=========================================================
      Ssq=STRAIN(INP)*STRAIN(INP)
      GEN(INP)=ABS(VIS(INP)-VISCOS)*Ssq


!=====PRODUCTION LIMITER FOR SST AND SAS MODELS:=======================
!.....10*bettainf=10*0.09=0.9 -> see below TODO BETTAST za Low-Re
      IF (SST.OR.SAS) THEN
!.....High-Re version.....................................................
      GEN(INP)=MIN(GEN(INP),0.9*DEN(INP)*TE(INP)*ED(INP))                !
      if (LowRe) then
!.....Low-Re version of Wilcox and SST k-omega.............................
      tmp=10.*0.09*(4./15.+(DEN(INP)*TE(INP)/(8.*VISCOS*ED(INP)))**4)  &  !
                  /(1.    +(DEN(INP)*TE(INP)/(8.*VISCOS*ED(INP)))**4)     !  
      GEN(INP)=MIN(GEN(INP),tmp*DEN(INP)*TE(INP)*ED(INP))                 !
!.........................................................................!
      end if 
      END IF

      IF(RNG) THEN
!=======================================================================
!     STRAIN = sqrt (Sij*Sij) za RNG k-epsilon ! FIXME Da li treba da se racuna GE sa ovim STRAIN!????
!     S = sqrt (2*Sij*Sij) za REALIZABLE k-epsilon !!!
!=======================================================================
       STRAIN(INP) = STRAIN(INP) * sqrt2r
!=======================================================================
      ENDIF

!
      ELSE IF(LASM) THEN
!=========================================================
!     EXACT PRODUCTION (ISKORISTI GRADIJENTE KOJE VEC IMAS DUDX = GRADU(1,INP)),...
!=========================================================

      DUDX = GRADU(1,INP)
      DUDY = GRADU(2,INP)
      DUDZ = GRADU(3,INP)

      DVDX = GRADV(1,INP)
      DVDY = GRADV(2,INP)
      DVDZ = GRADV(3,INP)

      DWDX = GRADW(1,INP)
      DWDY = GRADW(2,INP)
      DWDZ = GRADW(3,INP)

      GEN(INP)=-DEN(INP)*(UU(INP)*DUDX+UV(INP)*(DUDY+DVDX)+ &
                          UW(INP)*(DUDZ+DWDX)+VV(INP)*DVDY+ &
                          VW(INP)*(DVDZ+DWDY)+WW(INP)*DWDZ)

      END IF !![production calculation ]

      ENDDO
      ENDDO
      ENDDO

!
!----------------------------
!      CALL CALCSTRESS
!      CALL CALCHEATFLUX
!----------------------------

      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM

      INP=LK(K)+LI(I)+J
!
      GENP=MAX(GEN(INP),ZERO)
      GENN=MIN(GEN(INP),ZERO)
!
!=====================================
!.....[SOURCE TERMS]: isothermal
!=====================================
!.....ADD PRODUCTION TERM TO THE RHS (SAME FOR ALL MODELS):
      SU(INP)=GENP*VOL(INP)              

!.....ADD DESTRUCTION TERM TO THE LHS:

      IF(STDKEPS.OR.DURBIN.OR.RNG.OR.REALIZABLE) THEN
!======================================================================
      SP(INP)=ED(INP)*DEN(INP)*VOL(INP)/(TE(INP)+SMALL)

     !  Detached eddy simulation
     IF(LDES) THEN 
        ! K-eps length scale
        lengthrke = te(inp)**1.5_dp/ed(inp)

        ! dmax=max(dx,dy,dz)
        DimMax = max ( (x(inp)-x(inp-nj)), (y(inp)-y(inp-1)), (z(inp)-z(inp-nij)) ) 

        ! LES length scale, Cdes = 0.61
        lengthles = 0.61_dp*DimMax 

        ! DES lengthscale
        lengthdes = min(lengthrke,lengthles) 

        ! DES dissipation of TKE
        DesDissipTke = DEN(INP)*TE(INP)**1.5/lengthdes

        ! Add negative rhs source term to diagonal on LHS systrem matrix and divide by the unknown.
        SP(INP) = DesDissipTke*VOL(INP)/(TE(INP)+SMALL)
     ENDIF   

      ! Move negative prodction to LHS diagonal divided by the unknown.
      SP(INP)=SP(INP)-GENN*VOL(INP)/(TE(INP)+SMALL)  
!======================================================================
      ENDIF

      IF(Wilcox.OR.SST.OR.SAS.OR.EARSM_WJ.OR.EARSM_M)THEN
!======================================================================
!.....[SOURCE TERMS]: isothermal
!     Note there is possibility to add a source term to eliminate 
!     non-physical decay of turbulence variables in the freestream
!     for external aerodynamic problems
!     Reference:
!     Spalart, P. R. and Rumsey, C. L., "Effective Inflow Conditions for 
!     Turbulence Models in Aerodynamic Calculations," AIAA Journal,
!     Vol. 45, No. 10, 2007, pp. 2544 - 2553.
!======================================================================
!.....ADD SUSTAIN TERMS (ONLY FOR SST!):
!      SU(INP)=SU(INP)+0.09*TEIN*EDIN*DEN(INP)*VOL(INP)
!.....High-Re version.....................................................
      SP(INP)=BETTAST*ED(INP)*DEN(INP)*VOL(INP)    
      if(LowRe) then                                        
!.....Low-Re version of Wilcox and SST k-omega.............................
        tmp =   0.09*(4./15.+(DEN(INP)*TE(INP)/(8.*VISCOS*ED(INP)))**4) & !
                    /(1.    +(DEN(INP)*TE(INP)/(8.*VISCOS*ED(INP)))**4)   !
        SP(INP)=tmp*ED(INP)*DEN(INP)*VOL(INP)                             !           
!.........................................................................!
      endif
!.....IF PRODUCTION TERMS ARE NEGATIVE, MOVE THEM TO LHS:
      SP(INP)=SP(INP)-GENN*VOL(INP)/TE(INP)
!======================================================================
      ENDIF
!.....END: ADD DESTRUCTION TERM TO THE LHS

!
!=====================================
!.....UNSTEADY TERM
!=====================================
      IF(BDF) THEN
!=======================================================================
!    Three Level Implicit Time Integration Method:
!    in case that BTIME=0. --> Implicit Euler
!=======================================================================
      APOTIME=DEN(INP)*VOL(INP)/TIMESTEP
      SUT=APOTIME*((1+BTIME)*TEO(INP)-0.5*BTIME*TEOO(INP))
      SU(INP)=SU(INP)+SUT
      SP(INP)=SP(INP)+APOTIME*(1+0.5*BTIME)
!=======================================================================
      ENDIF
!
!=====================================
!.....[SOURCE TERMS]: buoyancy
!=====================================
      IF(LCAL(IEN).AND.LBUOY) THEN
!----
      IF(BOUSSINESQ) THEN
         UTTBUOY=-GRAVX*DEN(INP)*UTT(INP)*VOL(INP)*BETA
         VTTBUOY=-GRAVY*DEN(INP)*VTT(INP)*VOL(INP)*BETA
         WTTBUOY=-GRAVZ*DEN(INP)*WTT(INP)*VOL(INP)*BETA
      ELSE !IF(BOUSSINESQ.EQ.0) THEN
         UTTBUOY=-GRAVX*DEN(INP)*UTT(INP)*VOL(INP)/(T(INP)+273.)
         VTTBUOY=-GRAVY*DEN(INP)*VTT(INP)*VOL(INP)/(T(INP)+273.)
         WTTBUOY=-GRAVZ*DEN(INP)*WTT(INP)*VOL(INP)/(T(INP)+273.)
      END IF
!----
      UTP=MAX(UTTBUOY,ZERO)
      VTP=MAX(VTTBUOY,ZERO)
      WTP=MAX(WTTBUOY,ZERO)
      UTN=MIN(UTTBUOY,ZERO)
      VTN=MIN(VTTBUOY,ZERO)
      WTN=MIN(WTTBUOY,ZERO)
!----
      SU(INP)=SU(INP)+UTP+VTP+WTP
      SP(INP)=SP(INP)-(UTN+VTN+WTN)/(TE(INP)+SMALL)
!----
      END IF

!.....End of ITE volume source terms
      ENDDO
      ENDDO
      ENDDO

!****************************************
      ELSEIF(IFI.EQ.IED) THEN
!****************************************

!      magStrain =  volumeWeightedAverage(gradED)
!      write(66,*)'Volume weighted average of Dissipation Grad.',magStrain

!======================================================================
!.....Terms for Menter SST model:
      IF (SST.OR.SAS) call calculate_stuff_for_sst
!.....Terms for EARSM_WJ and EARSM_M model:
      IF (EARSM_WJ.OR.EARSM_M) call calculate_stuff_for_earsm
!======================================================================

      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM

      INP=LK(K)+LI(I)+J
!
      GENP=MAX(GEN(INP),ZERO)
      GENN=MIN(GEN(INP),ZERO)
!
!===============================================
!.....[Hybrid ALPHA model: Ce2=Ce1+0.48/alpha
!===============================================
      IF (ALPHAMODEL.and.DURBIN) THEN
      AL_RANS(INP)=TE(INP)**(3./2.)/(ED(INP)+SMALL)
      AL_LES(INP)=VOL(INP)**(1./3.)
      ALPH(INP)=MAX(1.d0,AL_RANS(INP)/AL_LES(INP))
      C2=C1+0.48/ALPH(INP)
      ENDIF

!======================================================================
!.....[Hybrid ALPHA model: Ce2=Ce1+0.48/alpha, samo k-omega sst verzija!
!     EXPERIMENTALNA VARIJANTA!
!     Ovde je DIFF je funkcija polozaja i ekvivalent je 0.48 kod k-eps
!      0.48=c2-c1 => u sst modelu: DIFF(inp)=BETTASST(inp)-ALPHASST(inp)
!     Aktivira se tako sto SST=True and SAS=True.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(SST.and.ALPHAMODEL) THEN
      AL_RANS(INP)=TE(INP)**(1/2.)/(ED(INP)*CMU)
      AL_LES(INP)=VOL(INP)**(1/3.)
      !AL_KOLM(INP)=((VISCOS**3.)/(ED(INP)*TE(INP)*CMU))**(1/4.)
      ALPH(INP)=MAX(1.d0,AL_RANS(INP)/AL_LES(INP))
      !ALPH(INP)=MAX(1.d0,0.5*AL_RANS(INP)/AL_LES(INP))
      BETTASST(INP)=ALPHASST(INP)+DIFF(INP)/ALPH(INP)
!======================================================================
      ENDIF


!=====CROSS DIFFUSION FOR SST MODEL====================================
      !IF (SST.OR.SAS.OR.EARSM_WJ.OR.EARSM_M) THEN
        DOMEGAP=MAX(DOMEGA(INP),ZERO)
        DOMEGAN=MIN(DOMEGA(INP),ZERO)
      !END IF
!======================================================================

!
!=====================================
!.....[SOURCE TERMS]: isothermal
!=====================================
      SV(INP)=0.0d0 ! This is why we use Sv in bcscalarm.f!!! So there we add only GENT calculated in wallbc.f

!
!=====Standard k-epsilon or KEPS+Durbin=============================
      IF (STDKEPS.OR.DURBIN) THEN
      SU(INP)=C1*GENP*ED(INP)*VOL(INP)/(TE(INP)+SMALL)
      SP(INP)=C2*DEN(INP)*ED(INP)*VOL(INP)/(TE(INP)+SMALL)
      SP(INP)=SP(INP)-C1*GENN*VOL(INP)/(TE(INP)+SMALL)
!=====END: Standard k-epsilon or KEPS+Durbin=============================
      ENDIF 

!=====REALIZABLE k-epsilon==============================================
      IF (REALIZABLE) THEN

      GENP=MAX(STRAIN(INP),ZERO)
      GENN=MIN(STRAIN(INP),ZERO)

      ETARLZB = STRAIN(INP)*TE(INP)/(ED(INP)+SMALL)
      C1 = MAX(0.43,ETARLZB/(ETARLZB+5.))
      SU(INP)=C1*GENP*ED(INP)*VOL(INP)
      SP(INP)=C2*DEN(INP)*ED(INP)*VOL(INP)/ &
                               ( TE(INP)+SQRT(VISCOS/DENSIT*ED(INP)) )
      SP(INP) = SP(INP) - C1*GENN*ED(INP)*VOL(INP)
!=====END:REALIZABLE k-epsilon==========================================
      ENDIF
!
!=====RENORMALIZATION GROUP (RNG) k-epsilon=============================
      IF (RNG) THEN
      SU(INP)=C1*GENP*ED(INP)*VOL(INP)/(TE(INP)+SMALL)
      SP(INP)=C2*DEN(INP)*ED(INP)*VOL(INP)/(TE(INP)+SMALL)
      SP(INP)=SP(INP)-C1*GENN*VOL(INP)/(TE(INP)+SMALL)

      ETARNG = STRAIN(INP)*TE(INP)/(ED(INP)+SMALL)
      RETA = CMU*ETARNG**3*(1-ETARNG/4.38)/(1+0.012*ETARNG**3)
      IF (ETARNG.LE.4.38d0) THEN
      SP(INP)=SP(INP)+RETA*DEN(INP)*ED(INP)*VOL(INP)/(TE(INP)+SMALL)
      ELSE
      SU(INP)=SU(INP)-RETA*DEN(INP)*ED(INP)*VOL(INP)/(TE(INP)+SMALL)
      END IF
!=====END:RENORMALIZATION GROUP (RNG) k-epsilon=========================
      ENDIF !RNG

!
!=====Low-Re Lam-Bremhorst k-epsilon====================================
      IF (LowRe_LB) THEN
      !
      ! Lam-Bremhorst :
      !

      wdis = wallDistance(inp)
      Re_y = DEN(INP)*sqrt(TE(INP))*wdis/(VISCOS+SMALL)
      RET(INP)=DEN(INP)*TE(INP)**2/(VISCOS*ED(INP)+SMALL) 
      f_mu = (1.-exp(-0.0165*Re_y))**2*(1.+20.5/RET(INP))
      f_1 = 1.+(0.05/f_mu)**3
      f_2 = 1.-exp(-RET(INP)**2) 
      ! 
      ! Launder-Sharma :
      !
      !RET(INP)=DEN(INP)*TE(INP)**2/(VISCOS*ED(INP)+SMALL) 
      !f_1 = 1.0d0
      !f_2 = 1.0d0-0.3d0*exp(-Ret(inp)**2)
      !
      ! For all:
      SU(INP) = f_1*C1*GENP*ED(INP)*VOL(INP)/(TE(INP)+SMALL)
      SP(INP) = f_2*C2*DEN(INP)*ED(INP)*VOL(INP)/(TE(INP)+SMALL)
      SP(INP) = SP(INP)-f_1*C1*GENN*VOL(INP)/(TE(INP)+SMALL)
!.....Za Low-Re modele - premesti f_mu u RET - da znas u MODVIS kada vidis RET to je zapravo f_mu  
      RET(INP) = f_mu
!=====END:Low-Re Lam-Bremhorst k-epsilon================================
      ENDIF ! Low-Re Lam-Bremhorst  


!
      IF (Wilcox) THEN
!=====Wilcox  k-omega===================================================
      SU(INP)=ALPHA*GENP*ED(INP)*VOL(INP)/TE(INP)
      if (LowRe) then
!.....Low-Re version......................................................
!.....Let's find alpha*                                                  !
      alphast=(0.024+(DENSIT*TE(INP))/(6.*VISCOS*ED(INP)))  &            !
             /(1.+(DENSIT*TE(INP))/(6.*VISCOS*ED(INP)))                  !
      tmp=ALPHA/alphast*(1./9.+(DENSIT*TE(INP))/(2.95*VISCOS*ED(INP))) & !
                        /((1.+(DENSIT*TE(INP))/(2.95*VISCOS*ED(INP))))   !
      SU(INP)=tmp*GENP*ED(INP)*VOL(INP)/TE(INP)                          !
!........................................................................!
      endif

!.....ADD DESTRUCTION TERM TO THE LHS:
      SP(INP)=BETTA*DEN(INP)*ED(INP)*VOL(INP) 
      SP(INP)=SP(INP)-ALPHA*GENN*VOL(INP)/(TE(INP)+SMALL)
!=====END: Wilcox  k-omega===============================================
      END IF
      
!=====Menter k-omega SST and SAS ========================================
!     Ovde ce sad uzeti limitiran GEN(INP) jer je gore usao u IF SST
!     ovde \alpha nije konstanta kao u STD vec se racuna u "calc_stuff_for_sst",
!     takodje ovde deli produkciju sa \nu_t umesto da mnozi sa \omega/k.
      IF (SST) THEN
      VIST = (VIS(INP)-VISCOS)/DENSIT
      SU(INP)=ALPHASST(INP)*GENP*VOL(INP)/VIST
      SU(INP)=SU(INP)+DOMEGAP*VOL(INP)
!.....ADD SUSTAIN TERMS
!      SU(INP)=SU(INP)+BETTASST(INP)*EDIN*EDIN*DEN(INP)*VOL(INP)
!.....ADD DESTRUCTION TERM TO THE LHS:
      VIST = (VIS(INP)-VISCOS)/DENSIT
      SP(INP)=BETTASST(INP)*DEN(INP)*ED(INP)*VOL(INP) 
      SP(INP)=SP(INP)-ALPHASST(INP)*GENN*VOL(INP)  &
                       /(VIST*ED(INP))
!      SP(INP)=SP(INP)-DOMEGAN*VOL(INP)/ED(INP)
!=====END: Menter k-omega SST ============================================
      END IF

!.....SAS ONLY. ADDS QSAS PRODUCTION
      IF (SAS) THEN
!=====SAS only ===============================================
      VIST = (VIS(INP)-VISCOS)/DEN(INP)
      SU(INP)=ALPHASST(INP)*GENP*VOL(INP)/VIST
      SU(INP)=SU(INP)+DOMEGAP*VOL(INP)
!.....ADD SUSTAIN TERM
!      SU(INP)=SU(INP)+BETTASST(INP)*EDIN*EDIN*DEN(INP)*VOL(INP)
!.....ADD SAS PRODUCTION TERM (!) :
      SU(INP)=SU(INP)+QSAS(INP)*VOL(INP)
!.....ADD DESTRUCTION TERM TO THE LHS:
      VIST = (VIS(INP)-VISCOS)/DENSIT
      SP(INP)=BETTASST(INP)*DEN(INP)*ED(INP)*VOL(INP) 
      SP(INP)=SP(INP)-ALPHASST(INP)*GENN*VOL(INP)  &
                       /(VIST*ED(INP))
!      SP(INP)=SP(INP)-DOMEGAN*VOL(INP)/ED(INP)
!=====END: SAS only ===============================================
      END IF


      IF (EARSM_WJ.OR.EARSM_M) THEN
!=====EARSM_WJ and EARSM_M models standalone or in SAS hybrid version ==========
      SU(INP)=ALPHASST(INP)*GENP*ED(INP)*VOL(INP)/(TE(INP)+SMALL)   !earsm
      SU(INP)=SU(INP)+DOMEGAP*VOL(INP)
!.....ADD SUSTAIN TERMS
!      SU(INP)=SU(INP)+BETTASST(INP)*EDIN*EDIN*DEN(INP)*VOL(INP)

!.....SAS ALSO ADDS QSAS PRODUCTION
      IF (SAS) THEN
!.....ADD SAS PRODUCTION TERM
      SU(INP)=SU(INP)+QSAS(INP)*VOL(INP)
      END IF

!.....ADD DESTRUCTION TERM TO THE LHS:
      SP(INP)=BETTASST(INP)*DEN(INP)*ED(INP)*VOL(INP) 
      SP(INP)=SP(INP)-ALPHASST(INP)*GENN*VOL(INP) &
                       /(TE(INP)+SMALL)              !earsm
      SP(INP)=SP(INP)-DOMEGAN*VOL(INP)/ED(INP)
!=====END: EARSM_WJ and EARSM_M models standalone or in SAS hybrid version =====
      END IF

!
!=====================================
!.....UNSTEADY TERM
!=====================================
      IF(BDF) THEN
!=======================================================================
!    Three Level Implicit Time Integration Method:
!    in case that BTIME=0. --> Implicit Euler
!=======================================================================
      APOTIME=DEN(INP)*VOL(INP)/TIMESTEP
      SUT=APOTIME*((1+BTIME)*EDO(INP)-0.5*BTIME*EDOO(INP))
      SU(INP)=SU(INP)+SUT
      SP(INP)=SP(INP)+APOTIME*(1+0.5*BTIME)
!=======================================================================
      ENDIF
!
!=====================================
!.....[SOURCE TERMS]: buoyancy
!=====================================
      IF(LCAL(IEN).AND.LBUOY) THEN
      CONST=C3*DEN(INP)*ED(INP)*VOL(INP)/(TE(INP)+SMALL)
!----
      IF(BOUSSINESQ) THEN
         UTTBUOY=-GRAVX*UTT(INP)*CONST*BETA
         VTTBUOY=-GRAVY*VTT(INP)*CONST*BETA
         WTTBUOY=-GRAVZ*WTT(INP)*CONST*BETA
      ELSE !IF(BOUSSINESQ.EQ.0) THEN
         UTTBUOY=-GRAVX*UTT(INP)*CONST/(T(INP)+273.)
         VTTBUOY=-GRAVY*VTT(INP)*CONST/(T(INP)+273.)
         WTTBUOY=-GRAVZ*WTT(INP)*CONST/(T(INP)+273.)
      END IF
!----
      UTP=MAX(UTTBUOY,ZERO)
      VTP=MAX(VTTBUOY,ZERO)
      WTP=MAX(WTTBUOY,ZERO)
      UTN=MIN(UTTBUOY,ZERO)
      VTN=MIN(VTTBUOY,ZERO)
      WTN=MIN(WTTBUOY,ZERO)
!----
      SU(INP)=SU(INP)+UTP+VTP+WTP
      SP(INP)=SP(INP)-(UTN+VTN+WTN)/(ED(INP)+SMALL)
!----
      END IF

!.....End of IED volume source terms
      ENDDO
      ENDDO
      ENDDO
!--------------------------------------
      END IF

!
!.....CALCULATE TERMS INTEGRATED OVER SURFACES
!.....ONLY INNER SURFACES
!.....EAST CELL - FACE
      CALL FLUXSCM(NJ,1,NIJ,FI,GRADFI,IFI, &
                   AR1X,AR1Y,AR1Z, &
                   FX,AE,AW,F1)

!.....NORTH CELL - FACE
      CALL FLUXSCM(1,NIJ,NJ,FI,GRADFI,IFI, &
                   AR2X,AR2Y,AR2Z, &
                   FY,AN,AS,F2)

!.....TOP   CELL - FACE
      CALL FLUXSCM(NIJ,NJ,1,FI,GRADFI,IFI, &
                   AR3X,AR3Y,AR3Z, &
                   FZ,AT,AB,F3)

      URFRS=URFR(IFI)
      URFMS=URFM(IFI)

      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM

      INP=LK(K)+LI(I)+J

      IF(CN) THEN
!.....Crank-Nicolson stuff - only once:
      AE(INP)=0.5D0*AE(INP)
      AN(INP)=0.5D0*AN(INP)
      AT(INP)=0.5D0*AT(INP)
      AW(INP)=0.5D0*AW(INP)
      AS(INP)=0.5D0*AS(INP)
      AB(INP)=0.5D0*AB(INP)
!.....Crank-Nicolson time stepping source terms
      APOTIME=DEN(INP)*VOL(INP)/TIMESTEP
      IF(IFI.eq.ITE) THEN
      SU(INP)=SU(INP)+(AE(INP)*TEO(INP+NJ)  + AW(INP)*TEO(INP-NJ)+    &
                       AN(INP)*TEO(INP+1)   + AS(INP)*TEO(INP-1)+     &
                       AT(INP)*TEO(INP+NIJ) + AB(INP)*TEO(INP-NIJ))+  &
              (APOTIME-AE(INP)-AW(INP)                                &
                      -AN(INP)-AS(INP)                                &
                      -AT(INP)-AB(INP))*TEO(INP) 
      ELSE ! IFI==IED   
      SU(INP)=SU(INP)+(AE(INP)*EDO(INP+NJ)  + AW(INP)*EDO(INP-NJ)+    &
                       AN(INP)*EDO(INP+1)   + AS(INP)*EDO(INP-1)+     &
                       AT(INP)*EDO(INP+NIJ) + AB(INP)*EDO(INP-NIJ))+  &
              (APOTIME-AE(INP)-AW(INP)                                &
                      -AN(INP)-AS(INP)                                &
                      -AT(INP)-AB(INP))*EDO(INP)      
      ENDIF ! Checking if this is tke or epsilon field      
      SP(INP)=SP(INP)+APOTIME
!.....End of Crank-Nicolson time stepping source terms
!.....End of Crank-Nicolson stuff:
      ENDIF

!.....Main diagonal term assembly
      AP(INP)=AE(INP)+AW(INP)+AN(INP)+AS(INP)+AT(INP)+AB(INP)+SP(INP)
!.....Underelaxation:
      AP(INP)=AP(INP)*URFRS
      SU(INP)=SU(INP)+URFMS*AP(INP)*FI(INP)
      ENDDO
      ENDDO
      ENDDO
!
!.....SOLVING LINEAR SYSTEM:
      CALL SIPSOL(FI,IFI)
      !CALL CGSTAB_SIP(FI,IFI)

!.....These field values cannot be negative
      IF(IFI.EQ.ITE.OR.IFI.EQ.IED) THEN
         DO INP=ICST,ICEN
         FI(INP)=MAX(FI(INP),SMALL)
         ENDDO
      ENDIF

      RETURN
      END
