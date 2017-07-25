!***********************************************************************
!
      SUBROUTINE calculate_stuff_for_earsm
!
!***********************************************************************
!     Adopting:
!     A subroutine to calculate the extra anisotropy components and
!     the effective eddy-viscosity coefficient used with the EARSM.
!     Date:      Author:     Affiliation:
!     23.11.2001 Ville H. /  Laboratory of Aerodynamics, Espoo, Finland 
!     13.12.2011 Nikola M. / Laboratory for Thermal Engineering and Energy,
!                            Institute of Nuclear Sciences "Vinca", Belgrade, Serbia  
!
!     There are two versions here: the original Wallin-Johansson and Menter et. al 2009
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
      USE BOUNDC
      USE OBSTACLE
      USE INLET
      USE GRADIENTS
      USE OMEGA_Turb_Models
      USE fieldManipulation

      IMPLICIT NONE
!
!***********************************************************************
!
      INTEGER :: I,J,K,INP,INB
      REAL(PREC) :: P3,P6,TT,C1E,C1P,C1P3,C1PSQ
      REAL(PREC) :: CT,TSTUR,TSVIS,TTS,HTTS                     ! Turbulence timescale related
      REAL(PREC) :: S11,S12,S13,S22,S23,S33,S21,S32,W12,W13,W23 ! Components of strain and vorticity tensors
      REAL(PREC) :: S11SQ,S22SQ,S33SQ,S12SQ,S13SQ,S23SQ         ! Strain rate tensor components squared
      REAL(PREC) :: W12SQ,W13SQ,W23SQ                           ! Vorticity tensor components squared
      REAL(PREC) :: SII,WII,WIIP3,SWWIV,SWWIVTT,SSWWV
      REAL(PREC) :: TERM3C11,TERM3C12,TERM3C13,TERM3C22,TERM3C23
      REAL(PREC) :: TERM4C11,TERM4C12,TERM4C13,TERM4C22,TERM4C23
      REAL(PREC) :: TERM6C11,TERM6C12,TERM6C13,TERM6C22,TERM6C23
      REAL(PREC) :: TERM9C11,TERM9C12,TERM9C13,TERM9C22,TERM9C23
      REAL(PREC) :: P1,P2,SQP2,PM,PMP,SQPM,FACOS,FNC,DE,FII1,FII2,FN,FNSQ,Q,PQ,PQ1,Q1
      REAL(PREC) :: BETA1,BETA3,BETA4,BETA6,BETA9
      REAL(PREC) :: TERMLR,TERMLR11,TERMLR12,TERMLR13,TERMLR22,TERMLR23
      !REAL(PREC) :: beta1eq

      REAL(PREC) :: WLDIST,DOMEGAPL,KSI,SIGMAD,F    ! SST related
      REAL(PREC) :: DELTA                                       ! SST and SAS Related
      REAL(PREC) :: DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ
      REAL(PREC) :: DTEDX,DTEDY,DTEDZ,DEDDX,DEDDY,DEDDZ
      REAL(PREC) :: FAC,DXF,DYF,DZF,ALF,BET,GAM,ARE,ARER,REYLR,FACLR,SIIEQ
      REAL(PREC) :: LMT,T1,T2,TMP1,TMP2

!     Often used numbers
      P3= 1./3.   ! One third
      P6= 0.5*P3  ! One sixth
      TT= 2.0*P3  ! Two thirds

!     Model coefficient and its variants
      C1E = 1.8
      C1P = 9.0*(C1E - 1.0)/4.0 ! C1_prime
      C1P3 = C1P*P3
      C1PSQ = C1P**2

!.....If needed get the Von Karman length scale for SAS model.
      if(SAS) lvk(:) = von_karman_lengthscale()

!=====Velocity Gradient Tensor computation=============================

!.....FIND GRADIENTS OF VELOCITIES U,V, AND W

      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM

      INP=LK(K)+LI(I)+J
      INB=(I-1)*NJ+J

!.....Turbulence Timescale
      CT = 6.0                                          ! Model coefficient
      TSTUR = 1./(CMU*ED(INP))                          ! BETTAST=0.09
      TSVIS = CT*SQRT(          VISCOS                  &
                    /( DENSIT*CMU*ED(INP)*TE(INP) )  )
      TTS = MAX(TSTUR,TSVIS)                            ! Limiter for turbulence time-scale

!.....Half of the time scale
      HTTS = 0.5*TTS

!
!--------------------------------
!      [VELOCITY GRADIENTS: ]
!--------------------------------
      DUDX = gradU(1,INP)
      DUDY = gradU(2,INP)
      DUDZ = gradU(3,INP)
      
      DVDX = gradV(1,INP)
      DVDY = gradV(2,INP)
      DVDZ = gradV(3,INP)

      DWDX = gradW(1,INP)
      DWDY = gradW(2,INP)
      DWDZ = gradW(3,INP)



!.....FIND SCALAR INVARIANT OF THE STRAIN RATE AND VORTICITY TENSOR
!     Found it in find_strain_rate subroutine
!      SIJMOD(INP) >>>       STRAIN(INP)
!      OMIJMOD(INP) = DSQRT(W12**2 + W13**2 + W23**2)

!     Strain-rate and vorticity tensor components
      S11=TTS*dUdx
      S22=TTS*dVdy
      S33=TTS*dWdz

      S12=HTTS*(dUdy + dVdx)
      S13=HTTS*(dUdz + dWdx)
      S23=HTTS*(dVdz + dWdy)
 
      S21=S12     
      S32=S23

      W12=HTTS*(dUdy - dVdx)
      W13=HTTS*(dUdz - dWdx)
      W23=HTTS*(dVdz - dWdy)

!     Squares of strain-rate and vorticity tensor components
      S11SQ = S11*S11
      S22SQ = S22*S22
      S33SQ = S33*S33

      S12SQ = S12*S12
      S13SQ = S13*S13
      S23SQ = S23*S23  

      W12SQ = W12*W12
      W13SQ = W13*W13
      W23SQ = W23*W23

!     Second invariants of the strain rate and vorticity tensors
      SII = S11SQ + S22SQ + S33SQ + 2.0*(S12SQ+S13SQ+S23SQ)
      WII =-2.0*(W12SQ + W13SQ + W23SQ)
      WIIP3 = P3*WII ! One third of the invariant

!     Third invanriant of the strain rate and vorticity tensors:
      SWWIV =-S11*(W12SQ + W13SQ) - S22*(W12SQ + W23SQ)     &
            - S33*(W13SQ + W23SQ)                           &
            + 2.0*(-S12*W13*W23 + S13*W12*W23 - S23*W12*W13)
      SWWIVTT= TT*SWWIV ! Two thirds of the invariant

!     Fourth invariant of the strain rate and vorticity tensors
      SSWWV = 2.0*(-(S12*S13 + S22*S23 + S23*S33)*W12*W13    &
                  +(S11*S13 + S12*S23 + S13*S33)*W12*W23     &
                  -(S11*S12 + S12*S22 + S13*S23)*W13*W23)    &
        - (S11SQ + S12SQ + S13SQ)*(W12SQ + W13SQ)            &
        - (S12SQ + S22SQ + S23SQ)*(W12SQ + W23SQ)            &
        - (S13SQ + S23SQ + S33SQ)*(W13SQ + W23SQ)

!     Tensor component terms for beta 3
      TERM3C11 = - W12SQ - W13SQ - WIIP3
      TERM3C22 = - W12SQ - W23SQ - WIIP3
      TERM3C12 = - W13*W23
      TERM3C13 =   W12*W23
      TERM3C23 = - W12*W13

!     Tensor component terms for beta 4
      TERM4C11 =-2.0*(S12*W12 + S13*W13)
      TERM4C22 = 2.0*(S12*W12 - S23*W23)
      TERM4C12 = (S11-S22)*W12       - S23*W13       - S13*W23
      TERM4C13 =     - S23*W12 + (S11-S33)*W13       + S12*W23
      TERM4C23 =       S13*W12       + S12*W13 + (S22-S33)*W23

!     tensor component terms for beta 6
      TERM6C11 = -2.0*((S12*W13 - S13*W12)*W23       &
                     + S11*(W12SQ + W13SQ)) - SWWIVTT- WII*S11
      TERM6C22 = -2.0*((S23*W12 + S12*W23)*W13       &
                     + S22*(W12SQ + W23SQ)) - SWWIVTT- WII*S22
      TERM6C12 = -S12*(2.0*W12SQ + W13SQ + W23SQ)    &
               - (S13*W13-S23*W23)*W12 - (S11+S22)*W13*W23- WII*S12
      TERM6C13 = -S13*(W12SQ + 2.0*W13SQ + W23SQ)    &
               - (S12*W12+S23*W23)*W13 + (S11+S33)*W12*W23- WII*S13
      TERM6C23 = -S23*(W12SQ + W13SQ + 2.0*W23SQ)    &
               + (S12*W12-S13*W13)*W23 - (S22+S33)*W12*W13- WII*S23

!     Tensor component terms for beta 9
      TERM9C11 =-2.0*(( S12*W12 + S13*W13 - S23*W23)*W12SQ  &
                     +( S12*W12 + S13*W13 + S23*W23)*W13SQ  &
                     +( S22-S33)*W12*W13*W23)
      TERM9C22 =-2.0*((-S12*W12 - S13*W13 + S23*W23)*W12SQ  &
                     +(-S12*W12 + S13*W13 + S23*W23)*W23SQ  &
                     +(-S11+S33)*W12*W13*W23)
      TERM9C12 = ((S11-S22)*W12 - 2.0*(S13*W23+S23*W13))*W12SQ  &
               + ((S11-S33)*W12 - 2.0* S13*W23)         *W13SQ  &
               + ((S33-S22)*W12 - 2.0* S23*W13)         *W23SQ
      TERM9C13 = ((S11-S22)*W13 + 2.0* S12*W23)         *W12SQ  &
               + ((S11-S33)*W13 + 2.0*(S12*W23-S23*W12))*W13SQ  &
               + ((S22-S33)*W13 - 2.0* S23*W12)         *W23SQ
      TERM9C23 = ((S22-S11)*W23 + 2.0* S12*W13)         *W12SQ  &
               + ((S11-S33)*W23 + 2.0* S13*W12)         *W13SQ  &
               + ((S22-S33)*W23 + 2.0*(S12*W13+S13*W12))*W23SQ

!     Solution of the third degree equation for N_c
!.....Correction to C1p (A3') appearing in Anti's thesis.......
!      beta1eq = -6./5.* (4.05/(16.4025 -2.*WII))              !
!      C1P = C1P + 9./4. * 2.2 * max(1+beta1eq*SII,0.)         !
!      C1PSQ = C1P**2                                          !
!      C1P3 = C1P*P3                                           !
!.............................................................!
      P1     = (C1PSQ/27.0 + 0.45*SII - TT*WII)*C1P
      P2     = P1**2 - (C1PSQ*P3**2 + 0.9*SII + TT*WII)**3
      IF (P2 .GE. 0.0) THEN
        SQP2 = SQRT(P2)
        PM = P1 - SQP2
        PMP = ABS(PM)**P3
        FNC = C1P3 + (P1+SQP2)**P3 + SIGN(PMP,PM)
      ELSE
        PM = P1**2 - P2
        SQPM = SQRT(PM)
        FACOS = P3*ACOS(P1/SQPM)
        FNC = C1P3 + 2.0*(PM**P6)*COS(FACOS)
      END IF

!.....Improvement of the approximation of the N...........................
!     Nonlinear EARSMko2005 Model with better approximation for 3D Flows !
      DE = 20.0*(FNC**4)*(FNC - 0.5*C1P)     &
        - WII*(10.0*FNC + 15.0*C1P)*FNC**2  &
        + 10.0*C1P*WII**2
      FII1 = SWWIV**2
      FII2 = SSWWV - 0.5*SII*WII
      FN = FNC + 162.0*(FII1 + FII2*FNC**2)/DE
!........................................................................!

!     The denominator of the betas
      FNSQ = FN**2
      Q = 5.0*P6*(FNSQ - 2.0*WII)*(2.0*FNSQ - WII)     !%$ <<< Wallin&Johansson A.Hellsten
      IF (EARSM_M) THEN
        Q = (1./1.245)*(FNSQ - 2.0*WII)               ! Menter 2009. 4% encrease
        Q1=Q*P6*(2.0*FNSQ - WII)                      !
        PQ1=1.0/Q1                                    !
      ENDIF
      PQ = 1.0/Q

!     Coefficients (betas) (WJ Original)
      BETA1 = - PQ*FN*(2.0*FNSQ - 7.0*WII)  !%$ <<< Wallin&Johansson A.Hellsten
      BETA3 = - PQ*12.0*SWWIV/FN            !%$ <<<
      BETA4 = - PQ*2.0*(FNSQ - 2.0*WII)     !%$ <<<
      BETA6 = - PQ*6.0*FN                   !%$ <<<
      BETA9 =   PQ*6.0                      !%$ <<<

!     Coefficients (betas) (Menter 2009)
      IF (EARSM_M) THEN
      BETA1 = - PQ*FN
      BETA3 = - PQ1*2.0*SWWIV/FN
      BETA4 = - PQ
      BETA6 = - PQ1*FN  
      BETA9 =   PQ1  
      ENDIF

!=====Low Re version===================================================
      IF (LowRe) THEN
!.....FIRST VECTOR-DISTANCE FROM GHOST CELL CENTER TO THIS CELL'S CENTER  
      FAC= -1.                                                 
      DXF=FAC*(XC(INB)-XC(INP))                                
      DYF=FAC*(YC(INB)-YC(INP))    
      DZF=FAC*(ZC(INB)-ZC(INP))
!.....SECOND X THIRD DEFINE  ALWAYS CF AREA
      ARE=DSQRT(AR1X(INB)**2+AR1Y(INB)**2+AR1Z(INB)**2)
!.....COMPONENTS OF THE UNIT NORMAL VECTOR
      ARER=1./ARE
      ALF=AR1X(INB)*ARER
      BET=AR1Y(INB)*ARER
      GAM=AR1Z(INB)*ARER
!.....NORMAL DISTANCE FROM SCALAR PRODUCT
      DELN=(DXF*ALF+DYF*BET+DZF*GAM)

!     Reynolds No. - Rey
      REYLR = Densit*sqrt(TE(inp))*DELN/Viscos             

!     factor for Low-Re(LR) correction - f1
      FACLR = 1.-exp(-0.092*sqrt(REYLR) - 1.2e-4*REYLR**2) 

!     405.*1.8**2/(216.*1.8-160.)
      SIIEQ = 5.735139863                                  

!     Additional term in anisotropy tensor expression 
      TERMLR = (1.-FACLR**2)*1.4/max(SII,SIIEQ)
      TERMLR11 = TERMLR * ((S11SQ+S12SQ+S13SQ)     - SII*P3)
      TERMLR22 = TERMLR * ((S21*S12+S22SQ+S23*S32) - SII*P3)
      TERMLR12 = TERMLR * (S11*S12+S12*S22+S13*S32)
      TERMLR13 = TERMLR * (S11*S13+S12*S23+S13*S33)
      TERMLR23 = TERMLR * (S21*S13+S22*S23+S23*S33)

!     Extra anisotropy components. Note that we use b_ij and Wallin
!     uses a_ij, which is a_ij = 2*b_ij.
      BIJ(1,INP) = 0.5*(TERMLR11                                     &           
           + FACLR**2*BETA3*TERM3C11                                 & 
           +(FACLR**2*BETA4-(1.-FACLR)*0.9/max(SII,SIIEQ))*TERM4C11  &
           + FACLR*   BETA6*TERM6C11                                 &
           + FACLR**2*BETA9*TERM9C11)
      BIJ(2,INP) = 0.5*(TERMLR12                                     &
           + FACLR**2*BETA3*TERM3C12                                 &
           +(FACLR**2*BETA4-(1.-FACLR)*0.9/max(SII,SIIEQ))*TERM4C12  &
           + FACLR*   BETA6*TERM6C12                                 &
           + FACLR**2*BETA9*TERM9C12)
      BIJ(3,INP) = 0.5*(TERMLR13                                     &
           + FACLR**2*BETA3*TERM3C13                                 &
           +(FACLR**2*BETA4-(1.-FACLR)*0.9/max(SII,SIIEQ))*TERM4C13  &
           + FACLR*   BETA6*TERM6C13                                 &
           + FACLR**2*BETA9*TERM9C13)
      BIJ(4,INP) = 0.5*(TERMLR22                                     &
           + FACLR**2*BETA3*TERM3C22                                 &
           +(FACLR**2*BETA4-(1.-FACLR)*0.9/max(SII,SIIEQ))*TERM4C22  &
           + FACLR*   BETA6*TERM6C22                                 &
           + FACLR**2*BETA9*TERM9C22)
      BIJ(5,INP) = 0.5*(TERMLR23                                     &
           + FACLR**2*BETA3*TERM3C23                                 &
           +(FACLR**2*BETA4-(1.-FACLR)*0.9/max(SII,SIIEQ))*TERM4C23  &
           + FACLR*   BETA6*TERM6C23                                 &
           + FACLR**2*BETA9*TERM9C23)
!=====End Of Low Re Version============================================
      ELSEIF (EARSM_WJ) THEN 
!     High-Re form:
!     Extra anisotropy components. Note that we use b_ij and Wallin
!     uses a_ij, which is a_ij = 2*b_ij.
      BIJ(1,INP) = 0.5*(BETA3*TERM3C11 + BETA4*TERM4C11  &
                      + BETA6*TERM6C11 + BETA9*TERM9C11)
      BIJ(2,INP) = 0.5*(BETA3*TERM3C12 + BETA4*TERM4C12  &
                      + BETA6*TERM6C12 + BETA9*TERM9C12)
      BIJ(3,INP) = 0.5*(BETA3*TERM3C13 + BETA4*TERM4C13  &
                      + BETA6*TERM6C13 + BETA9*TERM9C13)
      BIJ(4,INP) = 0.5*(BETA3*TERM3C22 + BETA4*TERM4C22  &
                      + BETA6*TERM6C22 + BETA9*TERM9C22)
      BIJ(5,INP) = 0.5*(BETA3*TERM3C23 + BETA4*TERM4C23  &
                      + BETA6*TERM6C23 + BETA9*TERM9C23)
      ELSEIF (EARSM_M) THEN 
      BIJ(1,INP) = 0.5*(BETA1*S11                        &
                      + BETA3*TERM3C11 + BETA4*TERM4C11  &
                      + BETA6*TERM6C11 + BETA9*TERM9C11)
      BIJ(2,INP) = 0.5*(BETA1*S12                        &
                      + BETA3*TERM3C12 + BETA4*TERM4C12  &
                      + BETA6*TERM6C12 + BETA9*TERM9C12)
      BIJ(3,INP) = 0.5*(BETA1*S13                        &
                      + BETA3*TERM3C13 + BETA4*TERM4C13  &
                      + BETA6*TERM6C13 + BETA9*TERM9C13)
      BIJ(4,INP) = 0.5*(BETA1*S22                        &
                      + BETA3*TERM3C22 + BETA4*TERM4C22  &
                      + BETA6*TERM6C22 + BETA9*TERM9C22)
      BIJ(5,INP) = 0.5*(BETA1*S23                        &
                      + BETA3*TERM3C23 + BETA4*TERM4C23  &
                      + BETA6*TERM6C23 + BETA9*TERM9C23)
      ENDIF

!.....Effective C_mu:
      CMUEFF(INP) = -0.5*(BETA1 + WII*BETA6) 
!.....Low Re correction for Effective Cmu:
      IF (LowRe) CMUEFF(INP) = CMUEFF(INP) * FACLR
!.....Use Limiter that Tom Gatski & Chris Rumsey use:
!      CMUEFF(INP) = max(CMUEFF(INP), 0.0005)

!----
      END DO !!I-loop
      END DO !!J-loop
      END DO !!K-loop
!=====END OF Velocity Gradient Tensor Invariants computation===========

!=====MAIN PART OF THE ROUTINE=======================================
      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM

      INP=LK(K)+LI(I)+J     

      WLDIST = wallDistance(inp)
!----------------------------------------------------------

!.....GRADIENT OF TURBULENCE KINETIC ENERGY
      DTEDX=gradTE(1,INP)
      DTEDY=gradTE(2,INP)
      DTEDZ=gradTE(3,INP)

!.....GRADIENT OF TURBULENCE KINETIC ENERGY SPECIFIC DISSIPATION RATE 
      DEDDX=gradED(1,INP)
      DEDDY=gradED(2,INP)
      DEDDZ=gradED(3,INP)

!=====FIND $D_{\omega}^{+}$ Dw+========================================
      TMP1=WLDIST**2/ED(INP) * (DTEDX*DEDDX+DTEDY*DEDDY+DTEDZ*DEDDZ)
      TMP2=200.*TEIN   ! 200*k_inf, koje moras definisati u MODINP.f i INC_INPUT
      DOMEGAPL=DMAX1(TMP1,TMP2)

!=====FIND KSI==========================================================
!      KSI=MIN(MAX(SQRT(TE(INP))/(0.09*WALL_DIST*ED(INP)),           &  !< SST radi uporedjivanja
!                 (500.*VISOB(INP)/DEN(INP))/(WALL_DIST**2*ED(INP))),&  !
!              4.*DEN(INP)*TE(INP)/(1.168*DOMEGAPL*WALL_DIST**2))       !
      IF (EARSM_WJ) THEN
      KSI=MIN(MAX(SQRT(TE(INP))/(0.09*WLDIST*ED(INP)),            &     !< EARSM 
                 (500.*VISOB(INP)/DEN(INP))/(WLDIST**2*ED(INP))), &     !
               20.*TE(INP)/DOMEGAPL)                                    !%$ WJ
      !ENDIF
      ELSE!IF (EARSM_M) THEN
      KSI=MIN(MAX(SQRT(TE(INP))/(0.09*WLDIST*ED(INP)),            &     !< EARSM 
                 (500.*VISOB(INP)/DEN(INP))/(WLDIST**2*ED(INP))), &     !
             2.*TE(INP)/DOMEGAPL)                                       ! Menter 2009
      ENDIF
!=====FIND F============================================================
      F=TANH(1.5*KSI**4)               ! A.Hellsten thesis - prolongs Omega region (!!!), usually factor is 1, here it's 1.5
      IF (EARSM_M) F=TANH(KSI**4)      ! Menter 2009

!=====NOW BLEND COEFFICIENTS WE'LL NEED LATER============================
!.....High-Re version.....................................................
      ALPHASST(INP) = F*0.518+(1.-F)*0.44      ! WJ and  A.Hellsten      !
      IF (EARSM_M) ALPHASST(INP) = F*0.553+(1.-F)*0.44   ! Menter 2009   !<
!.....Low-Re version of SST k-omega.......................................
!      alphast=(0.0249+(DENSIT*TE(IJK))/(6.*VISCOS*ED(IJK)))    &         !             
!             /(1.+(DENSIT*TE(IJK))/(6.*VISCOS*ED(IJK)))                  !
!      tmp=5./(9.*alphast)*                                     &         !                                     
!             (1./10.+ (DENSIT*TE(IJK))/(2.7*VISCOS*ED(IJK)))   &         !
!            /(1.   + (DENSIT*TE(IJK))/(2.7*VISCOS*ED(IJK)))              !
!      ALPHASST(INP) = F*tmp + (1.-F)*0.44                                !<
!........................................................................!
      BETTASST(INP) = F*0.0747+(1.-F)*0.0828  ! WJ and A.Hellsten
      IF (EARSM_M) BETTASST(INP) = F*0.075+(1.-F)*0.0828  ! Menter 2009


!.....EFFECTIVE DIFFUSIVITY (WJ and A.Hellsten original):
      IF (EARSM_WJ) THEN
      PRTINV_TE(INP)= F*1.1  + (1.-F)*1.1    ! sigma_k1=sigma_k2=1.1
      PRTINV_ED(INP)= F*0.53 + (1.-F)*1.0    ! sigma_omega1=0.53 sigma_omega2=1.0
      SIGMAD        = F*1.0  + (1.-F)*0.4    ! sigma_d
      !ENDIF

!.....EFFECTIVE DIFFUSIVITY (Menter 2009):
      ELSE!IF (EARSM_M) THEN
      PRTINV_TE(INP)= F*0.5 + (1.-F)*1.      ! sigma_k1=sigma_k2=1.1
      PRTINV_ED(INP)= F*0.5 + (1.-F)*0.856   ! sigma_omega1=0.53 sigma_omega2=1.0
      SIGMAD        =      2.*(1.-F)*0.856   ! sigma_d
      ENDIF

!=====FIND $D_{\omega}$ CROSS DIFFUSION MODIFICATION===================
      DOMEGA(INP)=SIGMAD * DEN(INP)/ED(INP) &
!                    *(DTEDX*DEDDX+DTEDY*DEDDY+DTEDZ*DEDDZ)
                  *max(DTEDX*DEDDX+DTEDY*DEDDY+DTEDZ*DEDDZ, 0.)

!=====ADDITIONAL SOURCE TERM FOR SAS MODEL: QSAS========================
      IF (SAS) THEN

!.....PROVIDE HIGH WAVE-NUMBER DUMPING BY LIMITING LVK
!     TODO CKECK IF NECESSARY PROVIDED THE FORM OF DISCRETIZATION WE USE 
!     FOR EVALUATING STRAIN.
!.....Do it like : Menter, Egorov, "Formulation of the Scale-Adaptive Simulation (SAS) Model during the DESIDER Project".
      DELTA=VOL(INP)**(1./3.)
      LVK(INP)=MAX(LVK(INP),0.26*DELTA)
!.....OR DO IT LIKE FLUENT DOES (COMMENT OUT IF UPPER APPROACH IS FAVOURED)
!#     LVK(INP)=MAX(LVK(INP), 0.26*SQRT(0.41*3.51/(BETTA/0.09-ALPHASST))*DELTA)  

!=====END von Karman LENGTH SCALE=======================================

!=====LENGTH SCALE OF MODELLED TURBULENCE===============================
!.....NOTE: C_MU**0.25=0.09**0.25=0.547722558
      LMT=DSQRT(DMAX1(TE(INP),ZERO))/(CMUEFF(INP)**0.25*ED(INP))

!=====FINALLY DEFINING QSAS TERM========================================
!.....NOTE: eta2*kappa=3.51*0.41=1.4391
      T1=1.4391*DEN(INP)*STRAIN(INP)**2*(LMT/LVK(INP))**2

      TMP1=(1./ED(INP)**2)*(DEDDX**2+DEDDY**2+DEDDZ**2)
      TMP2=(1./TE(INP)**2)*(DTEDX**2+DTEDY**2+DTEDZ**2)
!.....C*2/\sigma_{\Phi}=2*2/(2/3) = 6.
      T2=6.*DEN(INP)*TE(INP)*MAX(TMP1, TMP2)

      QSAS(INP)=MAX(T1-T2,0.)

!=====END QSAS TERM DEFINITION==========================================
      END IF

      END DO !!J-loop
      END DO !!I-loop
      END DO !!K-loop


      RETURN
      END
