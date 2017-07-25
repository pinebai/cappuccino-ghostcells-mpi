!***********************************************************************
!
      SUBROUTINE calculate_stuff_for_sst
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
      USE GRADIENTS
      USE OMEGA_Turb_Models
      USE fieldManipulation

      IMPLICIT NONE
!
!***********************************************************************
!
      INTEGER :: I,J,K,INP

      REAL(PREC) :: WLDIST,DOMEGAPL,KSI,F,TMP1,TMP2          ! SST related
      REAL(PREC) :: DELTA,LMT,T1,T2,tmp                      ! SST and SAS Related
      REAL(PREC) :: DTEDX,DTEDY,DTEDZ,DEDDX,DEDDY,DEDDZ
      REAL(PREC) :: alphast

!.....If needed get the Von Karman length scale for SAS model.
      if(SAS) lvk(:) = von_karman_lengthscale()


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
      TMP1=2.*DEN(INP)/(1.168*ED(INP))  &
             *(DTEDX*DEDDX+DTEDY*DEDDY+DTEDZ*DEDDZ)
      TMP2=1.E-10
      DOMEGAPL=MAX(TMP1,TMP2)

!=====FIND KSI==========================================================
      KSI=MIN(MAX(SQRT(TE(INP))/(0.09*WLDIST*ED(INP)),  &
                 (500.*VISCOS/DEN(INP))/(WLDIST**2*ED(INP))),  &
              4.*DEN(INP)*TE(INP)/(1.168*DOMEGAPL*WLDIST**2))

!=====FIND F============================================================
      F=TANH(KSI**4)


!=====NOW BLEND COEFFICIENTS WE'LL NEED LATER============================
!.....High-Re version....................................................
      ALPHASST(INP)=F*ALPHA1+(1.-F)*ALPHA2                               !<

      If(LowRe) then
!.....Low-Re version of SST k-omega......................................
      alphast=(0.024+(DENSIT*TE(INP))/(6.*VISCOS*ED(INP)))  &           !             
            /(1.+(DENSIT*TE(INP))/(6.*VISCOS*ED(INP)))                  !
      tmp=ALPHA1/alphast*                                   &           !                                     
             (1./9.+ (DENSIT*TE(INP))/(2.95*VISCOS*ED(INP)))&           !
            /(1.   + (DENSIT*TE(INP))/(2.95*VISCOS*ED(INP)))            !
      ALPHASST(INP)=F*tmp + (1.-F)*ALPHA2                               !<
!.......................................................................!
      endif

      BETTASST(INP)=F*BETAI1+(1.-F)*BETAI2

!=====Potrebno za Hybrid Seamless Alpha model:==========================
!     Ovde zelimo da razlika izmedju koeficijenata bude funkcija polozaja bas
!     kao sto su i sami koeficijenti.
      IF (ALPHAMODEL.and.SST) DIFF(INP) = BETTASST(INP)-ALPHASST(INP)

!=====EFFECTIVE DIFFUSIVITY:============================================
      PRTINV_TE(INP)=F*SIGMK1  + (1.-F)*SIGMK2
      PRTINV_ED(INP)=F*SIGMOM1 + (1.-F)*SIGMOM2

!=====FIND D_omega CROSS DIFFUSION MODIFICATION:========================
      DOMEGA(INP)=2.*(1.-F)*DEN(INP)/(1.168*ED(INP)) &
                  *(DTEDX*DEDDX+DTEDY*DEDDY+DTEDZ*DEDDZ)
      DOMEGA(INP) = max(DOMEGA(INP),0.)

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
      LMT=DSQRT(MAX(TE(INP),ZERO))/(CMU25*ED(INP))

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
