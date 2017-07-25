!***********************************************************************
!
      SUBROUTINE Additional_algebraic_stress_terms
!
!***********************************************************************
!     Calculates the additional agebraic stress terms for momentum eq.
!     and ads them to SU, SV, SW rhs. vectors.
!
!***********************************************************************
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE GEOMETRY
      USE COEF
      USE COEFB
      USE VARIABLES
      USE OBSTACLE
      USE GRADIENTS

      IMPLICIT NONE
!
!***********************************************************************
!
!
!     LOCAL VARIABLES
!
      INTEGER :: I, J, K, INP
      REAL(PREC) :: VOLR, &
                    TEE, TEW, TES, TEN, TET, TEB, &                               !#
                    DTEDX, DTEDY, DTEDZ, &                                        !
                    UVELE, UVELW, UVELN, UVELS, UVELT, UVELB, &                   ! VARS. FOR AFM
                    VVELE, VVELW, VVELN, VVELS, VVELT, VVELB, &                   !
                    WVELE, WVELW, WVELN, WVELS, WVELT, WVELB, &                   !#
                    TERM1E, TERM1W, TERM1N, TERM1S, TERM1T, TERM1B, DTERM1DX, &   !
                    TERM2E, TERM2W, TERM2N, TERM2S, TERM2T, TERM2B, DTERM2DY, &   !
                    TERM3E, TERM3W, TERM3N, TERM3S, TERM3T, TERM3B, DTERM3DZ      ! END VARS. FOR AFM




!     THE ADDITIONAL ALGEBRAIC STRESS TERMS
!
      two_thirds=2./3.
      INE = INP+NJ
      INW = INP-NJ
      INN = INP+1
      INS = INP-1
      INT = INP+NIJ
      INB = INP-NIJ
!-----------------------------------------------------
!     [U-VELOCITY COMPONENT: ]
!-----------------------------------------------------
      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM

      INP=LK(K)+LI(I)+J

      TERM1E=(DEN(INE)*UU(INE)+(gradU(1,INE)+gradU(1,INE))*(VIS(INE)-VISOB(INE)))*FX(INP)+ &
             (DEN(INP)*UU(INP)+(gradU(1,INP)+gradU(1,INP))*(VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM1W=(DEN(INP)*UU(INP)+(gradU(1,INP)+gradU(1,INP))*(VIS(INP)-VISOB(INP)))*FX(INW)+ &
             (DEN(INW)*UU(INW)+(gradU(1,INW)+gradU(1,INW))*(VIS(INW)-VISOB(INW)))*(1.-FX(INW))

      TERM1N=(DEN(INN)*UU(INN)+(gradU(1,INN)+gradU(1,INN))*(VIS(INN)-VISOB(INN)))*FY(INP)+ &
             (DEN(INP)*UU(INP)+(gradU(1,INP)+gradU(1,INP))*(VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM1S=(DEN(INP)*UU(INP)+(gradU(1,INP)+gradU(1,INP))*(VIS(INP)-VISOB(INP)))*FY(INS)+ &
             (DEN(INS)*UU(INS)+(gradU(1,INS)+gradU(1,INS))*(VIS(INS)-VISOB(INS)))*(1.-FY(INS))

      TERM1T=(DEN(INT)*UU(INT)+(gradU(1,INT)+gradU(1,INT))*(VIS(INT)-VISOB(INT)))*FZ(INP)+ &
             (DEN(INP)*UU(INP)+(gradU(1,INP)+gradU(1,INP))*(VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM1B=(DEN(INP)*UU(INP)+(gradU(1,INP)+gradU(1,INP))*(VIS(INP)-VISOB(INP)))*FZ(INB)+ &
             (DEN(INB)*UU(INB)+(gradU(1,INB)+gradU(1,INW))*(VIS(INB)-VISOB(INB)))*(1.-FZ(INB))

      DTERM1DX=((TERM1E-TERM1W)*AR1X(INP)+ &
                (TERM1N-TERM1S)*AR2X(INP)+ &
                (TERM1T-TERM1B)*AR3X(INP))
!-----------------------------------------------------------------
!
      TERM2E=(DEN(INE)*UV(INE)+(gradU(2,INE)+gradV(1,INE))*(VIS(INE)-VISOB(INE)))*FX(INP)+ &
             (DEN(INP)*UV(INP)+(gradU(2,INP)+gradV(1,INP))*(VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM2W=(DEN(INP)*UV(INP)+(gradU(2,INP)+gradV(1,INP))*(VIS(INP)-VISOB(INP)))*FX(INW)+ &
             (DEN(INW)*UV(INW)+(gradU(2,INW)+gradV(1,INW))*(VIS(INW)-VISOB(INW)))*(1.-FX(INW))

      TERM2N=(DEN(INN)*UV(INN)+(gradU(2,INN)+gradV(1,INN))*(VIS(INN)-VISOB(INN)))*FY(INP)+ &
             (DEN(INP)*UV(INP)+(gradU(2,INP)+gradV(1,INP))*(VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM2S=(DEN(INP)*UV(INP)+(gradU(2,INP)+gradV(1,INP))*(VIS(INP)-VISOB(INP)))*FY(INS)+ &
             (DEN(INS)*UV(INS)+(gradU(2,INS)+gradV(1,INS))*(VIS(INS)-VISOB(INS)))*(1.-FY(INS))

      TERM2T=(DEN(INT)*UV(INT)+(gradU(2,INT)+gradV(1,INT))*(VIS(INT)-VISOB(INT)))*FZ(INP)+ &
             (DEN(INP)*UV(INP)+(gradU(2,INP)+gradV(1,INP))*(VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM2B=(DEN(INP)*UV(INP)+(gradU(2,INP)+gradV(1,INP))*(VIS(INP)-VISOB(INP)))*FZ(INB)+ &
             (DEN(INB)*UV(INB)+(gradU(2,INB)+gradV(1,INW))*(VIS(INB)-VISOB(INB)))*(1.-FZ(INB))

      DTERM2DY=((TERM2E-TERM2W)*AR1Y(INP)+ &
                (TERM2N-TERM2S)*AR2Y(INP)+ &
                (TERM2T-TERM2B)*AR3Y(INP))
!------------------------------------------------------------------
!
      TERM3E=(DEN(INE)*UW(INE)+(gradU(3,INE)+gradW(1,INE))*(VIS(INE)-VISOB(INE)))*FX(INP)+ &
             (DEN(INP)*UW(INP)+(gradU(3,INP)+gradW(1,INP))*(VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM3W=(DEN(INP)*UW(INP)+(gradU(3,INP)+gradW(1,INP))*(VIS(INP)-VISOB(INP)))*FX(INW)+ &
             (DEN(INW)*UW(INW)+(gradU(3,INW)+gradW(1,INW))*(VIS(INW)-VISOB(INW)))*(1.-FX(INW))

      TERM3N=(DEN(INN)*UW(INN)+(gradU(3,INN)+gradW(1,INN))*(VIS(INN)-VISOB(INN)))*FY(INP)+ &
             (DEN(INP)*UW(INP)+(gradU(3,INP)+gradW(1,INP))*(VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM3S=(DEN(INP)*UW(INP)+(gradU(3,INP)+gradW(1,INP))*(VIS(INP)-VISOB(INP)))*FY(INS)+ &
             (DEN(INS)*UW(INS)+(gradU(3,INS)+gradW(1,INS))*(VIS(INS)-VISOB(INS)))*(1.-FY(INS))

      TERM3T=(DEN(INT)*UW(INT)+(gradU(3,INT)+gradW(1,INT))*(VIS(INT)-VISOB(INT)))*FZ(INP)+ &
             (DEN(INP)*UW(INP)+(gradU(3,INP)+gradW(1,INP))*(VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM3B=(DEN(INP)*UW(INP)+(gradU(3,INP)+gradW(1,INP))*(VIS(INP)-VISOB(INP)))*FZ(INB)+ &
             (DEN(INB)*UW(INB)+(gradU(3,INB)+gradW(1,INW))*(VIS(INB)-VISOB(INB)))*(1.-FZ(INB))

      DTERM3DZ=((TERM3E-TERM3W)*AR1Z(INP)+ &
                (TERM3N-TERM3S)*AR2Z(INP)+ &
                (TERM3T-TERM3B)*AR3Z(INP))
!
!------------------------------------------------------------------

      SU(INP)=SU(INP)-DTERM1DX-DTERM2DY-DTERM3DZ+two_thirds*gradTE(1,inp)

      END DO !!I-loop
      END DO !!J-loop
      END DO !!K-loop

!
!-----------------------------------------------------
!     [V-VELOCITY COMPONENT: ]
!-----------------------------------------------------
      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM

      INP=LK(K)+LI(I)+J

      TERM1E=(DEN(INE)*UV(INE)+(gradV(1,INE)+gradU(2,INE))*(VIS(INE)-VISOB(INE)))*FX(INP)+ &
             (DEN(INP)*UV(INP)+(gradV(1,INP)+gradU(2,INP))*(VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM1W=(DEN(INP)*UV(INP)+(gradV(1,INP)+gradU(2,INP))*(VIS(INP)-VISOB(INP)))*FX(INW)+ &
             (DEN(INW)*UV(INW)+(gradV(1,INW)+gradU(2,INW))*(VIS(INW)-VISOB(INW)))*(1.-FX(INW))

      TERM1N=(DEN(INN)*UV(INN)+(gradV(1,INN)+gradU(2,INN))*(VIS(INN)-VISOB(INN)))*FY(INP)+ &
             (DEN(INP)*UV(INP)+(gradV(1,INP)+gradU(2,INP))*(VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM1S=(DEN(INP)*UV(INP)+(gradV(1,INP)+gradU(2,INP))*(VIS(INP)-VISOB(INP)))*FY(INS)+ &
             (DEN(INS)*UV(INS)+(gradV(1,INS)+gradU(2,INS))*(VIS(INS)-VISOB(INS)))*(1.-FY(INS))

      TERM1T=(DEN(INT)*UV(INT)+(gradV(1,INT)+gradU(2,INT))*(VIS(INT)-VISOB(INT)))*FZ(INP)+ &
             (DEN(INP)*UV(INP)+(gradV(1,INP)+gradU(2,INP))*(VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM1B=(DEN(INP)*UV(INP)+(gradV(1,INP)+gradU(2,INP))*(VIS(INP)-VISOB(INP)))*FZ(INB)+ &
             (DEN(INB)*UV(INB)+(gradV(1,INB)+gradU(2,INW))*(VIS(INB)-VISOB(INB)))*(1.-FZ(INB))

      DTERM1DX=(TERM1E-TERM1W)*AR1X(INP)+ &
               (TERM1N-TERM1S)*AR2X(INP)+ &
               (TERM1T-TERM1B)*AR3X(INP)
!-----------------------------------------------------------------
!
      TERM2E=(DEN(INE)*VV(INE)+(gradV(2,INE)+gradV(2,INE))*(VIS(INE)-VISOB(INE)))*FX(INP)+ &
             (DEN(INP)*VV(INP)+(gradV(2,INP)+gradV(2,INP))*(VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM2W=(DEN(INP)*VV(INP)+(gradV(2,INP)+gradV(2,INP))*(VIS(INP)-VISOB(INP)))*FX(INW)+ &
             (DEN(INW)*VV(INW)+(gradV(2,INW)+gradV(2,INW))*(VIS(INW)-VISOB(INW)))*(1.-FX(INW))

      TERM2N=(DEN(INN)*VV(INN)+(gradV(2,INN)+gradV(2,INN))*(VIS(INN)-VISOB(INN)))*FY(INP)+ &
             (DEN(INP)*VV(INP)+(gradV(2,INP)+gradV(2,INP))*(VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM2S=(DEN(INP)*VV(INP)+(gradV(2,INP)+gradV(2,INP))*(VIS(INP)-VISOB(INP)))*FY(INS)+ &
             (DEN(INS)*VV(INS)+(gradV(2,INS)+gradV(2,INS))*(VIS(INS)-VISOB(INS)))*(1.-FY(INS))

      TERM2T=(DEN(INT)*VV(INT)+(gradV(2,INT)+gradV(2,INT))*(VIS(INT)-VISOB(INT)))*FZ(INP)+ &
             (DEN(INP)*VV(INP)+(gradV(2,INP)+gradV(2,INP))*(VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM2B=(DEN(INP)*VV(INP)+(gradV(2,INP)+gradV(2,INP))*(VIS(INP)-VISOB(INP)))*FZ(INB)+ &
             (DEN(INB)*VV(INB)+(gradV(2,INB)+gradV(2,INW))*(VIS(INB)-VISOB(INB)))*(1.-FZ(INB))

      DTERM2DY=(TERM2E-TERM2W)*AR1Y(INP)+ &
               (TERM2N-TERM2S)*AR2Y(INP)+ &
               (TERM2T-TERM2B)*AR3Y(INP)

!------------------------------------------------------------------
!
      TERM3E=(DEN(INE)*VW(INE)+(gradV(3,INE)+gradW(2,INE))*(VIS(INE)-VISOB(INE)))*FX(INP)+ &
             (DEN(INP)*VW(INP)+(gradU(3,INP)+gradW(1,INP))*(VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM3W=(DEN(INP)*VW(INP)+(gradV(3,INP)+gradW(2,INP))*(VIS(INP)-VISOB(INP)))*FX(INW)+ &
             (DEN(INW)*VW(INW)+(gradV(3,INW)+gradW(2,INW))*(VIS(INW)-VISOB(INW)))*(1.-FX(INW))

      TERM3N=(DEN(INN)*VW(INN)+(gradV(3,INN)+gradW(2,INN))*(VIS(INN)-VISOB(INN)))*FY(INP)+ &
             (DEN(INP)*VW(INP)+(gradV(3,INP)+gradW(2,INP))*(VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM3S=(DEN(INP)*VW(INP)+(gradV(3,INP)+gradW(2,INP))*(VIS(INP)-VISOB(INP)))*FY(INS)+ &
             (DEN(INS)*VW(INS)+(gradV(3,INS)+gradW(2,INS))*(VIS(INS)-VISOB(INS)))*(1.-FY(INS))

      TERM3T=(DEN(INT)*VW(INT)+(gradV(3,INT)+gradW(2,INT))*(VIS(INT)-VISOB(INT)))*FZ(INP)+ &
             (DEN(INP)*VW(INP)+(gradV(3,INP)+gradW(2,INP))*(VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM3B=(DEN(INP)*VW(INP)+(gradV(3,INP)+gradW(2,INP))*(VIS(INP)-VISOB(INP)))*FZ(INB)+ &
             (DEN(INB)*VW(INB)+(gradV(3,INB)+gradW(2,INW))*(VIS(INB)-VISOB(INB)))*(1.-FZ(INB))

      DTERM3DZ=(TERM3E-TERM3W)*AR1Z(INP)+ &
               (TERM3N-TERM3S)*AR2Z(INP)+ &
               (TERM3T-TERM3B)*AR3Z(INP)
!
!------------------------------------------------------------------

      SV(INP)=SV(INP)-DTERM1DX-DTERM2DY-DTERM3DZ+two_thirds*gradTE(2,inp)

      END DO !!I-loop
      END DO !!J-loop
      END DO !!K-loop
!
!-----------------------------------------------------
!     [W-VELOCITY COMPONENT: ]
!-----------------------------------------------------
      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM

      INP=LK(K)+LI(I)+J

      TERM1E=(DEN(INE)*UW(INE)+(gradW(1,INE)+gradU(3,INE))*(VIS(INE)-VISOB(INE)))*FX(INP)+ &
             (DEN(INP)*UW(INP)+(gradW(1,INP)+gradU(3,INP))*(VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM1W=(DEN(INP)*UW(INP)+(gradW(1,INP)+gradU(3,INP))*(VIS(INP)-VISOB(INP)))*FX(INW)+ &
             (DEN(INW)*UW(INW)+(gradW(1,INW)+gradU(3,INW))*VIS(INW)-VISOB(INW)))*(1.-FX(INW))

      TERM1N=(DEN(INN)*UW(INN)+(gradW(1,INN)+gradU(3,INN))*(VIS(INN)-VISOB(INN)))*FY(INP)+ &
             (DEN(INP)*UW(INP)+(gradW(1,INP)+gradU(3,INP))*(VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM1S=(DEN(INP)*UW(INP)+(gradW(1,INP)+gradU(3,INP))*(VIS(INP)-VISOB(INP)))*FY(INS)+ &
             (DEN(INS)*UW(INS)+(gradW(1,INS)+gradU(3,INS))*(VIS(INS)-VISOB(INS)))*(1.-FY(INS))

      TERM1T=(DEN(INT)*UW(INT)+(gradW(1,INT)+gradU(3,INT))*(VIS(INT)-VISOB(INT)))*FZ(INP)+ &
             (DEN(INP)*UW(INP)+(gradW(1,INP)+gradU(3,INP))*(VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM1B=(DEN(INP)*UW(INP)+(gradW(1,INP)+gradU(3,INP))*(VIS(INP)-VISOB(INP)))*FZ(INB)+ &
             (DEN(INB)*UW(INB)+(gradW(1,INB)+gradU(3,INW))*(VIS(INB)-VISOB(INB)))*(1.-FZ(INB))

      DTERM1DX=((TERM1E-TERM1W)*AR1X(INP)+ &
                (TERM1N-TERM1S)*AR2X(INP)+ &
                (TERM1T-TERM1B)*AR3X(INP))
!-----------------------------------------------------------------
!
      TERM2E=(DEN(INE)*VW(INE)+(gradW(2,INE)+gradV(3,INE))*(VIS(INE)-VISOB(INE)))*FX(INP)+ &
             (DEN(INP)*VW(INP)+(gradW(2,INP)+gradV(3,INP))*(VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM2W=(DEN(INP)*VW(INP)+(gradW(2,INP)+gradV(3,INP))*(VIS(INP)-VISOB(INP)))*FX(INW)+ &
             (DEN(INW)*VW(INW)+(gradW(2,INW)+gradV(3,INW))*(VIS(INW)-VISOB(INW)))*(1.-FX(INW))

      TERM2N=(DEN(INN)*VW(INN)+(gradW(2,INN)+gradV(3,INN))*(VIS(INN)-VISOB(INN)))*FY(INP)+ &
             (DEN(INP)*VW(INP)+(gradW(2,INP)+gradV(3,INP))*(VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM2S=(DEN(INP)*VW(INP)+(gradW(2,INP)+gradV(3,INP))*(VIS(INP)-VISOB(INP)))*FY(INS)+ &
             (DEN(INS)*VW(INS)+(gradW(2,INS)+gradV(3,INS))*(VIS(INS)-VISOB(INS)))*(1.-FY(INS))

      TERM2T=(DEN(INT)*VW(INT)+(gradW(2,INT)+gradV(3,INT))*(VIS(INT)-VISOB(INT)))*FZ(INP)+ &
             (DEN(INP)*VW(INP)+(gradW(2,INP)+gradV(3,INP))*(VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM2B=(DEN(INP)*VW(INP)+(gradW(2,INP)+gradV(3,INP))*(VIS(INP)-VISOB(INP)))*FZ(INB)+ &
             (DEN(INB)*VW(INB)+(gradW(2,INB)+gradV(3,INW))*(VIS(INB)-VISOB(INB)))*(1.-FZ(INB))


      DTERM2DY=(TERM2E-TERM2W)*AR1Y(INP)+ &
               (TERM2N-TERM2S)*AR2Y(INP)+ &
               (TERM2T-TERM2B)*AR3Y(INP)
!------------------------------------------------------------------
!
      TERM3E=(DEN(INE)*WW(INE)+(gradW(3,INE)+gradW(3,INE))*(VIS(INE)-VISOB(INE)))*FX(INP)+ &
             (DEN(INP)*WW(INP)+(gradW(3,INP)+gradW(3,INP))*(VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM3W=(DEN(INP)*WW(INP)+(gradW(3,INP)+gradW(3,INP))*(VIS(INP)-VISOB(INP)))*FX(INW)+ &
             (DEN(INW)*WW(INW)+(gradW(3,INW)+gradW(3,INW))*(VIS(INW)-VISOB(INW)))*(1.-FX(INW))

      TERM3N=(DEN(INN)*WW(INN)+(gradW(3,INN)+gradW(3,INN))*(VIS(INN)-VISOB(INN)))*FY(INP)+ &
             (DEN(INP)*WW(INP)+(gradW(3,INP)+gradW(3,INP))*(VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM3S=(DEN(INP)*WW(INP)+(gradW(3,INP)+gradW(3,INP))*(VIS(INP)-VISOB(INP)))*FY(INS)+ &
             (DEN(INS)*WW(INS)+(gradW(3,INS)+gradW(3,INS))*(VIS(INS)-VISOB(INS)))*(1.-FY(INS))

      TERM3T=(DEN(INT)*WW(INT)+(gradW(3,INT)+gradW(3,INT))*(VIS(INT)-VISOB(INT)))*FZ(INP)+ &
             (DEN(INP)*WW(INP)+(gradW(3,INP)+gradW(3,INP))*(VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM3B=(DEN(INP)*WW(INP)+(gradW(3,INP)+gradW(3,INP))*(VIS(INP)-VISOB(INP)))*FZ(INB)+ &
             (DEN(INB)*WW(INB)+(gradW(3,INB)+gradW(3,INW))*(VIS(INB)-VISOB(INB)))*(1.-FZ(INB))

      DTERM3DZ=(TERM3E-TERM3W)*AR1Z(INP)+ &
               (TERM3N-TERM3S)*AR2Z(INP)+ &
               (TERM3T-TERM3B)*AR3Z(INP)
!
!------------------------------------------------------------------

      SW(INP)=SW(INP)-DTERM1DX-DTERM2DY-DTERM3DZ+two_thirds*gradTE(3,inp)

      END DO !!J-loop
      END DO !!I-loop
      END DO !!K-loop


      RETURN
      END
