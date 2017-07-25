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



      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM

      INP=LK(K)+LI(I)+J

      VOLR=1./VOL(INP)

      UVELE=U(INP+NJ)*FX(INP)+U(INP)*(1.-FX(INP))
      UVELW=U(INP)*FX(INP-NJ)+U(INP-NJ)*(1.-FX(INP-NJ))
      UVELN=U(INP+1)*FY(INP)+U(INP)*(1.-FY(INP))
      UVELS=U(INP)*FY(INP-1)+U(INP-1)*(1.-FY(INP-1))
      UVELT=U(INP+NIJ)*FZ(INP)+U(INP)*(1.-FZ(INP))
      UVELB=U(INP)*FZ(INP-NIJ)+U(INP-NIJ)*(1.-FZ(INP-NIJ))

      VVELE=V(INP+NJ)*FX(INP)+V(INP)*(1.-FX(INP))
      VVELW=V(INP)*FX(INP-NJ)+V(INP-NJ)*(1.-FX(INP-NJ))
      VVELN=V(INP+1)*FY(INP)+V(INP)*(1.-FY(INP))
      VVELS=V(INP)*FY(INP-1)+V(INP-1)*(1.-FY(INP-1))
      VVELT=V(INP+NIJ)*FZ(INP)+V(INP)*(1.-FZ(INP))
      VVELB=V(INP)*FZ(INP-NIJ)+V(INP-NIJ)*(1.-FZ(INP-NIJ))

      WVELE=W(INP+NJ)*FX(INP)+W(INP)*(1.-FX(INP))
      WVELW=W(INP)*FX(INP-NJ)+W(INP-NJ)*(1.-FX(INP-NJ))
      WVELN=W(INP+1)*FY(INP)+W(INP)*(1.-FY(INP))
      WVELS=W(INP)*FY(INP-1)+W(INP-1)*(1.-FY(INP-1))
      WVELT=W(INP+NIJ)*FZ(INP)+W(INP)*(1.-FZ(INP))
      WVELB=W(INP)*FZ(INP-NIJ)+W(INP-NIJ)*(1.-FZ(INP-NIJ))
!---------------------------------------------------------
!     [AE] ----->>  DUDX
!     [AN] ----->>  DUDY
!     [AT] ----->>  DUDZ
!     [AW] ----->>  DVDX
!     [AS] ----->>  DVDY
!     [AB] ----->>  DVDZ
!     [BE] ----->>  DWDX
!     [BN] ----->>  DWDY
!     [BT] ----->>  DWDZ
!---------------------------------------------------------
      AE(INP)=((UVELE-UVELW)*AR1X(INP)+ &
               (UVELN-UVELS)*AR2X(INP)+ &
               (UVELT-UVELB)*AR3X(INP))*VOLR

      AN(INP)=((UVELE-UVELW)*AR1Y(INP)+ &
               (UVELN-UVELS)*AR2Y(INP)+ &
               (UVELT-UVELB)*AR3Y(INP))*VOLR

      AT(INP)=((UVELE-UVELW)*AR1Z(INP)+ &
               (UVELN-UVELS)*AR2Z(INP)+ &
               (UVELT-UVELB)*AR3Z(INP))*VOLR
!----
      AW(INP)=((VVELE-VVELW)*AR1X(INP)+ &
               (VVELN-VVELS)*AR2X(INP)+ &
               (VVELT-VVELB)*AR3X(INP))*VOLR

      AS(INP)=((VVELE-VVELW)*AR1Y(INP)+ &
               (VVELN-VVELS)*AR2Y(INP)+ &
               (VVELT-VVELB)*AR3Y(INP))*VOLR

      AB(INP)=((VVELE-VVELW)*AR1Z(INP)+ &
               (VVELN-VVELS)*AR2Z(INP)+ &
               (VVELT-VVELB)*AR3Z(INP))*VOLR
!----
      BE(INP)=((WVELE-WVELW)*AR1X(INP)+ &
               (WVELN-WVELS)*AR2X(INP)+ &
               (WVELT-WVELB)*AR3X(INP))*VOLR

      BN(INP)=((WVELE-WVELW)*AR1Y(INP)+ &
               (WVELN-WVELS)*AR2Y(INP)+ &
               (WVELT-WVELB)*AR3Y(INP))*VOLR

      BT(INP)=((WVELE-WVELW)*AR1Z(INP)+ &
               (WVELN-WVELS)*AR2Z(INP)+ &
               (WVELT-WVELB)*AR3Z(INP))*VOLR
!----
      END DO !!I-loop
      END DO !!J-loop
      END DO !!K-loop
!
!-----------------------------------------------------
!     [U-VELOCITY COMPONENT: ]
!-----------------------------------------------------
      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM

      INP=LK(K)+LI(I)+J

      TERM1E=(DEN(INP+NJ)*UU(INP+NJ)+(AE(INP+NJ)+AE(INP+NJ))* &
             (VIS(INP+NJ)-VISOB(INP+NJ)))*FX(INP)+ &
             (DEN(INP)*UU(INP)+(AE(INP)+AE(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM1W=(DEN(INP)*UU(INP)+(AE(INP)+AE(INP))* &
             (VIS(INP)-VISOB(INP)))*FX(INP-NJ)+ &
             (DEN(INP-NJ)*UU(INP-NJ)+(AE(INP-NJ)+AE(INP-NJ))* &
             (VIS(INP-NJ)-VISOB(INP-NJ)))*(1.-FX(INP-NJ))

      TERM1N=(DEN(INP+1)*UU(INP+1)+(AE(INP+1)+AE(INP+1))* &
             (VIS(INP+1)-VISOB(INP+1)))*FY(INP)+ &
             (DEN(INP)*UU(INP)+(AE(INP)+AE(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM1S=(DEN(INP)*UU(INP)+(AE(INP)+AE(INP))* &
             (VIS(INP)-VISOB(INP)))*FY(INP-1)+ &
             (DEN(INP-1)*UU(INP-1)+(AE(INP-1)+AE(INP-1))* &
             (VIS(INP-1)-VISOB(INP-1)))*(1.-FY(INP-1))

      TERM1T=(DEN(INP+NIJ)*UU(INP+NIJ)+(AE(INP+NIJ)+AE(INP+NIJ))* &
             (VIS(INP+NIJ)-VISOB(INP+NIJ)))*FZ(INP)+ &
             (DEN(INP)*UU(INP)+(AE(INP)+AE(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM1B=(DEN(INP)*UU(INP)+(AE(INP)+AE(INP))* &
             (VIS(INP)-VISOB(INP)))*FZ(INP-NIJ)+ &
             (DEN(INP-NIJ)*UU(INP-NIJ)+(AE(INP-NIJ)+AE(INP-NJ))* &
             (VIS(INP-NIJ)-VISOB(INP-NIJ)))*(1.-FZ(INP-NIJ))

      DTERM1DX=((TERM1E-TERM1W)*AR1X(INP)+ &
                (TERM1N-TERM1S)*AR2X(INP)+ &
                (TERM1T-TERM1B)*AR3X(INP))
!-----------------------------------------------------------------
!
      TERM2E=(DEN(INP+NJ)*UV(INP+NJ)+(AN(INP+NJ)+AW(INP+NJ))* &
             (VIS(INP+NJ)-VISOB(INP+NJ)))*FX(INP)+ &
             (DEN(INP)*UV(INP)+(AN(INP)+AW(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM2W=(DEN(INP)*UV(INP)+(AN(INP)+AW(INP))* &
             (VIS(INP)-VISOB(INP)))*FX(INP-NJ)+ &
             (DEN(INP-NJ)*UV(INP-NJ)+(AN(INP-NJ)+AW(INP-NJ))* &
             (VIS(INP-NJ)-VISOB(INP-NJ)))*(1.-FX(INP-NJ))

      TERM2N=(DEN(INP+1)*UV(INP+1)+(AN(INP+1)+AW(INP+1))* &
             (VIS(INP+1)-VISOB(INP+1)))*FY(INP)+ &
             (DEN(INP)*UV(INP)+(AN(INP)+AW(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM2S=(DEN(INP)*UV(INP)+(AN(INP)+AW(INP))* &
             (VIS(INP)-VISOB(INP)))*FY(INP-1)+ &
             (DEN(INP-1)*UV(INP-1)+(AN(INP-1)+AW(INP-1))* &
             (VIS(INP-1)-VISOB(INP-1)))*(1.-FY(INP-1))

      TERM2T=(DEN(INP+NIJ)*UV(INP+NIJ)+(AN(INP+NIJ)+AW(INP+NIJ))* &
             (VIS(INP+NIJ)-VISOB(INP+NIJ)))*FZ(INP)+ &
             (DEN(INP)*UV(INP)+(AN(INP)+AW(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM2B=(DEN(INP)*UV(INP)+(AN(INP)+AW(INP))* &
             (VIS(INP)-VISOB(INP)))*FZ(INP-NIJ)+ &
             (DEN(INP-NIJ)*UV(INP-NIJ)+(AN(INP-NIJ)+AW(INP-NJ))* &
             (VIS(INP-NIJ)-VISOB(INP-NIJ)))*(1.-FZ(INP-NIJ))

      DTERM2DY=((TERM2E-TERM2W)*AR1Y(INP)+ &
                (TERM2N-TERM2S)*AR2Y(INP)+ &
                (TERM2T-TERM2B)*AR3Y(INP))
!------------------------------------------------------------------
!
      TERM3E=(DEN(INP+NJ)*UW(INP+NJ)+(AT(INP+NJ)+BE(INP+NJ))* &
             (VIS(INP+NJ)-VISOB(INP+NJ)))*FX(INP)+ &
             (DEN(INP)*UW(INP)+(AT(INP)+BE(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM3W=(DEN(INP)*UW(INP)+(AT(INP)+BE(INP))* &
             (VIS(INP)-VISOB(INP)))*FX(INP-NJ)+ &
             (DEN(INP-NJ)*UW(INP-NJ)+(AT(INP-NJ)+BE(INP-NJ))* &
             (VIS(INP-NJ)-VISOB(INP-NJ)))*(1.-FX(INP-NJ))

      TERM3N=(DEN(INP+1)*UW(INP+1)+(AT(INP+1)+BE(INP+1))* &
             (VIS(INP+1)-VISOB(INP+1)))*FY(INP)+ &
             (DEN(INP)*UW(INP)+(AT(INP)+BE(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM3S=(DEN(INP)*UW(INP)+(AT(INP)+BE(INP))* &
             (VIS(INP)-VISOB(INP)))*FY(INP-1)+ &
             (DEN(INP-1)*UW(INP-1)+(AT(INP-1)+BE(INP-1))* &
             (VIS(INP-1)-VISOB(INP-1)))*(1.-FY(INP-1))

      TERM3T=(DEN(INP+NIJ)*UW(INP+NIJ)+(AT(INP+NIJ)+BE(INP+NIJ))* &
             (VIS(INP+NIJ)-VISOB(INP+NIJ)))*FZ(INP)+ &
             (DEN(INP)*UW(INP)+(AT(INP)+BE(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM3B=(DEN(INP)*UW(INP)+(AT(INP)+BE(INP))* &
             (VIS(INP)-VISOB(INP)))*FZ(INP-NIJ)+ &
             (DEN(INP-NIJ)*UW(INP-NIJ)+(AT(INP-NIJ)+BE(INP-NJ))* &
             (VIS(INP-NIJ)-VISOB(INP-NIJ)))*(1.-FZ(INP-NIJ))

      DTERM3DZ=((TERM3E-TERM3W)*AR1Z(INP)+ &
                (TERM3N-TERM3S)*AR2Z(INP)+ &
                (TERM3T-TERM3B)*AR3Z(INP))
!
!------------------------------------------------------------------
      TEE=TE(INP+NJ)*FX(INP)+TE(INP)*(1.-FX(INP))
      TEW=TE(INP)*FX(INP-NJ)+TE(INP-NJ)*(1.-FX(INP-NJ))
      TEN=TE(INP+1)*FY(INP)+TE(INP)*(1.-FY(INP))
      TES=TE(INP)*FY(INP-1)+TE(INP-1)*(1.-FY(INP-1))
      TET=TE(INP+NIJ)*FZ(INP)+TE(INP)*(1.-FZ(INP))
      TEB=TE(INP)*FZ(INP-NIJ)+TE(INP-NIJ)*(1.-FZ(INP-NIJ))

      DTEDX=((TEE-TEW)*AR1X(INP)+(TEN-TES)*AR2X(INP)+ &
             (TET-TEB)*AR3X(INP))
!------------------------------------------------------------------

      SU(INP)=SU(INP)-DTERM1DX-DTERM2DY-DTERM3DZ+2*DTEDX/3.

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

      TERM1E=(DEN(INP+NJ)*UV(INP+NJ)+(AW(INP+NJ)+AN(INP+NJ))* &
             (VIS(INP+NJ)-VISOB(INP+NJ)))*FX(INP)+ &
             (DEN(INP)*UV(INP)+(AW(INP)+AN(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM1W=(DEN(INP)*UV(INP)+(AW(INP)+AN(INP))* &
             (VIS(INP)-VISOB(INP)))*FX(INP-NJ)+ &
             (DEN(INP-NJ)*UV(INP-NJ)+(AW(INP-NJ)+AN(INP-NJ))* &
             (VIS(INP-NJ)-VISOB(INP-NJ)))*(1.-FX(INP-NJ))

      TERM1N=(DEN(INP+1)*UV(INP+1)+(AW(INP+1)+AN(INP+1))* &
             (VIS(INP+1)-VISOB(INP+1)))*FY(INP)+ &
             (DEN(INP)*UV(INP)+(AW(INP)+AN(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM1S=(DEN(INP)*UV(INP)+(AW(INP)+AN(INP))* &
             (VIS(INP)-VISOB(INP)))*FY(INP-1)+ &
             (DEN(INP-1)*UV(INP-1)+(AW(INP-1)+AN(INP-1))* &
             (VIS(INP-1)-VISOB(INP-1)))*(1.-FY(INP-1))

      TERM1T=(DEN(INP+NIJ)*UV(INP+NIJ)+(AW(INP+NIJ)+AN(INP+NIJ))* &
             (VIS(INP+NIJ)-VISOB(INP+NIJ)))*FZ(INP)+ &
             (DEN(INP)*UV(INP)+(AW(INP)+AN(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM1B=(DEN(INP)*UV(INP)+(AW(INP)+AN(INP))* &
             (VIS(INP)-VISOB(INP)))*FZ(INP-NIJ)+ &
             (DEN(INP-NIJ)*UV(INP-NIJ)+(AW(INP-NIJ)+AN(INP-NJ))* &
             (VIS(INP-NIJ)-VISOB(INP-NIJ)))*(1.-FZ(INP-NIJ))

      DTERM1DX=(TERM1E-TERM1W)*AR1X(INP)+ &
               (TERM1N-TERM1S)*AR2X(INP)+ &
               (TERM1T-TERM1B)*AR3X(INP)
!-----------------------------------------------------------------
!
      TERM2E=(DEN(INP+NJ)*VV(INP+NJ)+(AS(INP+NJ)+AS(INP+NJ))* &
             (VIS(INP+NJ)-VISOB(INP+NJ)))*FX(INP)+ &
             (DEN(INP)*VV(INP)+(AS(INP)+AS(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM2W=(DEN(INP)*VV(INP)+(AS(INP)+AS(INP))* &
             (VIS(INP)-VISOB(INP)))*FX(INP-NJ)+ &
             (DEN(INP-NJ)*VV(INP-NJ)+(AS(INP-NJ)+AS(INP-NJ))* &
             (VIS(INP-NJ)-VISOB(INP-NJ)))*(1.-FX(INP-NJ))

      TERM2N=(DEN(INP+1)*VV(INP+1)+(AS(INP+1)+AS(INP+1))* &
             (VIS(INP+1)-VISOB(INP+1)))*FY(INP)+ &
             (DEN(INP)*VV(INP)+(AS(INP)+AS(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM2S=(DEN(INP)*VV(INP)+(AS(INP)+AS(INP))* &
             (VIS(INP)-VISOB(INP)))*FY(INP-1)+ &
             (DEN(INP-1)*VV(INP-1)+(AS(INP-1)+AS(INP-1))* &
             (VIS(INP-1)-VISOB(INP-1)))*(1.-FY(INP-1))

      TERM2T=(DEN(INP+NIJ)*VV(INP+NIJ)+(AS(INP+NIJ)+AS(INP+NIJ))* &
             (VIS(INP+NIJ)-VISOB(INP+NIJ)))*FZ(INP)+ &
             (DEN(INP)*VV(INP)+(AS(INP)+AS(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM2B=(DEN(INP)*VV(INP)+(AS(INP)+AS(INP))* &
             (VIS(INP)-VISOB(INP)))*FZ(INP-NIJ)+ &
             (DEN(INP-NIJ)*VV(INP-NIJ)+(AS(INP-NIJ)+AS(INP-NJ))* &
             (VIS(INP-NIJ)-VISOB(INP-NIJ)))*(1.-FZ(INP-NIJ))

      DTERM2DY=(TERM2E-TERM2W)*AR1Y(INP)+ &
               (TERM2N-TERM2S)*AR2Y(INP)+ &
               (TERM2T-TERM2B)*AR3Y(INP)

!------------------------------------------------------------------
!
      TERM3E=(DEN(INP+NJ)*VW(INP+NJ)+(AB(INP+NJ)+BN(INP+NJ))* &
             (VIS(INP+NJ)-VISOB(INP+NJ)))*FX(INP)+ &
             (DEN(INP)*VW(INP)+(AT(INP)+BE(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM3W=(DEN(INP)*VW(INP)+(AB(INP)+BN(INP))* &
             (VIS(INP)-VISOB(INP)))*FX(INP-NJ)+ &
             (DEN(INP-NJ)*VW(INP-NJ)+(AB(INP-NJ)+BN(INP-NJ))* &
             (VIS(INP-NJ)-VISOB(INP-NJ)))*(1.-FX(INP-NJ))

      TERM3N=(DEN(INP+1)*VW(INP+1)+(AB(INP+1)+BN(INP+1))* &
             (VIS(INP+1)-VISOB(INP+1)))*FY(INP)+ &
             (DEN(INP)*VW(INP)+(AB(INP)+BN(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM3S=(DEN(INP)*VW(INP)+(AB(INP)+BN(INP))* &
             (VIS(INP)-VISOB(INP)))*FY(INP-1)+ &
             (DEN(INP-1)*VW(INP-1)+(AB(INP-1)+BN(INP-1))* &
             (VIS(INP-1)-VISOB(INP-1)))*(1.-FY(INP-1))

      TERM3T=(DEN(INP+NIJ)*VW(INP+NIJ)+(AB(INP+NIJ)+BN(INP+NIJ))* &
             (VIS(INP+NIJ)-VISOB(INP+NIJ)))*FZ(INP)+ &
             (DEN(INP)*VW(INP)+(AB(INP)+BN(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM3B=(DEN(INP)*VW(INP)+(AB(INP)+BN(INP))* &
             (VIS(INP)-VISOB(INP)))*FZ(INP-NIJ)+ &
             (DEN(INP-NIJ)*VW(INP-NIJ)+(AB(INP-NIJ)+BN(INP-NJ))* &
             (VIS(INP-NIJ)-VISOB(INP-NIJ)))*(1.-FZ(INP-NIJ))

      DTERM3DZ=(TERM3E-TERM3W)*AR1Z(INP)+ &
               (TERM3N-TERM3S)*AR2Z(INP)+ &
               (TERM3T-TERM3B)*AR3Z(INP)
!
!------------------------------------------------------------------
      TEE=TE(INP+NJ)*FX(INP)+TE(INP)*(1.-FX(INP))
      TEW=TE(INP)*FX(INP-NJ)+TE(INP-NJ)*(1.-FX(INP-NJ))
      TEN=TE(INP+1)*FY(INP)+TE(INP)*(1.-FY(INP))
      TES=TE(INP)*FY(INP-1)+TE(INP-1)*(1.-FY(INP-1))
      TET=TE(INP+NIJ)*FZ(INP)+TE(INP)*(1.-FZ(INP))
      TEB=TE(INP)*FZ(INP-NIJ)+TE(INP-NIJ)*(1.-FZ(INP-NIJ))

      DTEDY=((TEE-TEW)*AR1Y(INP)+(TEN-TES)*AR2Y(INP)+ &
             (TET-TEB)*AR3Y(INP))
!------------------------------------------------------------------

      SV(INP)=SV(INP)-DTERM1DX-DTERM2DY-DTERM3DZ+2*DTEDY/3.

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

      TERM1E=(DEN(INP+NJ)*UW(INP+NJ)+(BE(INP+NJ)+AT(INP+NJ))* &
             (VIS(INP+NJ)-VISOB(INP+NJ)))*FX(INP)+ &
             (DEN(INP)*UW(INP)+(BE(INP)+AT(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM1W=(DEN(INP)*UW(INP)+(BE(INP)+AT(INP))* &
             (VIS(INP)-VISOB(INP)))*FX(INP-NJ)+ &
             (DEN(INP-NJ)*UW(INP-NJ)+(BE(INP-NJ)+AT(INP-NJ))* &
             (VIS(INP-NJ)-VISOB(INP-NJ)))*(1.-FX(INP-NJ))

      TERM1N=(DEN(INP+1)*UW(INP+1)+(BE(INP+1)+AT(INP+1))* &
             (VIS(INP+1)-VISOB(INP+1)))*FY(INP)+ &
             (DEN(INP)*UW(INP)+(BE(INP)+AT(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM1S=(DEN(INP)*UW(INP)+(BE(INP)+AT(INP))* &
             (VIS(INP)-VISOB(INP)))*FY(INP-1)+ &
             (DEN(INP-1)*UW(INP-1)+(BE(INP-1)+AT(INP-1))* &
             (VIS(INP-1)-VISOB(INP-1)))*(1.-FY(INP-1))

      TERM1T=(DEN(INP+NIJ)*UW(INP+NIJ)+(BE(INP+NIJ)+AT(INP+NIJ))* &
             (VIS(INP+NIJ)-VISOB(INP+NIJ)))*FZ(INP)+ &
             (DEN(INP)*UW(INP)+(BE(INP)+AT(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM1B=(DEN(INP)*UW(INP)+(BE(INP)+AT(INP))* &
             (VIS(INP)-VISOB(INP)))*FZ(INP-NIJ)+ &
             (DEN(INP-NIJ)*UW(INP-NIJ)+(BE(INP-NIJ)+AT(INP-NJ))* &
             (VIS(INP-NIJ)-VISOB(INP-NIJ)))*(1.-FZ(INP-NIJ))

      DTERM1DX=((TERM1E-TERM1W)*AR1X(INP)+ &
                (TERM1N-TERM1S)*AR2X(INP)+ &
                (TERM1T-TERM1B)*AR3X(INP))
!-----------------------------------------------------------------
!
      TERM2E=(DEN(INP+NJ)*VW(INP+NJ)+(BN(INP+NJ)+AB(INP+NJ))* &
             (VIS(INP+NJ)-VISOB(INP+NJ)))*FX(INP)+ &
             (DEN(INP)*VW(INP)+(BN(INP)+AB(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM2W=(DEN(INP)*VW(INP)+(BN(INP)+AB(INP))* &
             (VIS(INP)-VISOB(INP)))*FX(INP-NJ)+ &
             (DEN(INP-NJ)*VW(INP-NJ)+(BN(INP-NJ)+AB(INP-NJ))* &
             (VIS(INP-NJ)-VISOB(INP-NJ)))*(1.-FX(INP-NJ))

      TERM2N=(DEN(INP+1)*VW(INP+1)+(BN(INP+1)+AB(INP+1))* &
             (VIS(INP+1)-VISOB(INP+1)))*FY(INP)+ &
             (DEN(INP)*VW(INP)+(BN(INP)+AB(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM2S=(DEN(INP)*VW(INP)+(BN(INP)+AB(INP))* &
             (VIS(INP)-VISOB(INP)))*FY(INP-1)+ &
             (DEN(INP-1)*VW(INP-1)+(BN(INP-1)+AB(INP-1))* &
             (VIS(INP-1)-VISOB(INP-1)))*(1.-FY(INP-1))

      TERM2T=(DEN(INP+NIJ)*VW(INP+NIJ)+(BN(INP+NIJ)+AB(INP+NIJ))* &
             (VIS(INP+NIJ)-VISOB(INP+NIJ)))*FZ(INP)+ &
             (DEN(INP)*VW(INP)+(BN(INP)+AB(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM2B=(DEN(INP)*VW(INP)+(BN(INP)+AB(INP))* &
             (VIS(INP)-VISOB(INP)))*FZ(INP-NIJ)+ &
             (DEN(INP-NIJ)*VW(INP-NIJ)+(BN(INP-NIJ)+AB(INP-NJ))* &
             (VIS(INP-NIJ)-VISOB(INP-NIJ)))*(1.-FZ(INP-NIJ))


      DTERM2DY=(TERM2E-TERM2W)*AR1Y(INP)+ &
               (TERM2N-TERM2S)*AR2Y(INP)+ &
               (TERM2T-TERM2B)*AR3Y(INP)
!------------------------------------------------------------------
!
      TERM3E=(DEN(INP+NJ)*WW(INP+NJ)+(BT(INP+NJ)+BT(INP+NJ))* &
             (VIS(INP+NJ)-VISOB(INP+NJ)))*FX(INP)+ &
             (DEN(INP)*WW(INP)+(BT(INP)+BT(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FX(INP))

      TERM3W=(DEN(INP)*WW(INP)+(BT(INP)+BT(INP))* &
             (VIS(INP)-VISOB(INP)))*FX(INP-NJ)+ &
             (DEN(INP-NJ)*WW(INP-NJ)+(BT(INP-NJ)+BT(INP-NJ))* &
             (VIS(INP-NJ)-VISOB(INP-NJ)))*(1.-FX(INP-NJ))

      TERM3N=(DEN(INP+1)*WW(INP+1)+(BT(INP+1)+BT(INP+1))* &
             (VIS(INP+1)-VISOB(INP+1)))*FY(INP)+ &
             (DEN(INP)*WW(INP)+(BT(INP)+BT(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FY(INP))

      TERM3S=(DEN(INP)*WW(INP)+(BT(INP)+BT(INP))* &
             (VIS(INP)-VISOB(INP)))*FY(INP-1)+ &
             (DEN(INP-1)*WW(INP-1)+(BT(INP-1)+BT(INP-1))* &
             (VIS(INP-1)-VISOB(INP-1)))*(1.-FY(INP-1))

      TERM3T=(DEN(INP+NIJ)*WW(INP+NIJ)+(BT(INP+NIJ)+BT(INP+NIJ))* &
             (VIS(INP+NIJ)-VISOB(INP+NIJ)))*FZ(INP)+ &
             (DEN(INP)*WW(INP)+(BT(INP)+BT(INP))* &
             (VIS(INP)-VISOB(INP)))*(1.-FZ(INP))

      TERM3B=(DEN(INP)*WW(INP)+(BT(INP)+BT(INP))* &
             (VIS(INP)-VISOB(INP)))*FZ(INP-NIJ)+ &
             (DEN(INP-NIJ)*WW(INP-NIJ)+(BT(INP-NIJ)+BT(INP-NJ))* &
             (VIS(INP-NIJ)-VISOB(INP-NIJ)))*(1.-FZ(INP-NIJ))

      DTERM3DZ=(TERM3E-TERM3W)*AR1Z(INP)+ &
               (TERM3N-TERM3S)*AR2Z(INP)+ &
               (TERM3T-TERM3B)*AR3Z(INP)
!
!------------------------------------------------------------------
      TEE=TE(INP+NJ)*FX(INP)+TE(INP)*(1.-FX(INP))
      TEW=TE(INP)*FX(INP-NJ)+TE(INP-NJ)*(1.-FX(INP-NJ))
      TEN=TE(INP+1)*FY(INP)+TE(INP)*(1.-FY(INP))
      TES=TE(INP)*FY(INP-1)+TE(INP-1)*(1.-FY(INP-1))
      TET=TE(INP+NIJ)*FZ(INP)+TE(INP)*(1.-FZ(INP))
      TEB=TE(INP)*FZ(INP-NIJ)+TE(INP-NIJ)*(1.-FZ(INP-NIJ))

      DTEDZ=((TEE-TEW)*AR1Z(INP)+(TEN-TES)*AR2Z(INP)+ &
             (TET-TEB)*AR3Z(INP))
!------------------------------------------------------------------

      SW(INP)=SW(INP)-DTERM1DX-DTERM2DY-DTERM3DZ+2*DTEDZ/3.

      END DO !!I-loop
      END DO !!J-loop
      END DO !!K-loop
!
!----------------------------
!   [clean coefficient matrices]
!----------------------------
      DO INP=ICST,ICEN
      AE(INP)=0.0D0
      AN(INP)=0.0D0
      AT(INP)=0.0D0
      AW(INP)=0.0D0
      AS(INP)=0.0D0
      AB(INP)=0.0D0
      BE(INP)=0.0D0
      BN(INP)=0.0D0
      BT(INP)=0.0D0
      END DO
!-------------------------------------------------------------------

      RETURN
      END
