!***********************************************************************
!
       SUBROUTINE WRITEHISTORY
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE GEOMETRY
      USE VARIABLES
      USE TITLE_MOD
      USE TIME_MOD
      USE OMEGA_Turb_Models

      IMPLICIT NONE
!
!***********************************************************************
!
     INTEGER :: I, J, K, IJK

!---------------------------------------------------------
!       RESULTS AT EACH TIME-STEP FOR TRANSIENT SIMUL.
!---------------------------------------------------------

      IF(LTRANSIENT) THEN
      DO IMON=1,MPOINTS
        READ(89,*) I,J,K
        IJK=LK(K)+LI(I)+J
        WRITE(91+IMON,'(2X,1P7E14.5,2X)') TIME,U(IJK),V(IJK),W(IJK),TE(IJK),ED(IJK)
      END DO
      REWIND 89
      END IF

      RETURN
      END
