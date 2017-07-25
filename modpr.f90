!***********************************************************************
!
      SUBROUTINE MODPR
!
!***********************************************************************
!
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
      USE VEC
      USE OBSTACLE
  
      IMPLICIT NONE 
!
!***********************************************************************
!
      INTEGER :: NSA, K, I, J, IJK 
      INTEGER :: IOSS, IOES, JOSS, JOES, KOSS, KOES ! GRANICE PARAMETARA ZA DO PETLJU
!                                                                       
!
!.....PRESSURE BOUND. BOUNDARY CONDITIONS
!
      IF(IOBST.EQ.0) RETURN
      DO NSA=1,NOBST
      IOSS=IOS(NSA)
      IOES=IOE(NSA)
      JOSS=JOS(NSA)
      JOES=JOE(NSA)
      KOSS=KOS(NSA)
      KOES=KOE(NSA)

      DO K=KOSS,KOES
      DO I=IOSS,IOES
      DO J=JOSS,JOES
      IJK=LK(K)+LI(I)+J
      SU(IJK)=0.0d0
      SP(IJK)=GREAT
      END DO
      END DO
      END DO

      END DO  ! NSA=1,NOBST

      RETURN
      END
