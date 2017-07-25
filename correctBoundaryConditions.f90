      SUBROUTINE correctBoundaryConditions
!  
!******************************************************************************
!
!     Updates values at boundaries with fixedGradient and zeroGradient B.C. and
!     periodicity after solution for $\phi$ is obtained.
! 
!******************************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE COEF
      USE GEOMETRY
      USE VARIABLES
      USE BC
      USE BOUNDC
      USE WALL

      IMPLICIT NONE
!
!     LOCAL VARIABLES
!
      INTEGER :: I, J, K, INP, INH

!.....UPDATE FAR EAST
      DO K=2,NKM
      DO J=2,NJM
      INP=LK(K)+LI(NI)+J !gde staviti
      INH=LK(K)+LI(2)+J  !odakle
        U(INP)=U(INH)
        V(INP)=V(INH)
        W(INP)=W(INH) 
        !TE(INP)=TE(INH)
        !ED(INP)=ED(INH)
        !VIS(INP)=VIS(INH)
        !F1(INP)=F1(INH)
      ENDDO
      ENDDO
!.....UPDATE FAR WEST
      DO K=2,NKM
      DO J=2,NJM
      INP=LK(K)+LI(1)+J    !gde staviti
      INH=LK(K)+LI(NJM)+J  !odakle
        U(INP)=U(INH)
        V(INP)=V(INH)
        W(INP)=W(INH) 
        !TE(INP)=TE(INH)
        !ED(INP)=ED(INH)
        !VIS(INP)=VIS(INH)
        F1(INP)=F1(INH)
      ENDDO
      ENDDO
!.....UPDATE NORTH
      DO K=2,NKM
      DO I=2,NIM
      INP=LK(K)+LI(I)+NJ
      !INH=INP-1 !<--ukoliko je symmetry
      INH=LK(K)+LI(I)+2 !<--periodicity
        U(INP)=U(INH)
        V(INP)=V(INH)
        W(INP)=W(INH) 
        !TE(INP)=TE(INH)
        !ED(INP)=ED(INH)
        !VIS(INP)=VIS(INH)
        !F2(INP)=F2(INH)
      ENDDO
      ENDDO
!.....UPDATE SOUTH
      DO K=2,NKM
      DO I=2,NIM
      INP=LK(K)+LI(I)+1
      !INH=INP+1 !<--ukoliko je symmetry
      INH=LK(K)+LI(I)+NJM
        U(INP)=U(INH)
        V(INP)=V(INH)
        W(INP)=W(INH) 
        !TE(INP)=TE(INH)
        !ED(INP)=ED(INH)
        !VIS(INP)=VIS(INH)
        F2(INP)=F2(INH)
      ENDDO
      ENDDO

      RETURN
      END
