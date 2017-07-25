!***********************************************************************
!
      SUBROUTINE Diagonally_Preconditioned_CG(FI,IFI)
!
!***********************************************************************
!
!    This routine incorporates the Diagonaly (Jacobi) Preconditioned 
!    Conjugate Gradient solver for symmetric matrices in 3D problems
!    with seven-diagonal matrix structure (see Sect. 5.3.6).
!
!    Array index IJK converted from 3D indices I, J, and K according to
!    Table 3.1. NS is the number of inner iterations (sweeps).
!
!    Original CG code written by Ismet Demirdzic, Sarajevo, 1991.
!    Modified by Nikola Mirkov, 28.01.2014. nmirkov@vinca.rs
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE COEF
      USE TITLE_MOD

      IMPLICIT NONE
!
!***********************************************************************
!
      INTEGER, INTENT(IN) :: IFI
      REAL(PREC), DIMENSION(NXYZA) :: FI 

!
!     LOCAL VARIABLES
!
      INTEGER :: I, J, K, IJK, NS, L
      REAL(PREC), DIMENSION(NXYZA) :: PK,ZK
      REAL(PREC) :: RSM, RESMAX, RES0, RESL
      REAL(PREC) :: S0, SK, ALF, BET, PKAPK

!.....MAX NO. OF INNER ITERS
      RESMAX = SOR(IFI)
!
!.....INITALIZE WORKING ARRAYS
!
      DO IJK=1,NIJK
        PK(IJK)=0.0D0
        ZK(IJK)=0.0D0
        RES(IJK)=0.0D0
      END DO
!
!.....CALCULATE INITIAL RESIDUAL VECTOR AND THE NORM
!
      RES0=0.0D0
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
             RES(IJK)=AE(IJK)*FI(IJK+NJ)+AW(IJK)*FI(IJK-NJ)+AN(IJK)* &
             FI(IJK+1)+AS(IJK)*FI(IJK-1)+AT(IJK)*FI(IJK+NIJ)+ &
             AB(IJK)*FI(IJK-NIJ)+SU(IJK)-AP(IJK)*FI(IJK)
            RES0=RES0+ABS(RES(IJK))
          END DO
        END DO
      END DO
!
!.....IF LTEST=True, PRINT THE NORM 
!
      IF(LTEST) WRITE(66,'(a,1PE10.3)') '                    RES0 = ',RES0
!
      S0=1.E20
!
!....START INNER ITERATIONS
!
      NS=NSW(IFI)
      DO L=1,NS
!
!.....DIAGONAL (JACOBI) PRECONDITIONING-SOLVE FOR ZK(IJK); CALCULATE SCALAR PRODUCT SK
!
      SK=0.0D0
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=RES(IJK)/(AP(IJK)+SMALL)
            SK=SK+RES(IJK)*ZK(IJK)
          END DO
        END DO
      END DO
!
!.....CALCULATE BETA
!
      BET=SK/S0
!
!.....CALCULATE NEW SEARCH VECTOR PK
!
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            PK(IJK)=ZK(IJK)+BET*PK(IJK)
          END DO
        END DO
      END DO
!
!.... CALCULATE SCALAR PRODUCT (PK.A PK) AND ALPHA (OVERWRITE ZK)
!
      PKAPK=0.0D0
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=AP(IJK)*PK(IJK)-AE(IJK)*PK(IJK+NJ)-                &
              AW(IJK)*PK(IJK-NJ)-AN(IJK)*PK(IJK+1)-AS(IJK)*PK(IJK-1)-  &
              AT(IJK)*PK(IJK+NIJ)-AB(IJK)*PK(IJK-NIJ)
            PKAPK=PKAPK+PK(IJK)*ZK(IJK)
          END DO
        END DO
      END DO

      ALF=SK/PKAPK
!
!.....CALCULATE VARIABLE CORRECTION, NEW RESIDUAL VECTOR, AND NORM
!
      RESL=0.0D0
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            FI(IJK)=FI(IJK)+ALF*PK(IJK)
            RES(IJK)=RES(IJK)-ALF*ZK(IJK)
            RESL=RESL+ABS(RES(IJK))
          END DO
        END DO
      END DO

      S0=SK
!
!.....CHECK CONVERGENCE
!
      IF(L.EQ.1) RESOR(IFI)=RES0
      RSM=RESL/(RESOR(IFI)+SMALL)
      IF(LTEST) WRITE(66,'(19x,3a,I4,a,1PE10.3,a,1PE10.3)') ' FI=',CHVAR(IFI),' SWEEP = ',L,' RESL = ',RESL,' RSM = ',RSM
      IF(RSM.LT.RESMAX) EXIT
!
!.....END OF ITERATION LOOP
!
      END DO

!.....Write linear solver report:
      write(66,'(3a,1PE10.3,a,1PE10.3,a,I0)') &
      'PCG(Jacobi):  Solving for ',trim(chvarSolver(IFI)), &
      ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L
!
      RETURN
      END
