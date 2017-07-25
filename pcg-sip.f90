!***********************************************************************
!
      SUBROUTINE PCGSIP(FI,IFI)
!
!***********************************************************************
!
!    This routine incorporates the SIP Preconditioned 
!    Conjugate Gradient solver for symmetric matrices in 3D problems
!    with seven-diagonal matrix structure.
!
!    Written by Nikola Mirkov, 28.01.2014. nmirkov@vinca.rs
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE COEF
      USE COEFB
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
      REAL(PREC) :: RSM, RESMAX, RES0, RESL, P1, P2, P3
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

!.....CALCULATE COEFFICIENTS OF  L  AND  U  MATRICES USING SIP
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
          IJK=LK(K)+LI(I)+J
          BB(IJK)=-AB(IJK)/(1.+ALFA*(BN(IJK-NIJ)+BE(IJK-NIJ)))
          BW(IJK)=-AW(IJK)/(1.+ALFA*(BN(IJK-NJ)+BT(IJK-NJ)))
          BS(IJK)=-AS(IJK)/(1.+ALFA*(BE(IJK-1)+BT(IJK-1)))
          P1=ALFA*(BB(IJK)*BN(IJK-NIJ)+BW(IJK)*BN(IJK-NJ))
          P2=ALFA*(BB(IJK)*BE(IJK-NIJ)+BS(IJK)*BE(IJK-1))
          P3=ALFA*(BW(IJK)*BT(IJK-NJ)+BS(IJK)*BT(IJK-1))
          BP(IJK)=1./(AP(IJK)+P1+P2+P3 &
            -BB(IJK)*BT(IJK-NIJ) &
            -BW(IJK)*BE(IJK-NJ) &
            -BS(IJK)*BN(IJK-1)+SMALL)
           BN(IJK)=(-AN(IJK)-P1)*BP(IJK)
           BE(IJK)=(-AE(IJK)-P2)*BP(IJK)
           BT(IJK)=(-AT(IJK)-P3)*BP(IJK)
          END DO
        END DO
      END DO
!
      S0=1.E20
!
!....START INNER ITERATIONS
!
      NS=NSW(IFI)
      DO L=1,NS
!
!.....SOLVE FOR ZK(IJK) -- FORWARD ELIMINATION
!
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=(RES(IJK)-BB(IJK)*ZK(IJK-NIJ)-BW(IJK)*ZK(IJK-NJ)- &
            BS(IJK)*ZK(IJK-1))*BP(IJK)
          END DO
        END DO
      END DO

!
!..... BACKWARD SUBSTITUTION; CALCULATE SCALAR PRODUCT SK
!
      SK=0.0D0
      DO K=NKMM,3,-1
        DO I=NIMM,3,-1
          DO J=NJMM,3,-1
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=ZK(IJK)-BN(IJK)*ZK(IJK+1)-BE(IJK)*ZK(IJK+NJ)- &
                    BT(IJK)*ZK(IJK+NIJ)
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
      'PCG(SIP):  Solving for ',trim(chvarSolver(IFI)), &
      ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L
!
      RETURN
      END
