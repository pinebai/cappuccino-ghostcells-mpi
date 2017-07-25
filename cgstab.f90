      SUBROUTINE CGSTAB(FI,IFI)
!
!***********************************************************************
!
!    This routine incorporates the CGSTAB solver for seven-diagonal,
!    non-symmetric coefficient matrices (suitable for convection/
!    diffusion problems). See Sect. 5.3.7 for details. Array index
!    IJK calculated from indices I, J, and K according to Table 3.1.
!
!    Writen by Samir Muzaferija, Institut fuer Schiffbau, Hamburg, 1995.
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
      REAL(PREC), DIMENSION(NXYZA) :: RESO,PK,UK,ZK,VK,D
      REAL(PREC) :: RSM, RESMAX, RES0, RESL
      REAL(PREC) :: ALF, BETO, GAM, BET, OM, VRES, VV, UKRESO

!.....MAX NO. OF INNER ITERS
      RESMAX = SOR(IFI) 
!
!.....CALCULATE INITIAL RESIDUAL VECTOR
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
      IF(LTEST) WRITE(66,'(a,1PE10.3)') '                    RES0 = ',RES0
!
!.....CALCULATE ELEMENTS OF PRECONDITIONING MATRIX DIAGONAL
!
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            D(IJK)=1./(AP(IJK) - AW(IJK)*D(IJK-NJ)*AE(IJK-NJ) &
                   - AS(IJK)*D(IJK-1)*AN(IJK-1)               &
                   - AB(IJK)*D(IJK-NIJ)*AT(IJK-NIJ)) 
          END DO
        END DO
      END DO
!
!.....INITIALIZE WORKING ARRAYS AND CONSTANTS
!
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            RESO(IJK)=RES(IJK)
            PK(IJK)=0.0D0
            UK(IJK)=0.0D0
            ZK(IJK)=0.0D0
            VK(IJK)=0.0D0
          END DO
        END DO
      END DO
      ALF=1.0D0
      BETO=1.0D0
      GAM=1.0D0
!
!.....START INNER ITERATIONS
!
      NS=NSW(IFI)
      DO L=1,NS
!
!..... CALCULATE BETA AND OMEGA
!
      BET=0.0D0
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            BET=BET+RES(IJK)*RESO(IJK)
          END DO
        END DO
      END DO
      OM=BET*GAM/(ALF*BETO+SMALL)
      BETO=BET
!
!..... CALCULATE PK
!
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            PK(IJK)=RES(IJK)+OM*(PK(IJK)-ALF*UK(IJK))
          END DO
        END DO
      END DO
!
!.....SOLVE (M ZK = PK) - FORWARD SUBSTITUTION
!
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=(PK(IJK)+AW(IJK)*ZK(IJK-NJ)  &
                    +AS(IJK)*ZK(IJK-1)+AB(IJK)*ZK(IJK-NIJ))*D(IJK)
          END DO
        END DO
      END DO

      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=ZK(IJK)/(D(IJK)+SMALL)
          END DO
        END DO
      END DO
!
!..... BACKWARD SUBSTITUTION
!
      DO K=NKMM,3,-1
        DO I=NIMM,3,-1
          DO J=NJMM,3,-1
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=(ZK(IJK)+AE(IJK)*ZK(IJK+NJ)  &
                    +AN(IJK)*ZK(IJK+1)+AT(IJK)*ZK(IJK+NIJ))*D(IJK)
          END DO
        END DO
      END DO
!
!.....CALCULATE UK = A.PK
!
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            UK(IJK)=AP(IJK)*ZK(IJK)-AE(IJK)*ZK(IJK+NJ) -    &
                     AW(IJK)*ZK(IJK-NJ)-AN(IJK)*ZK(IJK+1)-  &
                     AS(IJK)*ZK(IJK-1)-AT(IJK)*ZK(IJK+NIJ)- &
                     AB(IJK)*ZK(IJK-NIJ)
          END DO
        END DO
      END DO
!
!..... CALCULATE SCALAR PRODUCT UK.RESO AND GAMMA
!
      UKRESO=0.0D0
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            UKRESO=UKRESO+UK(IJK)*RESO(IJK)
          END DO
        END DO
      END DO
      GAM=BET/UKRESO
!
!.....UPDATE (FI) AND CALCULATE W (OVERWRITE RES - IT IS RES-UPDATE)
!
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            FI(IJK)=FI(IJK)+GAM*ZK(IJK)   
            RES(IJK)=RES(IJK)-GAM*UK(IJK) !W
          END DO
        END DO
      END DO
!
!.....SOLVE (M Y = W); Y OVERWRITES ZK; FORWARD SUBSTITUTION
!
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=(RES(IJK)+AW(IJK)*ZK(IJK-NJ)+  &
                    AS(IJK)*ZK(IJK-1)+AB(IJK)*ZK(IJK-NIJ))*D(IJK)
           END DO
         END DO
      END DO

      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=ZK(IJK)/(D(IJK)+SMALL)
          END DO
        END DO
      END DO
!
!.....BACKWARD SUBSTITUTION
!
      DO K=NKMM,3,-1
        DO I=NIMM,3,-1
          DO J=NJMM,3,-1
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=(ZK(IJK)+AE(IJK)*ZK(IJK+NJ)+  &
                    AN(IJK)*ZK(IJK+1)+AT(IJK)*ZK(IJK+NIJ))*D(IJK)
          END DO
        END DO
      END DO
!
!.....CALCULATE V = A.Y (VK = A.ZK)
!
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            VK(IJK)=AP(IJK)*ZK(IJK)   -AE(IJK)*ZK(IJK+NJ)-  &
                    AW(IJK)*ZK(IJK-NJ)-AN(IJK)*ZK(IJK+1)-   &
                    AS(IJK)*ZK(IJK-1) -AT(IJK)*ZK(IJK+NIJ)- &
                    AB(IJK)*ZK(IJK-NIJ)
          END DO
        END DO
      END DO
!
!..... CALCULATE ALPHA (ALF)
!
      VRES=0.0D0
      VV=0.0D0
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            VRES=VRES+VK(IJK)*RES(IJK)
            VV=VV+VK(IJK)*VK(IJK)
          END DO
        END DO
      END DO

      ALF=VRES/(VV+SMALL)
!
!.....UPDATE VARIABLE (FI) AND RESIDUAL (RES) VECTORS
!
      RESL=0.0D0
      DO K=3,NKMM
        DO I=3,NIMM
          DO J=3,NJMM 
            IJK=LK(K)+LI(I)+J
            FI(IJK)=FI(IJK)+ALF*ZK(IJK)
            RES(IJK)=RES(IJK)-ALF*VK(IJK)
            RESL=RESL+ABS(RES(IJK))
          END DO
        END DO
      END DO
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
      'BiCGStab(ILU(0)):  Solving for ',trim(chvarSolver(IFI)), &
      ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L
!
      RETURN
      END

