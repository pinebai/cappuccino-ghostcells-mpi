!***********************************************************************
!
      SUBROUTINE CALC_STATISTICS
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE GEOMETRY
      USE VARIABLES
      USE TIME_MOD
      USE STATISTICS

      IMPLICIT NONE
!
!***********************************************************************
!
      INTEGER :: I, J, K, INP
      REAL(PREC) :: U_NSAMPLE,V_NSAMPLE,W_NSAMPLE!,TE_NSAMPLE !,CON_NSAMPLE


      N_SAMPLE=N_SAMPLE+1
      
      ! LOOP ALL THE CELLS INCLUDING GHOST CELLS
      DO K=2,NKM
      DO I=2,NIM
      DO J=2,NJM

      INP=LK(K)+LI(I)+J

!.....Velocity field
      U_AVER(INP)=U_AVER(INP)+U(INP) 
      V_AVER(INP)=V_AVER(INP)+V(INP) 
      W_AVER(INP)=W_AVER(INP)+W(INP) 

!      CON_AVER(INP)=CON_AVER(INP)+CON(INP)

!.....Ensemble average over N samples
      U_NSAMPLE=U_AVER(INP)/N_SAMPLE
      V_NSAMPLE=V_AVER(INP)/N_SAMPLE
      W_NSAMPLE=W_AVER(INP)/N_SAMPLE
!      CON_NSAMPLE=CON_AVER(INP)/N_SAMPLE

!.....Reynolds stress components
      UU_AVER(INP)=UU_AVER(INP)+(U(INP)-U_NSAMPLE)**2
      VV_AVER(INP)=VV_AVER(INP)+(V(INP)-V_NSAMPLE)**2
      WW_AVER(INP)=WW_AVER(INP)+(W(INP)-W_NSAMPLE)**2

      UV_AVER(INP)=UV_AVER(INP)+ &
                   ((U(INP)-U_NSAMPLE)*(V(INP)-V_NSAMPLE))
      UW_AVER(INP)=UW_AVER(INP)+ &
                   ((U(INP)-U_NSAMPLE)*(W(INP)-W_NSAMPLE))
      VW_AVER(INP)=VW_AVER(INP)+ &
                   ((V(INP)-V_NSAMPLE)*(W(INP)-W_NSAMPLE))

!.....Turbulence kinetic energy
      TE_AVER(INP)=(UU_AVER(INP)+VV_AVER(INP)+WW_AVER(INP)) &
                        /(2*N_SAMPLE)

!.....Concentration
!      CONCON_AVER(INP)=CONCON_AVER(INP)+ &
!                   (CON(INP)-CON_NSAMPLE)**2 

!.....Concentration flux     
!      UCON_AVER(INP)=UCON_AVER(INP)+ &
!                   ((U(INP)-U_NSAMPLE)*(CON(INP)-CON_NSAMPLE))
!      VCON_AVER(INP)=VCON_AVER(INP)+ &
!                   ((V(INP)-V_NSAMPLE)*(CON(INP)-CON_NSAMPLE))
!      WCON_AVER(INP)=WCON_AVER(INP)+ &
!                   ((W(INP)-W_NSAMPLE)*(CON(INP)-CON_NSAMPLE))

      END DO    
      END DO    
      END DO    
 
      RETURN
      END
