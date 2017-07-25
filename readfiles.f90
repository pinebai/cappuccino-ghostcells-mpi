!***********************************************************************
!
      SUBROUTINE READFILES
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE GEOMETRY
      USE VARIABLES
      USE TITLE_MOD
      USE BUOY
      USE TIME_MOD
      USE STATISTICS

      IMPLICIT NONE  
!
!***********************************************************************
! 
      INTEGER :: INP

      OPEN(UNIT=3,FILE=restart_file,FORM='UNFORMATTED')
      REWIND 3

      READ(3) ITIME,TIME
      if(const_mflux) READ(3) gradPcmf
      READ(3)(F1(INP),INP=ICST,ICEN)
      READ(3)(F2(INP),INP=ICST,ICEN)
      READ(3)(F3(INP),INP=ICST,ICEN)
      READ(3)(U(INP),INP=ICST,ICEN)
      READ(3)(V(INP),INP=ICST,ICEN)
      READ(3)(W(INP),INP=ICST,ICEN)
      READ(3)(P(INP),INP=ICST,ICEN)
      READ(3)(TE(INP),INP=ICST,ICEN)
      READ(3)(ED(INP),INP=ICST,ICEN)
      READ(3)(T(INP),INP=ICST,ICEN)
      READ(3)(VIS(INP),INP=ICST,ICEN)
      !READ(3)(VISOB(INP),INP=ICST,ICEN)
      !READ(3)(VART(INP),INP=ICST,ICEN)
      !READ(3)(EDD(INP),INP=ICST,ICEN)
      !READ(3)(RET(INP),INP=ICST,ICEN)
      !READ(3)(DEN(INP),INP=ICST,ICEN)
      !READ(3)(UTT(INP),INP=ICST,ICEN)
      !READ(3)(VTT(INP),INP=ICST,ICEN)
      !READ(3)(WTT(INP),INP=ICST,ICEN)
      READ(3)(UU(INP),INP=ICST,ICEN)
      READ(3)(VV(INP),INP=ICST,ICEN)
      READ(3)(WW(INP),INP=ICST,ICEN)
      READ(3)(UV(INP),INP=ICST,ICEN)
      READ(3)(UW(INP),INP=ICST,ICEN)
      READ(3)(VW(INP),INP=ICST,ICEN)
      READ(3)(UO(INP),INP=ICST,ICEN)
      READ(3)(VO(INP),INP=ICST,ICEN)
      READ(3)(WO(INP),INP=ICST,ICEN)
      !READ(3)(TO(INP),INP=ICST,ICEN)
      READ(3)(TEO(INP),INP=ICST,ICEN)
      READ(3)(EDO(INP),INP=ICST,ICEN)
      !READ(3)(VARTO(INP),INP=ICST,ICEN)
      !READ(3)(CON(INP),INP=ICST,ICEN)
      !READ(3)(CONO(INP),INP=ICST,ICEN)
      !READ(3)(ALPH(INP),INP=ICST,ICEN)
      REWIND 3
      CLOSE (3)

      if (ltransient) then
!------------------------------------------------
!     [read statistics after first collection: ]
!------------------------------------------------
     OPEN(UNIT=85,FILE=trim(out_folder_path)//'/statistics1')   ! <- N_SAMPLE is here, statistics restart file 1
     OPEN(UNIT=86,FILE=trim(out_folder_path)//'/statistics2')   ! <- U_AVER, V_AVER,... are here, statistics restart file 2
     REWIND 85
     REWIND 86

     READ(85,*) N_SAMPLE
     READ(86,*) U_AVER,V_AVER,W_AVER, &
                UU_AVER,VV_AVER,WW_AVER, &
                UV_AVER,UW_AVER,VW_AVER, &
                TE_AVER
     CLOSE (85)
     CLOSE (86)
      endif

      RETURN
      END
