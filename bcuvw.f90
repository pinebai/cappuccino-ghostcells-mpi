!***********************************************************************
!
      SUBROUTINE BCUVW
!
!***********************************************************************
!     UVW-MOMENTUM BOUNDARY CONDITIONS
!***********************************************************************
! 

      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE COEF
      USE GEOMETRY
      USE VARIABLES
      USE BC
      USE BOUNDC
      USE TITLE_MOD
      USE BUOY
      USE TIME_MOD
      USE VEC
      USE WALL
      USE OBSTACLE

      IMPLICIT NONE
!
!***********************************************************************
! 
      INTEGER :: I, J, K, INP, INBC, IK, JK, IJ, IJK, &
                 NSA, IOSS, IOES, JOSS, JOES, KOSS, KOES, &
                 I1, K1
!-----------------
!.....WEST
!-----------------
      DO 300 K=2,NKM
      DO 300 J=2,NJM
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
!
      INP=LK(K)+LI(2)+J
      INBC=(J-1)*NK+K
      LBHLP=LBW(INBC+JKS)
      IF(LBHLP.EQ.1) THEN
!.....INLET
if (periodic_boundary) cycle !!!<-za periodic !
      CALL INLBC(INP,NJ,1,NIJ,FY,FZ,F1(INP-NJ))
      ELSEIF(LBHLP.EQ.2) THEN
!.....OUTLET, CROSS-DERIVATES
      CALL OUTCS(INP,NJ,1,NIJ,FY,FZ)
      ELSEIF(LBHLP.EQ.4.) THEN
!.....WALL
      CALL WALLBC(INP,NJ,1,NIJ)
      GENTW(INBC)=GENT
      SUEDW(INBC)=SUED
      ELSEIF(LBHLP.EQ.3) THEN
!.....SYMMETRY
      CALL SYMBC(INP,NJ,1,NIJ)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SVU
      SW(INP)=SW(INP)+SWU
      AP(INP)=AP(INP)+SUP
      SPV(INP)=SPV(INP)+SVP
      SP(INP)=SP(INP)+SWP
      AW(INP)=0.
  300 CONTINUE
!-----------------
!.....EAST
!-----------------
      DO 305 K=2,NKM
      DO 305 J=2,NJM
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
      INP=LK(K)+LI(NIM)+J
      INBC=(J-1)*NK+K
      LBHLP=LBE(INBC+JKS)
      IF(LBHLP.EQ.1) THEN
!.....INLET
      CALL INLBC(INP,-NJ,1,NIJ,FY,FZ,-F1(INP))
      ELSEIF(LBHLP.EQ.2) THEN
!.....OUTLET, CROSS-DERIVATES
if (periodic_boundary) cycle  !!! <-za periodic - do nothing just cycle
      CALL OUTCS(INP,-NJ,1,NIJ,FY,FZ) 
      ELSEIF(LBHLP.EQ.4) THEN
!.....WALL
      CALL WALLBC(INP,-NJ,1,NIJ)
      GENTE(INBC)=GENT
      SUEDE(INBC)=SUED
      ELSEIF(LBHLP.EQ.3) THEN
!.....SYMMETRY
      CALL SYMBC(INP,-NJ,1,NIJ)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SVU
      SW(INP)=SW(INP)+SWU
      AP(INP)=AP(INP)+SUP
      SPV(INP)=SPV(INP)+SVP
      SP(INP)=SP(INP)+SWP
      AE(INP)=0.
  305 CONTINUE
!-----------------
!.....SOUTH
!-----------------
      DO 310 K=2,NKM
      DO 310 I=2,NIM
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
!
      INP=LK(K)+LI(I)+2
      INBC=(I-1)*NK+K
      LBHLP=LBS(INBC+IKS)
      IF(LBHLP.EQ.1) THEN
!.....INLET
      CALL INLBC(INP,1,NIJ,NJ,FZ,FX,F2(INP-1))
      ELSEIF(LBHLP.EQ.2) THEN
!.....OUTLET, CROSS-DERIVATES
      CALL OUTCS(INP,1,NIJ,NJ,FZ,FX)
      ELSEIF(LBHLP.EQ.4) THEN
!.....WALL
      CALL WALLBC(INP,1,NIJ,NJ)
      GENTS(INBC)=GENT
      SUEDS(INBC)=SUED
      ELSEIF(LBHLP.EQ.3) THEN
!.....SYMMETRY
if (periodic_boundary) cycle  !!! <-za periodic
      CALL SYMBC(INP,1,NIJ,NJ)
      ENDIF
      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SVU
      SW(INP)=SW(INP)+SWU
      AP(INP)=AP(INP)+SUP
      SPV(INP)=SPV(INP)+SVP
      SP(INP)=SP(INP)+SWP
      AS(INP)=0.
  310 CONTINUE
!-----------------
!.....NORTH
!-----------------
      DO 315 K=2,NKM
      DO 315 I=2,NIM
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
      INP=LK(K)+LI(I)+NJM
      INBC=(I-1)*NK+K
      LBHLP=LBN(INBC+IKS)
      IF(LBHLP.EQ.1) THEN
!.....INLET
      CALL INLBC(INP,-1,NIJ,NJ,FZ,FX,-F2(INP))
      ELSEIF(LBHLP.EQ.2) THEN
!.....OUTLET, CROSS-DERIVATES
      CALL OUTCS(INP,-1,NIJ,NJ,FZ,FX)
      ELSEIF(LBHLP.EQ.4) THEN
!.....WALL
      CALL WALLBC(INP,-1,NIJ,NJ)
      GENTN(INBC)=GENT
      SUEDN(INBC)=SUED
      ELSEIF(LBHLP.EQ.3) THEN
!.....SYMMETRY
if (periodic_boundary) cycle  !!!<-za periodic
      CALL SYMBC(INP,-1,NIJ,NJ)
      ENDIF
      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SVU
      SW(INP)=SW(INP)+SWU
      AP(INP)=AP(INP)+SUP
      SPV(INP)=SPV(INP)+SVP
      SP(INP)=SP(INP)+SWP
      AN(INP)=0.
  315 CONTINUE
!-----------------
!.....BOTTOM
!-----------------
      DO 320 I=2,NIM
      DO 320 J=2,NJM
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
!
      INP=LK(2)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBB(INBC+IJS)
      IF(LBHLP.EQ.1) THEN
!.....INLET
      CALL INLBC(INP,NIJ,NJ,1,FX,FY,F3(INP-NIJ))
      ELSEIF(LBHLP.EQ.2) THEN
!.....OUTLET, CROSS-DERIVATES
      CALL OUTCS(INP,NIJ,NJ,1,FX,FY)
      ELSEIF(LBHLP.EQ.4) THEN
!.....WALL
      CALL WALLBC(INP,NIJ,NJ,1)
      GENTB(INBC)=GENT
      SUEDB(INBC)=SUED
      ELSEIF(LBHLP.EQ.3) THEN
!.....SYMMETRY
      CALL SYMBC(INP,NIJ,NJ,1)
      ENDIF
      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SVU
      SW(INP)=SW(INP)+SWU
      AP(INP)=AP(INP)+SUP
      SPV(INP)=SPV(INP)+SVP
      SP(INP)=SP(INP)+SWP
      AB(INP)=0.
  320 CONTINUE
!-----------------
!.....TOP
!-----------------
      DO 325 I=2,NIM
      DO 325 J=2,NJM
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
      INP=LK(NKM)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBT(INBC+IJS)
      IF(LBHLP.EQ.1) THEN
!.....INLET
      CALL INLBC(INP,-NIJ,NJ,1,FX,FY,-F3(INP))
      ELSEIF(LBHLP.EQ.2) THEN
!.....OUTLET, CROSS-DERIVATES
      CALL OUTCS(INP,-NIJ,NJ,1,FX,FY)
      ELSEIF(LBHLP.EQ.4) THEN
!.....WALL
!--------------------------------------
!.....[LID DRIVEN CAVITY: ]
!--------------------------------------
!      U(INP+NIJ)=1.
!--------------------------------------
      CALL WALLBC(INP,-NIJ,NJ,1)
      GENTT(INBC)=GENT
      SUEDT(INBC)=SUED
      ELSEIF(LBHLP.EQ.3) THEN
!.....SYMMETRY
      CALL SYMBC(INP,-NIJ,NJ,1)
      ENDIF
!--------------------------------------
      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SVU
      SW(INP)=SW(INP)+SWU
      AP(INP)=AP(INP)+SUP
      SPV(INP)=SPV(INP)+SVP
      SP(INP)=SP(INP)+SWP
      AT(INP)=0.
  325 CONTINUE
!
!##########################################
!     [OBSTACLE MODIFICATIONS: ]
!##########################################
!
      IF(IOBST.EQ.0) RETURN
      DO NSA=1,NOBST
      IOSS=IOS(NSA)
      IOES=IOE(NSA)
      JOSS=JOS(NSA)
      JOES=JOE(NSA)
      KOSS=KOS(NSA)
      KOES=KOE(NSA)
!##########################################
!     [SOUTH WALL OF SUBMERGED BODY: ]
!##########################################
      DO I=IOSS,IOES
      I1=(I-IOSS)*(KOES-KOSS+1)
      DO K=KOSS,KOES
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
      IK=I1+K-KOSS+1
      IF(JOSS.NE.2.AND.TYPOBS(NSA,IK).EQ.1) THEN
      IJK=LI(I)+LK(K)+JOSS-1
      DX1=X(IJK)-X(IJK-NIJ-NJ)
      DY1=Y(IJK)-Y(IJK-NIJ-NJ)
      DZ1=Z(IJK)-Z(IJK-NIJ-NJ)
      DX2=X(IJK-NJ)-X(IJK-NIJ)
      DY2=Y(IJK-NJ)-Y(IJK-NIJ)
      DZ2=Z(IJK-NJ)-Z(IJK-NIJ)
      LW=IJK
      DELN=DNOS(NSA,IK)
      CALL WALLBCOBST
      SU(IJK)=SU(IJK)+SUU
      SV(IJK)=SV(IJK)+SVU
      SW(IJK)=SW(IJK)+SWU
      AP(IJK)=AP(IJK)+SUP
      SPV(IJK)=SPV(IJK)+SVP
      SP(IJK)=SP(IJK)+SWP
      GENOS(NSA,IK)=GENT
      SUEDOS(NSA,IK)=SUED
      AN(IJK)=0.
      END IF
!##########################################
!     [NORTH WALL OF SUBMERGED BODY: ]
!##########################################
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
      IF(JOES.NE.NJM.AND.TYPOBN(NSA,IK).EQ.1) THEN
      IJK=LI(I)+LK(K)+JOES+1
      DX1=X(IJK-NJ-1)-X(IJK-NIJ-1)
      DY1=Y(IJK-NJ-1)-Y(IJK-NIJ-1)
      DZ1=Z(IJK-NJ-1)-Z(IJK-NIJ-1)
      DX2=X(IJK-1)-X(IJK-NIJ-NJ-1)
      DY2=Y(IJK-1)-Y(IJK-NIJ-NJ-1)
      DZ2=Z(IJK-1)-Z(IJK-NIJ-NJ-1)
      LW=IJK
      DELN=DNON(NSA,IK)
      CALL WALLBCOBST
      SU(IJK)=SU(IJK)+SUU
      SV(IJK)=SV(IJK)+SVU
      SW(IJK)=SW(IJK)+SWU
      AP(IJK)=AP(IJK)+SUP
      SPV(IJK)=SPV(IJK)+SVP
      SP(IJK)=SP(IJK)+SWP
      GENON(NSA,IK)=GENT
      SUEDON(NSA,IK)=SUED
      AS(IJK)=0.
      END IF
      END DO
      END DO
!##########################################
!     [WEST WALL OF SUBMERGED BODY: ]
!##########################################
      DO K=KOSS,KOES
      K1=(K-KOSS)*(JOES-JOSS+1)
      DO J=JOSS,JOES
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
      JK=K1+J-JOSS+1
      IF(IOSS.NE.2.AND.TYPOBW(NSA,JK).EQ.1) THEN
      IJK=LK(K)+LI(IOSS-1)+J
      DX1=X(IJK-1)-X(IJK-NIJ)
      DY1=Y(IJK-1)-Y(IJK-NIJ)
      DZ1=Z(IJK-1)-Z(IJK-NIJ)
      DX2=X(IJK)-X(IJK-NIJ-1)
      DY2=Y(IJK)-Y(IJK-NIJ-1)
      DZ2=Z(IJK)-Z(IJK-NIJ-1)
      LW=IJK
      DELN=DNOW(NSA,JK)
      CALL WALLBCOBST
      SU(IJK)=SU(IJK)+SUU
      SV(IJK)=SV(IJK)+SVU
      SW(IJK)=SW(IJK)+SWU
      AP(IJK)=AP(IJK)+SUP
      SPV(IJK)=SPV(IJK)+SVP
      SP(IJK)=SP(IJK)+SWP
      GENOW(NSA,JK)=GENT
      SUEDOW(NSA,JK)=SUED
      AE(IJK)=0.
      END IF
!##########################################
!     [EAST WALL OF SUBMERGED BODY: ]
!##########################################
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
      IF(IOES.NE.NIM.AND.TYPOBE(NSA,JK).EQ.1) THEN
      IJK=LK(K)+LI(IOES+1)+J
      DX1=X(IJK-NJ)-X(IJK-NIJ-NJ-1)
      DY1=Y(IJK-NJ)-Y(IJK-NIJ-NJ-1)
      DZ1=Z(IJK-NJ)-Z(IJK-NIJ-NJ-1)
      DX2=X(IJK-NJ-1)-X(IJK-NIJ-NJ)
      DY2=Y(IJK-NJ-1)-Y(IJK-NIJ-NJ)
      DZ2=Z(IJK-NJ-1)-Z(IJK-NIJ-NJ)
      LW=IJK
      DELN=DNOE(NSA,JK)
      CALL WALLBCOBST
      SU(IJK)=SU(IJK)+SUU
      SV(IJK)=SV(IJK)+SVU
      SW(IJK)=SW(IJK)+SWU
      AP(IJK)=AP(IJK)+SUP
      SPV(IJK)=SPV(IJK)+SVP
      SP(IJK)=SP(IJK)+SWP
      GENOE(NSA,JK)=GENT
      SUEDOE(NSA,JK)=SUED
      AW(IJK)=0.
      END IF
      END DO
      END DO
!##########################################
!     [TOP WALL OF SUBMERGED BODY: ]
!##########################################
      DO I=IOSS,IOES
      I1=(I-IOSS)*(JOES-JOSS+1)
      DO J=JOSS,JOES
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
      IJ=I1+J-JOSS+1
      IF(KOES.NE.NKM.AND.TYPOBT(NSA,IJ).EQ.1) THEN
      IJK=LI(I)+LK(KOES+1)+J
      DX1=X(IJK-NIJ-1)-X(IJK-NIJ-NJ)
      DY1=Y(IJK-NIJ-1)-Y(IJK-NIJ-NJ)
      DZ1=Z(IJK-NIJ-1)-Z(IJK-NIJ-NJ)
      DX2=X(IJK-NIJ)-X(IJK-NIJ-NJ-1)
      DY2=Y(IJK-NIJ)-Y(IJK-NIJ-NJ-1)
      DZ2=Z(IJK-NIJ)-Z(IJK-NIJ-NJ-1)
      LW=IJK
      DELN=DNOT(NSA,IJ)
      CALL WALLBCOBST
      SU(IJK)=SU(IJK)+SUU
      SV(IJK)=SV(IJK)+SVU
      SW(IJK)=SW(IJK)+SWU
      AP(IJK)=AP(IJK)+SUP
      SPV(IJK)=SPV(IJK)+SVP
      SP(IJK)=SP(IJK)+SWP
      GENOT(NSA,IJ)=GENT
      SUEDOT(NSA,IJ)=SUED
      AB(IJK)=0.
      END IF
!##########################################
!     [BOTTOM WALL OF SUBMERGED BODY: ]
!##########################################
      SUU=0.
      SVU=0.
      SWU=0.
      SUP=0.
      SVP=0.
      SWP=0.
      IF(KOSS.NE.2.AND.TYPOBB(NSA,IJ).EQ.1) THEN
      IJK=LI(I)+LK(KOSS-1)+J
      DX1=X(IJK)-X(IJK-NJ-1)
      DY1=Y(IJK)-Y(IJK-NJ-1)
      DZ1=Z(IJK)-Z(IJK-NJ-1)
      DX2=X(IJK-1)-X(IJK-NJ)
      DY2=Y(IJK-1)-Y(IJK-NJ)
      DZ2=Z(IJK-1)-Z(IJK-NJ)
      LW=IJK
      DELN=DNOB(NSA,IJ)
      CALL WALLBCOBST
      SU(IJK)=SU(IJK)+SUU
      SV(IJK)=SV(IJK)+SVU
      SW(IJK)=SW(IJK)+SWU
      AP(IJK)=AP(IJK)+SUP
      SPV(IJK)=SPV(IJK)+SVP
      SP(IJK)=SP(IJK)+SWP
      GENOB(NSA,IJ)=GENT
      SUEDOB(NSA,IJ)=SUED
      AT(IJK)=0.
      END IF
      END DO
      END DO
!#######################################################
      DO K=KOSS,KOES
      DO I=IOSS,IOES
      DO J=JOSS,JOES
      IJK=LK(K)+LI(I)+J
      SU(IJK)=0.
      SV(IJK)=0.
      SW(IJK)=0.
      AP(IJK)=GREAT
      SPV(IJK)=GREAT
      SP(IJK)=GREAT
      END DO
      END DO
      END DO

      END DO

      RETURN
      END
