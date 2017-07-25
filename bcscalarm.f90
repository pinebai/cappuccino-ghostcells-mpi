!***********************************************************************
!
      SUBROUTINE BCSCALARM
!
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
      INTEGER :: I, J, K, INP, INBC, &
                 NSA, IOSS, IOES, JOSS, JOES, KOSS, KOES, &
                 I1, IK, IJK, K1, JK, IJ
      REAL(PREC) :: fac,are,arer,dxf,dyf,dzf,alf,bet,gam,cbeps,d1sq
!
!##########################################################
!.....BOUNDARY CONDITIONS FOR TURBULENT KINETIC ENERGY
!##########################################################
      IF(IDIR.EQ.ITE) THEN
!
!================
!.....[WEST : ]
!================
!
      DO 400 K=2,NKM
      DO 400 J=2,NJM
      SUU=0.
      SUP=0.
!
      INP=LK(K)+LI(2)+J
      INBC=(J-1)*NK+K
      LBHLP=LBW(INBC+JKS)
      IF(LBHLP.EQ.1) THEN
!-------------------
!.....INLET
!-------------------
if (periodic_boundary) cycle 
      CALL INLSC(INP,NJ,1,NIJ,FY,FZ,F1(INP-NJ),ITE,TE)
!-------------------
!.....WALL
!-------------------
      ELSEIF(LBHLP.EQ.4) THEN
!.....[High Re :]
      GEN(INP)=GENTW(INBC)
      SU(INP)=SV(INP)+GEN(INP)*VOL(INP)
!-------------------
!.....SYMMETRY
!-------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,NJ,1,NIJ,ITE,TE)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SUU
      SP(INP)=SP(INP)+SUP

      AW(INP)=0.

  400 CONTINUE
!
!================
!.....[EAST : ]
!================
!
      DO 405 K=2,NKM
      DO 405 J=2,NJM
      SUU=0.
      SUP=0.
      INP=LK(K)+LI(NIM)+J
      INBC=(J-1)*NK+K
      LBHLP=LBE(INBC+JKS)
      IF(LBHLP.EQ.1) THEN
!-------------------
!.....INLET
!-------------------
      CALL INLSC(INP,-NJ,1,NIJ,FY,FZ,-F1(INP),ITE,TE)
!-------------------
!.....WALL
!-------------------
      ELSEIF(LBHLP.EQ.4) THEN
      GEN(INP)=GENTE(INBC)
      SU(INP)=SV(INP)+GEN(INP)*VOL(INP)
!-------------------
!.....SYMMETRY
!-------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-NJ,1,NIJ,ITE,TE)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AE(INP)=0.

  405 CONTINUE
!
!===================
!.....[SOUTH : ]
!===================
!
      DO 410 K=2,NKM
      DO 410 I=2,NIM
      SUU=0.
      SUP=0.
!
      INP=LK(K)+LI(I)+2
      INBC=(I-1)*NK+K
      LBHLP=LBS(INBC+IKS)
!-------------------
!.....INLET
!-------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,1,NIJ,NJ,FZ,FX,F2(INP-1),ITE,TE)
!-------------------
!.....WALL
!-------------------
      ELSEIF(LBHLP.EQ.4) THEN
      GEN(INP)=GENTS(INBC)
      SU(INP)=SV(INP)+GEN(INP)*VOL(INP)
!-------------------
!.....SYMMETRY
!-------------------
      ELSEIF(LBHLP.EQ.3) THEN
if (periodic_boundary) cycle 
      CALL SYMSC(INP,1,NIJ,NJ,ITE,TE)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AS(INP)=0.

  410 CONTINUE
!
!==========================
!.....[NORTH : ]
!==========================
!
      DO 415 K=2,NKM
      DO 415 I=2,NIM
      SUU=0.
      SUP=0.
      INP=LK(K)+LI(I)+NJM
      INBC=(I-1)*NK+K
      LBHLP=LBN(INBC+IKS)
!-----------------------------
!.....INLET
!-----------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-1,NIJ,NJ,FZ,FX,-F2(INP),ITE,TE)
!-----------------------------
!.....WALL
!-----------------------------
      ELSEIF(LBHLP.EQ.4) THEN
      GEN(INP)=GENTN(INBC)
      SU(INP)=SV(INP)+GEN(INP)*VOL(INP)
!-----------------------------
!.....SYMMETRY
!-----------------------------
      ELSEIF(LBHLP.EQ.3) THEN
if (periodic_boundary) cycle 
      CALL SYMSC(INP,-1,NIJ,NJ,ITE,TE)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AN(INP)=0.

  415 CONTINUE
!
!==========================
!.....BOTTOM
!==========================
!
      DO 420 I=2,NIM
      DO 420 J=2,NJM
      SUU=0.
      SUP=0.
!
      INP=LK(2)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBB(INBC+IJS)
!-----------------------------
!.....INLET
!-----------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,NIJ,NJ,1,FX,FY,F3(INP-NIJ),ITE,TE)
!-----------------------------
!.....WALL
!-----------------------------
      ELSEIF(LBHLP.EQ.4) THEN
!+++++++++++++++++++++++++++++++++++++++++++++
!.....[STANDARD WALL FUNCTION: ]
!+++++++++++++++++++++++++++++++++++++++++++++
!      IF(.NOT.LowRe_LB)THEN
      GEN(INP)=GENTB(INBC)
      SU(INP)=SV(INP)+GEN(INP)*VOL(INP)
!      ENDIF
!+++++++++++++++++++++++++++++++++++++++++++++
!.....[NO WALL FUNCTION: ]
!+++++++++++++++++++++++++++++++++++++++++++++
!      IF(LowRe_LB)THEN
!      TE(INP-NIJ)=0.
!      SUU=0.
!      SUP=0.
!      ENDIF
!+++++++++++++++++++++++++++++++++++++++++++++
!.....[BUOYANCY WALL FUNCTION: ]
!+++++++++++++++++++++++++++++++++++++++++++++
!23456
!      DTDZ=(T(INP)-T(INP-NIJ))/(ZC(INP)-ZC(INP-NIJ))
!      PART1=DABS((VIS(INP)-VISCOS)*DTDZ/PRANT)
!      PART2=BETA*CAPPA*GRAVZ*(ZC(INP)-ZC(INP-NIJ))/CMU75
!      TE_FIRST=DABS(PART1*PART2)**(2/3.)
!      SUU=GREAT*TE_FIRST
!      SUP=GREAT

!      WRITE(6,*)'TE_FIRST: ',TE_FIRST

!
!-----------------------------
!.....SYMMETRY
!-----------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,NIJ,NJ,1,ITE,TE)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AB(INP)=0.

  420 CONTINUE
!
!=============================
!.....[TOP :]
!=============================
!
      DO 425 I=2,NIM
      DO 425 J=2,NJM
      SUU=0.
      SUP=0.
      INP=LK(NKM)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBT(INBC+IJS)
!------------------------------
!.....INLET
!------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-NIJ,NJ,1,FX,FY,-F3(INP),ITE,TE)
!------------------------------
!.....WALL
!------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
      GEN(INP)=GENTT(INBC)
      SU(INP)=SV(INP)+GEN(INP)*VOL(INP)
!+++++++++++++++++++++++++++++++++++++++++++++
!.....[NO WALL FUNCTION: ]
!+++++++++++++++++++++++++++++++++++++++++++++
!      IF(LowRe_LB)THEN
!      TE(INP+NIJ)=0.
!      SUU=0.
!      SUP=0.
!      ENDIF
!------------------------------
!.....SYMMETRY
!------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-NIJ,NJ,1,ITE,TE)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SV(INP)=SV(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AT(INP)=0.

  425 CONTINUE
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
      IK=I1+K-KOSS+1
      IF(JOSS.NE.2.AND.TYPOBS(NSA,IK).EQ.1) THEN
      IJK=LI(I)+LK(K)+JOSS-1
      SU(IJK)=GENOS(NSA,IK)*VOL(IJK)
      AN(IJK)=0.
      END IF
!##########################################
!     [NORTH WALL OF SUBMERGED BODY: ]
!##########################################
      IF(JOES.NE.NJM.AND.TYPOBN(NSA,IK).EQ.1) THEN
      IJK=LI(I)+LK(K)+JOES+1
      SU(IJK)=GENON(NSA,IK)*VOL(IJK)
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
      JK=K1+J-JOSS+1
      IF(IOSS.NE.2.AND.TYPOBW(NSA,JK).EQ.1) THEN
      IJK=LK(K)+LI(IOSS-1)+J
      SU(IJK)=GENOW(NSA,JK)*VOL(IJK)
      AE(IJK)=0.
!      write(6,*)'genow: ',genow(nsa,jk)
      END IF
!##########################################
!     [EAST WALL OF SUBMERGED BODY: ]
!##########################################
      IF(IOES.NE.NIM.AND.TYPOBE(NSA,JK).EQ.1) THEN
      IJK=LK(K)+LI(IOES+1)+J
      SU(IJK)=GENOE(NSA,JK)*VOL(IJK)
      AW(IJK)=0.
!      write(6,*)'genoe: ',genoe(nsa,jk)
      END IF
      END DO
      END DO
!##########################################
!     [TOP WALL OF SUBMERGED BODY: ]
!##########################################
      DO I=IOSS,IOES
      I1=(I-IOSS)*(JOES-JOSS+1)
      DO J=JOSS,JOES
      IJ=I1+J-JOSS+1
      IF(KOES.NE.NKM.AND.TYPOBT(NSA,IJ).EQ.1) THEN
      IJK=LI(I)+LK(KOES+1)+J
      SU(IJK)=GENOT(NSA,IJ)*VOL(IJK)
      AB(IJK)=0.
      END IF
!##########################################
!     [BOTTOM WALL OF SUBMERGED BODY: ]
!##########################################
      IF(KOSS.NE.2.AND.TYPOBB(NSA,IJ).EQ.1) THEN
      IJK=LI(I)+LK(KOSS-1)+J
      SU(IJK)=GENOB(NSA,IJ)*VOL(IJK)
      AT(IJK)=0.
!      write(6,*)'genob: ',genob(nsa,ij)
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
      SP(IJK)=0.
      END DO
      END DO
      END DO

      END DO
      RETURN
!
!###############################################################
!.....BOUNDARY CONDITIONS FOR DISSIPATION OF TURB.KIN.ENERGY
!###############################################################
      ELSEIF(IDIR.EQ.IED) THEN
!================
!.....WEST
!================
      DO 500 K=2,NKM
      DO 500 J=2,NJM
      SUU=0.
      SUP=0.
!
      INP=LK(K)+LI(2)+J
      INBC=(J-1)*NK+K
      LBHLP=LBW(INBC+JKS)
!------------------------------
!.....INLET
!------------------------------
      IF(LBHLP.EQ.1) THEN
if (periodic_boundary) cycle 
      CALL INLSC(INP,NJ,1,NIJ,FY,FZ,F1(INP-NJ),IED,ED)
!------------------------------
!.....WALL
!------------------------------
      ELSEIF(LBHLP.EQ.4.) THEN

      TE(INP)=DABS(TE(INP))
      SU(INP)=SUEDW(INBC)*TE(INP)*DSQRT(TE(INP))
      sp(inp)=1.0d0; ae(inp)=0.0d0; aw(inp)=0.0d0; an(inp)=0.0d0; as(inp)=0.0d0; at(inp)=0.0d0; ab(inp)=0.0d0
!------------------------------
!.....SYMMETRY
!------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,NJ,1,NIJ,IED,ED)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AW(INP)=0.

  500 CONTINUE
!
!===============================
!.....[E A S T : ]
!===============================
!
      DO 505 K=2,NKM
      DO 505 J=2,NJM
      SUU=0.
      SUP=0.
      INP=LK(K)+LI(NIM)+J
      INBC=(J-1)*NK+K
      LBHLP=LBE(INBC+JKS)
!---------------------------
!.....INLET
!---------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-NJ,1,NIJ,FY,FZ,-F1(INP),IED,ED)
!---------------------------
!.....WALL
!---------------------------
      ELSEIF(LBHLP.EQ.4) THEN
      TE(INP)=DABS(TE(INP))
      SU(INP)=SUEDE(INBC)*TE(INP)*DSQRT(TE(INP))
      sp(inp)=1.0d0; ae(inp)=0.0d0; aw(inp)=0.0d0; an(inp)=0.0d0; as(inp)=0.0d0; at(inp)=0.0d0; ab(inp)=0.0d0
!---------------------------
!.....SYMMETRY
!---------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-NJ,1,NIJ,IED,ED)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AE(INP)=0.

  505 CONTINUE
!
!=====================
!.....[S O U T H : ]
!=====================
!
      DO 510 K=2,NKM
      DO 510 I=2,NIM
      SUU=0.
      SUP=0.
!
      INP=LK(K)+LI(I)+2
      INBC=(I-1)*NK+K
      LBHLP=LBS(INBC+IKS)
!--------------------------
!.....INLET
!--------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,1,NIJ,NJ,FZ,FX,F2(INP-1),IED,ED)
!--------------------------
!.....WALL
!--------------------------
      ELSEIF(LBHLP.EQ.4) THEN
      TE(INP)=DABS(TE(INP))
      SU(INP)=SUEDS(INBC)*TE(INP)*DSQRT(TE(INP))
      sp(inp)=1.0d0; ae(inp)=0.0d0; aw(inp)=0.0d0; an(inp)=0.0d0; as(inp)=0.0d0; at(inp)=0.0d0; ab(inp)=0.0d0
!--------------------------
!.....SYMMETRY
!--------------------------
      ELSEIF(LBHLP.EQ.3) THEN
if (periodic_boundary) cycle 
      CALL SYMSC(INP,1,NIJ,NJ,IED,ED)
      ENDIF
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AS(INP)=0.
  510 CONTINUE
!
!===============================
!.....[N O R T H : ]
!===============================
!
      DO 515 K=2,NKM
      DO 515 I=2,NIM
      SUU=0.
      SUP=0.
      INP=LK(K)+LI(I)+NJM
      INBC=(I-1)*NK+K
      LBHLP=LBN(INBC+IKS)
!-------------------------------
!.....INLET
!-------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-1,NIJ,NJ,FZ,FX,-F2(INP),IED,ED)
!-------------------------------
!.....WALL
!-------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
      TE(INP)=DABS(TE(INP))
      SU(INP)=SUEDN(INBC)*TE(INP)*DSQRT(TE(INP))
      sp(inp)=1.0d0; ae(inp)=0.0d0; aw(inp)=0.0d0; an(inp)=0.0d0; as(inp)=0.0d0; at(inp)=0.0d0; ab(inp)=0.0d0
!-------------------------------
!.....SYMMETRY
!-------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
if (periodic_boundary) cycle 
      CALL SYMSC(INP,-1,NIJ,NJ,IED,ED)
      ENDIF
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AN(INP)=0.
  515 CONTINUE
!==========================
!.....[ B O T T O M : ]
!==========================
      DO 520 I=2,NIM
      DO 520 J=2,NJM
      SUU=0.
      SUP=0.
!
      INP=LK(2)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBB(INBC+IJS)
!---------------------------
!.....INLET
!---------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,NIJ,NJ,1,FX,FY,F3(INP-NIJ),IED,ED)
!---------------------------
!.....WALL
!---------------------------
      ELSEIF(LBHLP.EQ.4) THEN


      IF (Wilcox.or.SST.or.SAS.or.EARSM_WJ.or.EARSM_M) THEN ! models based on omega eq.
        call adaptive_wallbc_omega(inp,nij)
      ELSEIF(LowRe_LB) THEN 
        ! Lam-Bremhorst or any other Low-Re-k-epsilon model
        ! cmu*k^2/eps=k/omega; omega_wall_menter-sst=60*nu/(0.075*deln**2) => epsilon_wall_lowRE!
        ! or
        ! cmu*k^2/eps=k/omega; omega_wall_wicox=6*nu/(0.075*deln**2) => epsilon_wall_lowRE! 
        fac = sign(1.,float(-nij))
        dxf = fac*(xc(inp-nij)-xc(inp))
        dyf = fac*(yc(inp-nij)-yc(inp))
        dzf = fac*(zc(inp-nij)-zc(inp))
        are = dsqrt(ar1x(inp-nij)**2+ar1y(inp-nij)**2+ar1z(inp-nij)**2)
        arer = 1./are
        alf = ar1x(inp-nij)*arer
        bet = ar1y(inp-nij)*arer 
        gam = ar1z(inp-nij)*arer
        deln = dxf*alf+dyf*bet+dzf*gam
        d1sq = deln*deln
        cbeps = 6./(0.075*cmu)
        ED(INP)=cbeps*(viscos/den(inp))/(te(inp)*d1sq) ! <-? proveriti
        SUU=0.0d0
        SUP=0.0d0
      ELSE ! High-Re model
        TE(INP)=DABS(TE(INP))
        SU(INP)=SUEDB(INBC)*TE(INP)*DSQRT(TE(INP))
        sp(inp)=1.0d0; ae(inp)=0.0d0; aw(inp)=0.0d0; an(inp)=0.0d0; as(inp)=0.0d0; at(inp)=0.0d0; ab(inp)=0.0d0
      ENDIF

!---------------------------
!.....SYMMETRY
!---------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,NIJ,NJ,1,IED,ED)
      ENDIF
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AB(INP)=0.
  520 CONTINUE
!========================
!.....[ T O P : ]
!========================
      DO 525 I=2,NIM
      DO 525 J=2,NJM
      SUU=0.
      SUP=0.
      INP=LK(NKM)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBT(INBC+IJS)
!-------------------------------
!.....INLET
!-------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-NIJ,NJ,1,FX,FY,-F3(INP),IED,ED)
!-------------------------------
!.....WALL
!-------------------------------
      ELSEIF(LBHLP.EQ.4) THEN

      IF (Wilcox.or.SST.or.SAS.or.EARSM_WJ.or.EARSM_M) THEN ! models based on omega eq.
        call adaptive_wallbc_omega(inp,-nij)
      ELSEIF(LowRe_LB) THEN 
        ! Lam-Bremhorst or any other Low-Re-k-epsilon model
        ! cmu*k^2/eps=k/omega; omega_wall_menter-sst=60*nu/(0.075*deln**2) => epsilon_wall_lowRE!
        ! or
        ! cmu*k^2/eps=k/omega; omega_wall_wicox=6*nu/(0.075*deln**2) => epsilon_wall_lowRE! 
        fac = sign(1.,float(nij))
        dxf = fac*(xc(inp+nij)-xc(inp))
        dyf = fac*(yc(inp+nij)-yc(inp))
        dzf = fac*(zc(inp+nij)-zc(inp))
        are = dsqrt(ar1x(inp+nij)**2+ar1y(inp+nij)**2+ar1z(inp+nij)**2)
        arer = 1./are
        alf = ar1x(inp+nij)*arer
        bet = ar1y(inp+nij)*arer 
        gam = ar1z(inp+nij)*arer
        deln = dxf*alf+dyf*bet+dzf*gam
        d1sq = deln*deln
        cbeps = 6./(0.075*cmu)
        ED(INP)=cbeps*(viscos/den(inp))/(te(inp)*d1sq)
        SUU=0.0d0
        SUP=0.0d0
      ELSE ! it is a k-epsilon model
        TE(INP)=DABS(TE(INP))
        SU(INP)=SUEDT(INBC)*TE(INP)*DSQRT(TE(INP))
        sp(inp)=1.0d0; ae(inp)=0.0d0; aw(inp)=0.0d0; an(inp)=0.0d0; as(inp)=0.0d0; at(inp)=0.0d0; ab(inp)=0.0d0
      ENDIF
!-------------------------------
!.....SYMMETRY
!-------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-NIJ,NJ,1,IED,ED)
      ENDIF
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AT(INP)=0.
  525 CONTINUE
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
      IK=I1+K-KOSS+1
      IF(JOSS.NE.2.AND.TYPOBS(NSA,IK).EQ.1) THEN
      IJK=LI(I)+LK(K)+JOSS-1
      SU(IJK)=GREAT*SUEDOS(NSA,IK)*TE(IJK)*DSQRT(TE(IJK))
      SP(IJK)=GREAT
      AN(IJK)=0.
      END IF
!##########################################
!     [NORTH WALL OF SUBMERGED BODY: ]
!##########################################
      IF(JOES.NE.NJM.AND.TYPOBN(NSA,IK).EQ.1) THEN
      IJK=LI(I)+LK(K)+JOES+1
      SU(IJK)=GREAT*SUEDON(NSA,IK)*TE(IJK)*DSQRT(TE(IJK))
      SP(IJK)=GREAT
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
      JK=K1+J-JOSS+1
      IF(IOSS.NE.2.AND.TYPOBW(NSA,JK).EQ.1) THEN
      IJK=LK(K)+LI(IOSS-1)+J
      SU(IJK)=GREAT*SUEDOW(NSA,JK)*TE(IJK)*DSQRT(TE(IJK))
      SP(IJK)=GREAT
      AE(IJK)=0.
!      write(6,*)'suedow: ',suedow(nsa,jk)
      END IF
!##########################################
!     [EAST WALL OF SUBMERGED BODY: ]
!##########################################
      IF(IOES.NE.NIM.AND.TYPOBE(NSA,JK).EQ.1) THEN
      IJK=LK(K)+LI(IOES+1)+J
      SU(IJK)=GREAT*SUEDOE(NSA,JK)*TE(IJK)*DSQRT(TE(IJK))
      SP(IJK)=GREAT
      AW(IJK)=0.
!      write(6,*)'suedoe: ',suedoe(nsa,jk)
      END IF
      END DO
      END DO
!##########################################
!     [TOP WALL OF SUBMERGED BODY: ]
!##########################################
      DO I=IOSS,IOES
      I1=(I-IOSS)*(JOES-JOSS+1)
      DO J=JOSS,JOES
      IJ=I1+J-JOSS+1
      IF(KOES.NE.NKM.AND.TYPOBT(NSA,IJ).EQ.1) THEN
      IJK=LI(I)+LK(KOES+1)+J
      SU(IJK)=GREAT*SUEDOT(NSA,IJ)*TE(IJK)*DSQRT(TE(IJK))
      SP(IJK)=GREAT
      AB(IJK)=0.
      END IF
!##########################################
!     [BOTTOM WALL OF SUBMERGED BODY: ]
!##########################################
      IF(KOSS.NE.2.AND.TYPOBB(NSA,IJ).EQ.1) THEN
      IJK=LI(I)+LK(KOSS-1)+J
      SU(IJK)=GREAT*SUEDOB(NSA,IJ)*TE(IJK)*DSQRT(TE(IJK))
      SP(IJK)=GREAT
      AT(IJK)=0.
!      write(6,*)'suedob: ',suedob(nsa,ij)
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
      SP(IJK)=0.
      END DO
      END DO
      END DO

      END DO
      RETURN
!-----------------------------------------------

      ENDIF
      RETURN
      END
