!***********************************************************************
!
      SUBROUTINE BCSCALART
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
      INTEGER :: I, J, K, INP, INBC, LBHLPT, LBHLPC, &
                 NSA, IOSS, IOES, JOSS, JOES, KOSS, KOES, &
                 I1, IK, IJK, K1, JK, IJ
      REAL(PREC) :: Q_FLUX
!
!###############################################
!.....BOUNDARY CONDITIONS FOR TEMPERATURE
!###############################################
      IF(IDIR.EQ.IEN) THEN
!
!######################
!.....[ W E S T ]
!######################
!
      DO 600 K=2,NKM
      DO 600 J=2,NJM
      SUU=0.
      SUP=0.
!
      INP=LK(K)+LI(2)+J
      INBC=(J-1)*NK+K
      LBHLP=LBW(INBC+JKS)
      LBHLPT=LBWT(INBC+JKS)
!--------------------------------------
!.....INLET
!--------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,NJ,1,NIJ,FY,FZ,F1(INP-NJ),IEN,T)
!---------------------------------------------------
!.....WALL: thermally active [41] or adiabatic [40]
!.... [42] heat flux wall conditions: Q_flux [mK/s]
!---------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
              IF(LBHLPT.EQ.41) THEN
                 T(INP-NJ)=TWEST(J,K)
                 CALL WALLSC(INP,NJ,1,NIJ)
              END IF
              IF(LBHLPT.EQ.40) THEN
                CALL SYMSC(INP,NJ,1,NIJ,IEN,T)
              END IF
              IF(LBHLPT.EQ.42) THEN
                Q_flux=QFLXWEST(J,K)
                CALL WALLHFLX(INP,NJ,1,NIJ,Q_flux)
                T(INP-NJ)=T(INP)+T_HFLX
             END IF
!--------------------------------------
!.....SYMMETRY
!--------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,NJ,1,NIJ,IEN,T)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AW(INP)=0.

  600 CONTINUE
!
!######################
!.....[ E A S T ]
!######################
!
      DO 605 K=2,NKM
      DO 605 J=2,NJM
      SUU=0.
      SUP=0.
      INP=LK(K)+LI(NIM)+J
      INBC=(J-1)*NK+K
      LBHLP=LBE(INBC+JKS)
      LBHLPT=LBET(INBC+JKS)
!-------------------------------------
!.....INLET
!-------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-NJ,1,NIJ,FY,FZ,-F1(INP),IEN,T)
!--------------------------------------------------
!.....WALL: thermally active [41] or adiabatic [40]
!--------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
             IF(LBHLPT.EQ.41) THEN
                T(INP+NJ)=TEAST(J,K)
                CALL WALLSC(INP,-NJ,1,NIJ)
             END IF
             IF(LBHLPT.EQ.40) THEN
                CALL SYMSC(INP,-NJ,1,NIJ,IEN,T)
             END IF
             IF(LBHLPT.EQ.42) THEN
                Q_flux=QFLXEAST(J,K)
                CALL WALLHFLX(INP,-NJ,1,NIJ,Q_flux)
                T(INP+NJ)=T(INP)+T_HFLX
             END IF
!-------------------------------------
!.....SYMMETRY
!-------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-NJ,1,NIJ,IEN,T)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AE(INP)=0.

  605 CONTINUE
!
!######################
!.....[ S O U T H ]
!######################
!
      DO 610 K=2,NKM
      DO 610 I=2,NIM
      SUU=0.
      SUP=0.
!
      INP=LK(K)+LI(I)+2
      INBC=(I-1)*NK+K
      LBHLP=LBS(INBC+IKS)
      LBHLPT=LBST(INBC+IKS)
!-------------------------------------
!.....INLET
!-------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,1,NIJ,NJ,FZ,FX,F2(INP-1),IEN,T)
!--------------------------------------------------
!.....WALL: thermally active [41] or adiabatic [40]
!--------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
             IF(LBHLPT.EQ.41) THEN
                T(INP-1)=TSOUTH(I,K)
                CALL WALLSC(INP,1,NIJ,NJ)
             END IF
             IF(LBHLPT.EQ.40) THEN
                CALL SYMSC(INP,1,NIJ,NJ,IEN,T)
             END IF
             IF(LBHLPT.EQ.42) THEN
                Q_flux=QFLXSOUTH(I,K)
                CALL WALLHFLX(INP,1,NIJ,NJ,Q_flux)
                T(INP-1)=T(INP)+T_HFLX
             END IF
!-------------------------------------
!.....SYMMETRY
!-------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,1,NIJ,NJ,IEN,T)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AS(INP)=0.

  610 CONTINUE
!
!######################
!.....[ N O R T H ]
!######################
!
      DO 615 K=2,NKM
      DO 615 I=2,NIM
      SUU=0.
      SUP=0.
      INP=LK(K)+LI(I)+NJM
      INBC=(I-1)*NK+K
      LBHLP=LBN(INBC+IKS)
      LBHLPT=LBNT(INBC+IKS)
!------------------------------------
!.....INLET
!------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-1,NIJ,NJ,FZ,FX,-F2(INP),IEN,T)
!--------------------------------------------------
!.....WALL: thermally active [41] or adiabatic [40]
!--------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
             IF(LBHLPT.EQ.41) THEN
                T(INP+1)=TNORTH(I,K)
                CALL WALLSC(INP,-1,NIJ,NJ)
             END IF
             IF(LBHLPT.EQ.40) THEN
                CALL SYMSC(INP,-1,NIJ,NJ,IEN,T)
             END IF
             IF(LBHLPT.EQ.42) THEN
                Q_flux=QFLXNORTH(I,K)
                CALL WALLHFLX(INP,-1,NIJ,NJ,Q_flux)
                T(INP+1)=T(INP)+T_HFLX
             END IF
!------------------------------------
!.....SYMMETRY
!------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-1,NIJ,NJ,IEN,T)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AN(INP)=0.

  615 CONTINUE
!
!#####################
!.....[ B O T T O M ]
!#####################
!
      DO 620 I=2,NIM
      DO 620 J=2,NJM
      SUU=0.
      SUP=0.
!
      INP=LK(2)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBB(INBC+IJS)
      LBHLPT=LBBT(INBC+IJS)
!--------------------------------------
!.....INLET
!--------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,NIJ,NJ,1,FX,FY,F3(INP-NIJ),IEN,T)
!--------------------------------------------------
!.....WALL: thermally active [41] or adiabatic [40]
!.....      heat flux conditions [42] Q_flux [mK/s]
!--------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
             IF(LBHLPT.EQ.41) THEN
                T(INP-NIJ)=TBOTTOM(I,J)
                CALL WALLSC(INP,NIJ,NJ,1)
             END IF
             IF(LBHLPT.EQ.40) THEN
                CALL SYMSC(INP,NIJ,NJ,1,IEN,T)
             END IF
             IF(LBHLPT.EQ.42) THEN
                Q_flux=QFLXBOTTOM(I,J)
                CALL WALLHFLX(INP,NIJ,NJ,1,Q_flux)
                T(INP-NIJ)=T(INP)+T_HFLX
             END IF
!--------------------------------------
!.....SYMMETRY
!--------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,NIJ,NJ,1,IEN,T)
      ENDIF

!      IF(SUU.LT.0.OR.SUP.LT.O) THEN
!      WRITE(6,*) 'WARNING!!!: NEGATIVE COEFFICIENTS!!!!!!! '
!      WRITE(6,*) 'I= ',I,' J= ',J,' SUU= ',SUU,' SUP= ',SUP
!      END IF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AB(INP)=0.

  620 CONTINUE
!
!#####################
!.....[ T O P ]
!#####################
!
      DO 625 I=2,NIM
      DO 625 J=2,NJM
      SUU=0.
      SUP=0.
      INP=LK(NKM)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBT(INBC+IJS)
      LBHLPT=LBTT(INBC+IJS)
!------------------------------------
!.....INLET
!------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-NIJ,NJ,1,FX,FY,-F3(INP),IEN,T)
!--------------------------------------------------
!.....WALL: thermally active [41] or adiabatic [40]
!--------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
             IF(LBHLPT.EQ.41) THEN
                T(INP+NIJ)=TTOP(I,J)
                CALL WALLSC(INP,-NIJ,NJ,1)
             END IF
             IF(LBHLPT.EQ.40) THEN
                CALL SYMSC(INP,-NIJ,NJ,1,IEN,T)
             END IF
             IF(LBHLPT.EQ.42) THEN
                Q_flux=QFLXTOP(I,J)
                CALL WALLHFLX(INP,-NIJ,NJ,1,Q_flux)
                T(INP+NIJ)=T(INP)+T_HFLX
             END IF
!------------------------------------
!.....SYMMETRY
!------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-NIJ,NJ,1,IEN,T)
!      T(INP+NIJ)=TEMPERATURATOP
      CALL INLSC(INP,-NIJ,NJ,1,FX,FY,-F3(INP),IEN,T)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AT(INP)=0.

  625 CONTINUE
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
      SUP=0.
      IK=I1+K-KOSS+1
!----------------------
!     [adiabatic wall]
!----------------------
      IF(JOSS.NE.2.AND.TYPOBTS(NSA,IK).EQ.0) THEN
         IJK=LI(I)+LK(K)+JOSS-1
         T(IJK+1)=T(IJK)
         AN(IJK)=0.
      END IF
!----------------------------
!     [thermally active wall]
!----------------------------
      IF(JOSS.NE.2.AND.TYPOBTS(NSA,IK).EQ.1) THEN
         IJK=LI(I)+LK(K)+JOSS-1
         DX1=X(IJK)-X(IJK-NIJ-NJ)
         DY1=Y(IJK)-Y(IJK-NIJ-NJ)
         DZ1=Z(IJK)-Z(IJK-NIJ-NJ)
         DX2=X(IJK-NJ)-X(IJK-NIJ)
         DY2=Y(IJK-NJ)-Y(IJK-NIJ)
         DZ2=Z(IJK-NJ)-Z(IJK-NIJ)
         LW=IJK
         LWW=IJK+1
         T(IJK+1)=TOBSOUTH(NSA,IK)
         DELN=DNOS(NSA,IK)
         CALL WALLSCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AN(IJK)=0.
      END IF
!##########################################
!     [NORTH WALL OF SUBMERGED BODY: ]
!##########################################
      SUU=0.
      SUP=0.
!----------------------
!     [adiabatic wall]
!----------------------
      IF(JOES.NE.NJM.AND.TYPOBTN(NSA,IK).EQ.0) THEN
         IJK=LI(I)+LK(K)+JOES+1
         T(IJK-1)=T(IJK)
         AS(IJK)=0.
      END IF
!----------------------------
!     [thermally active wall]
!----------------------------
      IF(JOES.NE.NJM.AND.TYPOBTN(NSA,IK).EQ.1) THEN
         IJK=LI(I)+LK(K)+JOES+1
         DX1=X(IJK-NJ-1)-X(IJK-NIJ-1)
         DY1=Y(IJK-NJ-1)-Y(IJK-NIJ-1)
         DZ1=Z(IJK-NJ-1)-Z(IJK-NIJ-1)
         DX2=X(IJK-1)-X(IJK-NIJ-NJ-1)
         DY2=Y(IJK-1)-Y(IJK-NIJ-NJ-1)
         DZ2=Z(IJK-1)-Z(IJK-NIJ-NJ-1)
         LW=IJK
         LWW=IJK-1
         T(IJK-1)=TOBNORTH(NSA,IK)
         DELN=DNON(NSA,IK)
         CALL WALLSCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AS(IJK)=0.
      END IF
!----------------------------
      END DO
      END DO
!##########################################
!     [WEST WALL OF SUBMERGED BODY: ]
!##########################################
      DO K=KOSS,KOES
      K1=(K-KOSS)*(JOES-JOSS+1)
      DO J=JOSS,JOES
      SUU=0.
      SUP=0.
      JK=K1+J-JOSS+1
!----------------------
!     [adiabatic wall]
!----------------------
      IF(IOSS.NE.2.AND.TYPOBTW(NSA,JK).EQ.0) THEN
         IJK=LK(K)+LI(IOSS-1)+J
         T(IJK+NJ)=T(IJK)
         AE(IJK)=0.
      END IF
!----------------------------
!     [thermally active wall]
!----------------------------
      IF(IOSS.NE.2.AND.TYPOBTW(NSA,JK).EQ.1) THEN
         IJK=LK(K)+LI(IOSS-1)+J
         DX1=X(IJK-1)-X(IJK-NIJ)
         DY1=Y(IJK-1)-Y(IJK-NIJ)
         DZ1=Z(IJK-1)-Z(IJK-NIJ)
         DX2=X(IJK)-X(IJK-NIJ-1)
         DY2=Y(IJK)-Y(IJK-NIJ-1)
         DZ2=Z(IJK)-Z(IJK-NIJ-1)
         LW=IJK
         LWW=IJK+NJ
         T(IJK+NJ)=TOBWEST(NSA,JK)
         DELN=DNOW(NSA,JK)
         CALL WALLSCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AE(IJK)=0.
      END IF
!##########################################
!     [EAST WALL OF SUBMERGED BODY: ]
!##########################################
      SUU=0.
      SUP=0.
!----------------------
!     [adiabatic wall]
!----------------------
      IF(IOES.NE.NIM.AND.TYPOBTE(NSA,JK).EQ.0) THEN
         IJK=LK(K)+LI(IOES+1)+J
         T(IJK-NJ)=T(IJK)
         AW(IJK)=0.
      END IF
!----------------------------
!     [thermally active wall]
!----------------------------
      IF(IOES.NE.NIM.AND.TYPOBTE(NSA,JK).EQ.1) THEN
         IJK=LK(K)+LI(IOES+1)+J
         DX1=X(IJK-NJ)-X(IJK-NIJ-NJ-1)
         DY1=Y(IJK-NJ)-Y(IJK-NIJ-NJ-1)
         DZ1=Z(IJK-NJ)-Z(IJK-NIJ-NJ-1)
         DX2=X(IJK-NJ-1)-X(IJK-NIJ-NJ)
         DY2=Y(IJK-NJ-1)-Y(IJK-NIJ-NJ)
         DZ2=Z(IJK-NJ-1)-Z(IJK-NIJ-NJ)
         LW=IJK
         LWW=IJK-NJ
         T(IJK-NJ)=TOBEAST(NSA,JK)
         DELN=DNOE(NSA,JK)
         CALL WALLSCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AW(IJK)=0.
!----------------------------
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
      SUP=0.
      IJ=I1+J-JOSS+1
!----------------------
!     [adiabatic wall]
!----------------------
      IF(KOES.NE.NKM.AND.TYPOBTT(NSA,IJ).EQ.0) THEN
         IJK=LI(I)+LK(KOES+1)+J
         T(IJK-NIJ)=T(IJK)
         AB(IJK)=0.
      END IF
!----------------------------
!     [thermally active wall]
!----------------------------
      IF(KOES.NE.NKM.AND.TYPOBTT(NSA,IJ).EQ.1) THEN
         IJK=LI(I)+LK(KOES+1)+J
         DX1=X(IJK-NIJ-1)-X(IJK-NIJ-NJ)
         DY1=Y(IJK-NIJ-1)-Y(IJK-NIJ-NJ)
         DZ1=Z(IJK-NIJ-1)-Z(IJK-NIJ-NJ)
         DX2=X(IJK-NIJ)-X(IJK-NIJ-NJ-1)
         DY2=Y(IJK-NIJ)-Y(IJK-NIJ-NJ-1)
         DZ2=Z(IJK-NIJ)-Z(IJK-NIJ-NJ-1)
         LW=IJK
         LWW=IJK-NIJ
         T(IJK-NIJ)=TOBTOP(NSA,IJ)
         DELN=DNOT(NSA,IJ)
         CALL WALLSCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AB(IJK)=0.
      END IF
!----------------------------
!##########################################
!     [BOTTOM WALL OF SUBMERGED BODY: ]
!##########################################
      SUU=0.
      SUP=0.
!----------------------
!     [adiabatic wall]
!----------------------
      IF(KOSS.NE.2.AND.TYPOBTB(NSA,IJ).EQ.0) THEN
         IJK=LI(I)+LK(KOSS-1)+J
         T(IJK+NIJ)=T(IJK)
         AT(IJK)=0.
      END IF
!----------------------------
!     [thermally active wall]
!----------------------------
      IF(KOSS.NE.2.AND.TYPOBTB(NSA,IJ).EQ.1) THEN
         IJK=LI(I)+LK(KOSS-1)+J
         DX1=X(IJK)-X(IJK-NJ-1)
         DY1=Y(IJK)-Y(IJK-NJ-1)
         DZ1=Z(IJK)-Z(IJK-NJ-1)
         DX2=X(IJK-1)-X(IJK-NJ)
         DY2=Y(IJK-1)-Y(IJK-NJ)
         DZ2=Z(IJK-1)-Z(IJK-NJ)
         LW=IJK
         LWW=IJK+NIJ
         T(IJK+NIJ)=TOBBOTTOM(NSA,IJ)
         DELN=DNOB(NSA,IJ)
         CALL WALLSCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AT(IJK)=0.
      END IF
!----------------------------
      END DO
      END DO
!#######################################################
      DO K=KOSS,KOES
      DO J=JOSS,JOES
      DO I=IOSS,IOES
      IJK=LK(K)+LI(I)+J
      SU(IJK)=0.
      SP(IJK)=0.
      END DO
      END DO
      END DO
      END DO

      RETURN
!
!###################################################
!.....BOUNDARY CONDITIONS FOR TEMPERATURE VARIANCE
!###################################################
      ELSEIF(IDIR.EQ.IVART) THEN
!
!######################
!.....[ W E S T ]
!######################
!
      DO 630 K=2,NKM
      DO 630 J=2,NJM
      SUU=0.
      SUP=0.
!
      INP=LK(K)+LI(2)+J
      INBC=(J-1)*NK+K
      LBHLP=LBW(INBC+JKS)
!--------------------------------------
!.....INLET
!--------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,NJ,1,NIJ,FY,FZ,F1(INP-NJ),IVART,VART)
!--------------------------------------
!.....WALL: adiabatic or thermal active
!--------------------------------------
      ELSEIF(LBHLP.EQ.4.) THEN
      VART(INP-NJ)=0.
      SUU=0.
      SUP=0.
!--------------------------------------
!.....SYMMETRY
!--------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,NJ,1,NIJ,IVART,VART)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AW(INP)=0.

  630 CONTINUE
!
!######################
!.....[ E A S T ]
!######################
!
      DO 635 K=2,NKM
      DO 635 J=2,NJM
      SUU=0.
      SUP=0.
      INP=LK(K)+LI(NIM)+J
      INBC=(J-1)*NK+K
      LBHLP=LBE(INBC+JKS)
!-------------------------------------
!.....INLET
!-------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-NJ,1,NIJ,FY,FZ,-F1(INP),IVART,VART)
!--------------------------------------
!.....WALL: adiabatic or thermal active
!--------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
      VART(INP+NJ)=0.
      SUU=0.
      SUP=0.
!-------------------------------------
!.....SYMMETRY
!-------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-NJ,1,NIJ,IVART,VART)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AE(INP)=0.

  635 CONTINUE
!
!######################
!.....[ S O U T H ]
!######################
!
      DO 640 K=2,NKM
      DO 640 I=2,NIM
      SUU=0.
      SUP=0.
!
      INP=LK(K)+LI(I)+2
      INBC=(I-1)*NK+K
      LBHLP=LBS(INBC+IKS)
!-------------------------------------
!.....INLET
!-------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,1,NIJ,NJ,FZ,FX,F2(INP-1),IVART,VART)
!--------------------------------------
!.....WALL: adiabatic or thermal active
!--------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
      VART(INP-1)=0.
      SUU=0.
      SUP=0.
!-------------------------------------
!.....SYMMETRY
!-------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,1,NIJ,NJ,IVART,VART)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AS(INP)=0.

  640 CONTINUE
!
!######################
!.....[ N O R T H ]
!######################
!
      DO 645 K=2,NKM
      DO 645 I=2,NIM
      SUU=0.
      SUP=0.
      INP=LK(K)+LI(I)+NJM
      INBC=(I-1)*NK+K
      LBHLP=LBN(INBC+IKS)
!------------------------------------
!.....INLET
!------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-1,NIJ,NJ,FZ,FX,-F2(INP),IVART,VART)
!---------------------------------------
!.....WALL: adiabatic or thermal active
!---------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
      VART(INP+1)=0.
      SUU=0.
      SUP=0.
!------------------------------------
!.....SYMMETRY
!------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-1,NIJ,NJ,IVART,VART)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AN(INP)=0.

  645 CONTINUE
!
!#####################
!.....[ B O T T O M ]
!#####################
!
      DO 650 I=2,NIM
      DO 650 J=2,NJM
      SUU=0.
      SUP=0.
!
      INP=LK(2)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBB(INBC+IJS)
!--------------------------------------
!.....INLET
!--------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,NIJ,NJ,1,FX,FY,F3(INP-NIJ),IVART,VART)
!--------------------------------------
!.....WALL: adiabatic or thermal active
!--------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
      VART(INP-NIJ)=0.
      SUU=0.
      SUP=0.
!--------------------------------------
!.....SYMMETRY
!--------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,NIJ,NJ,1,IVART,VART)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AB(INP)=0.

  650 CONTINUE
!
!#####################
!.....[ T O P ]
!#####################
!
      DO 655 I=2,NIM
      DO 655 J=2,NJM
      SUU=0.
      SUP=0.
      INP=LK(NKM)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBT(INBC+IJS)
!------------------------------------
!.....INLET
!------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-NIJ,NJ,1,FX,FY,-F3(INP),IVART,VART)
!--------------------------------------
!.....WALL: adiabatic or thermal active
!--------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
      VART(INP+NIJ)=0.
      SUU=0.
      SUP=0.
!------------------------------------
!.....SYMMETRY
!------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-NIJ,NJ,1,IVART,VART)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AT(INP)=0.

  655 CONTINUE
!
!###################################################
!.....BOUNDARY CONDITIONS FOR CONCENTRATION
!###################################################
      ELSEIF(IDIR.EQ.ICON) THEN
!
!######################
!.....[ W E S T ]
!######################
!
      DO 660 K=2,NKM
      DO 660 J=2,NJM
      SUU=0.
      SUP=0.
!
      INP=LK(K)+LI(2)+J
      INBC=(J-1)*NK+K
      LBHLP=LBW(INBC+JKS)
      LBHLPC=LBWC(INBC+JKS)
!--------------------------------------
!.....INLET
!--------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,NJ,1,NIJ,FY,FZ,F1(INP-NJ),ICON,CON)
!---------------------------------------------------
!.....WALL: active [41] or passive [40]
!---------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
              IF(LBHLPC.EQ.41) THEN
                 CON(INP-NJ)=CONWEST(J,K)
                 CALL WALLSCC(INP,NJ,1,NIJ)
              END IF
              IF(LBHLPC.EQ.40) THEN
                CALL SYMSC(INP,NJ,1,NIJ,ICON,CON)
              END IF
!--------------------------------------
!.....SYMMETRY
!--------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,NJ,1,NIJ,ICON,CON)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AW(INP)=0.

  660 CONTINUE
!
!######################
!.....[ E A S T ]
!######################
!
      DO 665 K=2,NKM
      DO 665 J=2,NJM
      SUU=0.
      SUP=0.
      INP=LK(K)+LI(NIM)+J
      INBC=(J-1)*NK+K
      LBHLP=LBE(INBC+JKS)
      LBHLPC=LBEC(INBC+JKS)
!-------------------------------------
!.....INLET
!-------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-NJ,1,NIJ,FY,FZ,-F1(INP),ICON,CON)
!--------------------------------------------------
!.....WALL: active [41] or passive [40]
!--------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
             IF(LBHLPC.EQ.41) THEN
                CON(INP+NJ)=CONEAST(J,K)
                CALL WALLSCC(INP,-NJ,1,NIJ)
             END IF
             IF(LBHLPC.EQ.40) THEN
                CALL SYMSC(INP,-NJ,1,NIJ,ICON,CON)
             END IF
!-------------------------------------
!.....SYMMETRY
!-------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-NJ,1,NIJ,ICON,CON)
      ENDIF
!
      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AE(INP)=0.

  665 CONTINUE
!
!######################
!.....[ S O U T H ]
!######################
!
      DO 670 K=2,NKM
      DO 670 I=2,NIM
      SUU=0.
      SUP=0.
!
      INP=LK(K)+LI(I)+2
      INBC=(I-1)*NK+K
      LBHLP=LBS(INBC+IKS)
      LBHLPC=LBSC(INBC+IKS)
!-------------------------------------
!.....INLET
!-------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,1,NIJ,NJ,FZ,FX,F2(INP-1),ICON,CON)
!--------------------------------------------------
!.....WALL: active [41] or passive [40]
!--------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
             IF(LBHLPC.EQ.41) THEN
                CON(INP-1)=CONSOUTH(I,K)
                CALL WALLSCC(INP,1,NIJ,NJ)
             END IF
             IF(LBHLPC.EQ.40) THEN
                CALL SYMSC(INP,1,NIJ,NJ,ICON,CON)
             END IF
!-------------------------------------
!.....SYMMETRY
!-------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,1,NIJ,NJ,ICON,CON)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AS(INP)=0.

  670 CONTINUE
!
!######################
!.....[ N O R T H ]
!######################
!
      DO 675 K=2,NKM
      DO 675 I=2,NIM
      SUU=0.
      SUP=0.
      INP=LK(K)+LI(I)+NJM
      INBC=(I-1)*NK+K
      LBHLP=LBN(INBC+IKS)
      LBHLPC=LBNC(INBC+IKS)
!------------------------------------
!.....INLET
!------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-1,NIJ,NJ,FZ,FX,-F2(INP),ICON,CON)
!--------------------------------------------------
!.....WALL: active [41] or passive [40]
!--------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
             IF(LBHLPC.EQ.41) THEN
                CON(INP+1)=CONNORTH(I,K)
                CALL WALLSCC(INP,-1,NIJ,NJ)
             END IF
             IF(LBHLPC.EQ.40) THEN
                CALL SYMSC(INP,-1,NIJ,NJ,ICON,CON)
             END IF
!------------------------------------
!.....SYMMETRY
!------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-1,NIJ,NJ,ICON,CON)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AN(INP)=0.

  675 CONTINUE
!
!#####################
!.....[ B O T T O M ]
!#####################
!
      DO 680 I=2,NIM
      DO 680 J=2,NJM
      SUU=0.
      SUP=0.
!
      INP=LK(2)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBB(INBC+IJS)
      LBHLPC=LBBC(INBC+IJS)
!--------------------------------------
!.....INLET
!--------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,NIJ,NJ,1,FX,FY,F3(INP-NIJ),ICON,CON)
!--------------------------------------------------
!.....WALL: active [41] or passive [40]
!--------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
             IF(LBHLPC.EQ.41) THEN
                CON(INP-NIJ)=CONBOTTOM(I,J)
                CALL WALLSCC(INP,NIJ,NJ,1)
             END IF
             IF(LBHLPT.EQ.40) THEN
                CALL SYMSC(INP,NIJ,NJ,1,ICON,CON)
             END IF
!--------------------------------------
!.....SYMMETRY
!--------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,NIJ,NJ,1,ICON,CON)
      ENDIF

!      IF(SUU.LT.0.OR.SUP.LT.O) THEN
!      WRITE(6,*) 'WARNING!!!: NEGATIVE COEFFICIENTS!!!!!!! '
!      WRITE(6,*) 'I= ',I,' J= ',J,' SUU= ',SUU,' SUP= ',SUP
!      END IF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AB(INP)=0.

  680 CONTINUE
!
!#####################
!.....[ T O P ]
!#####################
!
      DO 685 I=2,NIM
      DO 685 J=2,NJM
      SUU=0.
      SUP=0.
      INP=LK(NKM)+LI(I)+J
      INBC=(I-1)*NJ+J
      LBHLP=LBT(INBC+IJS)
      LBHLPC=LBTC(INBC+IJS)
!------------------------------------
!.....INLET
!------------------------------------
      IF(LBHLP.EQ.1) THEN
      CALL INLSC(INP,-NIJ,NJ,1,FX,FY,-F3(INP),ICON,CON)
!--------------------------------------------------
!.....WALL: active [41] or passive [40]
!--------------------------------------------------
      ELSEIF(LBHLP.EQ.4) THEN
             IF(LBHLPC.EQ.41) THEN
                CON(INP+NIJ)=CONTOP(I,J)
                CALL WALLSCC(INP,-NIJ,NJ,1)
             END IF
             IF(LBHLPT.EQ.40) THEN
                CALL SYMSC(INP,-NIJ,NJ,1,ICON,CON)
             END IF
!------------------------------------
!.....SYMMETRY
!------------------------------------
      ELSEIF(LBHLP.EQ.3) THEN
      CALL SYMSC(INP,-NIJ,NJ,1,ICON,CON)
      ENDIF

      SU(INP)=SU(INP)+SUU
      SP(INP)=SP(INP)+SUP
      AT(INP)=0.

  685 CONTINUE
!
!-----------------------------------------
!     [point sources: ]
!-----------------------------------------
!      I=10
!      J=10
!      K=5
!      INP=LK(K)+LI(I)+J
!      write(6,*) 'inp: ',inp
!      SU(INP)=GREAT*1.
!      SP(INP)=GREAT
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
      SUP=0.
      IK=I1+K-KOSS+1
!----------------------
!     [passive wall]
!----------------------
      IF(JOSS.NE.2.AND.TYPOBCS(NSA,IK).EQ.0) THEN
         IJK=LI(I)+LK(K)+JOSS-1
         CON(IJK+1)=CON(IJK)
         AN(IJK)=0.
      END IF
!----------------------------
!     [active wall]
!----------------------------
      IF(JOSS.NE.2.AND.TYPOBCS(NSA,IK).EQ.1) THEN
         IJK=LI(I)+LK(K)+JOSS-1
         DX1=X(IJK)-X(IJK-NIJ-NJ)
         DY1=Y(IJK)-Y(IJK-NIJ-NJ)
         DZ1=Z(IJK)-Z(IJK-NIJ-NJ)
         DX2=X(IJK-NJ)-X(IJK-NIJ)
         DY2=Y(IJK-NJ)-Y(IJK-NIJ)
         DZ2=Z(IJK-NJ)-Z(IJK-NIJ)
         LW=IJK
         LWW=IJK+1
         CON(IJK+1)=COBSOUTH(NSA,IK)
         DELN=DNOS(NSA,IK)
         CALL WALLSCCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AN(IJK)=0.
      END IF
!##########################################
!     [NORTH WALL OF SUBMERGED BODY: ]
!##########################################
      SUU=0.
      SUP=0.
!----------------------
!     [passive wall]
!----------------------
      IF(JOES.NE.NJM.AND.TYPOBCN(NSA,IK).EQ.0) THEN
         IJK=LI(I)+LK(K)+JOES+1
         CON(IJK-1)=CON(IJK)
         AS(IJK)=0.
      END IF
!----------------------------
!     [active wall]
!----------------------------
      IF(JOES.NE.NJM.AND.TYPOBCN(NSA,IK).EQ.1) THEN
         IJK=LI(I)+LK(K)+JOES+1
         DX1=X(IJK-NJ-1)-X(IJK-NIJ-1)
         DY1=Y(IJK-NJ-1)-Y(IJK-NIJ-1)
         DZ1=Z(IJK-NJ-1)-Z(IJK-NIJ-1)
         DX2=X(IJK-1)-X(IJK-NIJ-NJ-1)
         DY2=Y(IJK-1)-Y(IJK-NIJ-NJ-1)
         DZ2=Z(IJK-1)-Z(IJK-NIJ-NJ-1)
         LW=IJK
         LWW=IJK-1
         CON(IJK-1)=COBNORTH(NSA,IK)
         DELN=DNON(NSA,IK)
         CALL WALLSCCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AS(IJK)=0.
      END IF
!----------------------------
      END DO
      END DO
!##########################################
!     [WEST WALL OF SUBMERGED BODY: ]
!##########################################
      DO K=KOSS,KOES
      K1=(K-KOSS)*(JOES-JOSS+1)
      DO J=JOSS,JOES
      SUU=0.
      SUP=0.
      JK=K1+J-JOSS+1
!----------------------
!     [passive wall]
!----------------------
      IF(IOSS.NE.2.AND.TYPOBCW(NSA,JK).EQ.0) THEN
         IJK=LK(K)+LI(IOSS-1)+J
         CON(IJK+NJ)=CON(IJK)
         AE(IJK)=0.
      END IF
!----------------------------
!     [active wall]
!----------------------------
      IF(IOSS.NE.2.AND.TYPOBCW(NSA,JK).EQ.1) THEN
         IJK=LK(K)+LI(IOSS-1)+J
         DX1=X(IJK-1)-X(IJK-NIJ)
         DY1=Y(IJK-1)-Y(IJK-NIJ)
         DZ1=Z(IJK-1)-Z(IJK-NIJ)
         DX2=X(IJK)-X(IJK-NIJ-1)
         DY2=Y(IJK)-Y(IJK-NIJ-1)
         DZ2=Z(IJK)-Z(IJK-NIJ-1)
         LW=IJK
         LWW=IJK+NJ
         CON(IJK+NJ)=COBWEST(NSA,JK)
         DELN=DNOW(NSA,JK)
         CALL WALLSCCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AE(IJK)=0.
      END IF
!##########################################
!     [EAST WALL OF SUBMERGED BODY: ]
!##########################################
      SUU=0.
      SUP=0.
!----------------------
!     [passive wall]
!----------------------
      IF(IOES.NE.NIM.AND.TYPOBCE(NSA,JK).EQ.0) THEN
         IJK=LK(K)+LI(IOES+1)+J
         CON(IJK-NJ)=CON(IJK)
         AW(IJK)=0.
      END IF
!----------------------------
!     [active wall]
!----------------------------
      IF(IOES.NE.NIM.AND.TYPOBCE(NSA,JK).EQ.1) THEN
         IJK=LK(K)+LI(IOES+1)+J
         DX1=X(IJK-NJ)-X(IJK-NIJ-NJ-1)
         DY1=Y(IJK-NJ)-Y(IJK-NIJ-NJ-1)
         DZ1=Z(IJK-NJ)-Z(IJK-NIJ-NJ-1)
         DX2=X(IJK-NJ-1)-X(IJK-NIJ-NJ)
         DY2=Y(IJK-NJ-1)-Y(IJK-NIJ-NJ)
         DZ2=Z(IJK-NJ-1)-Z(IJK-NIJ-NJ)
         LW=IJK
         LWW=IJK-NJ
         CON(IJK-NJ)=COBEAST(NSA,JK)
         DELN=DNOE(NSA,JK)
         CALL WALLSCCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AW(IJK)=0.
!----------------------------
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
      SUP=0.
      IJ=I1+J-JOSS+1
!----------------------
!     [passive wall]
!----------------------
      IF(KOES.NE.NKM.AND.TYPOBCT(NSA,IJ).EQ.0) THEN
         IJK=LI(I)+LK(KOES+1)+J
         CON(IJK-NIJ)=CON(IJK)
         AB(IJK)=0.
      END IF
!----------------------------
!     [ active wall]
!----------------------------
      IF(KOES.NE.NKM.AND.TYPOBCT(NSA,IJ).EQ.1) THEN
         IJK=LI(I)+LK(KOES+1)+J
         DX1=X(IJK-NIJ-1)-X(IJK-NIJ-NJ)
         DY1=Y(IJK-NIJ-1)-Y(IJK-NIJ-NJ)
         DZ1=Z(IJK-NIJ-1)-Z(IJK-NIJ-NJ)
         DX2=X(IJK-NIJ)-X(IJK-NIJ-NJ-1)
         DY2=Y(IJK-NIJ)-Y(IJK-NIJ-NJ-1)
         DZ2=Z(IJK-NIJ)-Z(IJK-NIJ-NJ-1)
         LW=IJK
         LWW=IJK-NIJ
         CON(IJK-NIJ)=COBTOP(NSA,IJ)
         DELN=DNOT(NSA,IJ)
         CALL WALLSCCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AB(IJK)=0.
      END IF
!----------------------------
!##########################################
!     [BOTTOM WALL OF SUBMERGED BODY: ]
!##########################################
      SUU=0.
      SUP=0.
!----------------------
!     [passive wall]
!----------------------
      IF(KOSS.NE.2.AND.TYPOBCB(NSA,IJ).EQ.0) THEN
         IJK=LI(I)+LK(KOSS-1)+J
         CON(IJK+NIJ)=CON(IJK)
         AT(IJK)=0.
      END IF
!----------------------------
!     [active wall]
!----------------------------
      IF(KOSS.NE.2.AND.TYPOBCB(NSA,IJ).EQ.1) THEN
         IJK=LI(I)+LK(KOSS-1)+J
         DX1=X(IJK)-X(IJK-NJ-1)
         DY1=Y(IJK)-Y(IJK-NJ-1)
         DZ1=Z(IJK)-Z(IJK-NJ-1)
         DX2=X(IJK-1)-X(IJK-NJ)
         DY2=Y(IJK-1)-Y(IJK-NJ)
         DZ2=Z(IJK-1)-Z(IJK-NJ)
         LW=IJK
         LWW=IJK+NIJ
         CON(IJK+NIJ)=COBBOTTOM(NSA,IJ)
         DELN=DNOB(NSA,IJ)
         CALL WALLSCCOBST
         SU(IJK)=SU(IJK)+SUU
         SP(IJK)=SP(IJK)+SUP
         AT(IJK)=0.
      END IF
!----------------------------
      END DO
      END DO
!#######################################################
      DO K=KOSS,KOES
      DO J=JOSS,JOES
      DO I=IOSS,IOES
      IJK=LK(K)+LI(I)+J
      SU(IJK)=0.
      SP(IJK)=0.
      END DO
      END DO
      END DO
      END DO
!#######################################################
!
      ENDIF
      RETURN
      END
