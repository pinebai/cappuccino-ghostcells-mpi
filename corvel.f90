      SUBROUTINE CORVEL
!#######################################################
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE VARIABLES
      USE BC
      USE BUOY
      USE TIME_MOD
      USE OBSTACLE
      
      IMPLICIT NONE

      INTEGER :: K, INP

!.....CORNER LINES
!.....IN X-WAY
!      DO 100 I=2,NIM
!      INP=LK(1)+LI(I)+1
!      U(INP)=U(INP+1)
!      TE(INP)=TE(INP+1)
!      ED(INP)=ED(INP+1)
!      T(INP)=T(INP+1)
!      VART(INP)=VART(INP+1)
!      VIS(INP)=VIS(INP+1)
!
!      INP=LK(1)+LI(I)+NJ
!      U(INP)=U(INP-1)
!      TE(INP)=TE(INP-1)
!      ED(INP)=ED(INP-1)
!      T(INP)=T(INP-1)
!      VART(INP)=VART(INP-1)
!      VIS(INP)=VIS(INP-1)
!
!      INP=LK(NK)+LI(I)+NJ
!      U(INP)=U(INP-1)
!      TE(INP)=TE(INP-1)
!      ED(INP)=ED(INP-1)
!      T(INP)=T(INP-1)
!      VART(INP)=T(INP-1)
!      VIS(INP)=VIS(INP-1)
!.....
!      INBC=(I-1)*NK+NKM
!     IF(LBS(INBC+IKS).EQ.3) THEN
!      INP=LK(NK)+LI(I)+1
!      U(INP)=U(INP-NIJ)
!      TE(INP)=TE(INP-NIJ)
!      ED(INP)=ED(INP-NIJ)
!      T(INP)=T(INP-NIJ)
!c     VIS(INP)=VIS(INP-NIJ)
!c     ENDIF
!  100 CONTINUE
!.....IN Y-WAY
!      DO 200 J=2,NJM
!      INP=LK(1)+LI(1)+J
!      U(INP)=U(INP+NIJ)
!      V(INP)=V(INP+NIJ)
!      TE(INP)=TE(INP+NIJ)
!      ED(INP)=ED(INP+NIJ)
!      T(INP)=T(INP+NIJ)
!c     VIS(INP)=VIS(INP+NIJ)
!.....
!      INP=LK(1)+LI(NI)+J
!      U(INP)=U(INP+NIJ)
!      V(INP)=V(INP+NIJ)
!      TE(INP)=TE(INP+NIJ)
!      ED(INP)=ED(INP+NIJ)
!      T(INP)=T(INP+NIJ)
!c     VIS(INP)=VIS(INP+NIJ)
!.....
!      INP=LK(NK)+LI(1)+J
!      U(INP)=U(INP-NIJ)
!      V(INP)=V(INP-NIJ)
!      TE(INP)=TE(INP-NIJ)
!      ED(INP)=ED(INP-NIJ)
!      T(INP)=T(INP-NIJ)
!c     VIS(INP)=VIS(INP-NIJ)
!.....
!      INP=LK(NK)+LI(NI)+J
!      U(INP)=U(INP-NIJ)
!      V(INP)=V(INP-NIJ)
!      TE(INP)=TE(INP-NIJ)
!      ED(INP)=ED(INP-NIJ)
!      T(INP)=T(INP-NIJ)
!c     VIS(INP)=VIS(INP-NIJ)
!  200 CONTINUE
!.....IN Z-WAY
      DO 300 K=2,NKM
      INP=LK(K)+LI(2)+1
      U(INP)=U(INP+1)
      W(INP)=W(INP+1)
      TE(INP)=TE(INP+1)
      ED(INP)=ED(INP+1)
      T(INP)=T(INP+1)
      VIS(INP)=VIS(INP+1)
      VART(INP)=VART(INP+1)
!
      INP=LK(K)+LI(2)+NJ
      U(INP)=U(INP-1)
      W(INP)=W(INP-1)
      TE(INP)=TE(INP-1)
      ED(INP)=ED(INP-1)
      T(INP)=T(INP-1)
      VART(INP)=VART(INP-1)
      VIS(INP)=VIS(INP-1)
!
      INP=LK(K)+LI(NIM)+1
      U(INP)=U(INP+1)
      W(INP)=W(INP+1)
      TE(INP)=TE(INP+1)
      ED(INP)=ED(INP+1)
      T(INP)=T(INP+1)
      VART(INP)=VART(INP+1)
      VIS(INP)=VIS(INP+1)
!
      INP=LK(K)+LI(NIM)+NJ
      U(INP)=U(INP-1)
      W(INP)=W(INP-1)
      TE(INP)=TE(INP-1)
      ED(INP)=ED(INP-1)
      T(INP)=T(INP-1)
      VART(INP)=VART(INP-1)
      VIS(INP)=VIS(INP-1)
  300 CONTINUE
!
!.....CORNERS
      INP=LK( 2)+LI( 2)+2
      U(INP)=U(INP+1)
      TE(INP)=TE(INP+1)
      ED(INP)=ED(INP+1)
      T(INP)=T(INP+1)
      VIS(INP)=VIS(INP+1)
!
      INP=LK( 2)+LI(NIM)+2
      U(INP)=U(INP+1)
      TE(INP)=TE(INP+1)
      ED(INP)=ED(INP+1)
      T(INP)=T(INP+1)
      VIS(INP)=VIS(INP+1)
!
      INP=LK( 2)+LI( 2)+NJM
      U(INP)=U(INP-1)
      TE(INP)=TE(INP-1)
      ED(INP)=ED(INP-1)
      T(INP)=T(INP-1)
      VIS(INP)=VIS(INP-1)
!
      INP=LK( 2)+LI(NIM)+NJM
      U(INP)=U(INP-1)
      TE(INP)=TE(INP-1)
      ED(INP)=ED(INP-1)
      T(INP)=T(INP-1)
      VIS(INP)=VIS(INP-1)
!.....TOP
      INP=LK(NKM)+LI( 2)+2
      U(INP)=U(INP+1)
      TE(INP)=TE(INP+1)
      ED(INP)=ED(INP+1)
      T(INP)=T(INP+1)
      VIS(INP)=VIS(INP+1)
!
      INP=LK(NKM)+LI(NIM)+2
      U(INP)=U(INP+1)
      TE(INP)=TE(INP+1)
      ED(INP)=ED(INP+1)
      T(INP)=T(INP+1)
!     VIS(INP)=VIS(INP+1)
!
      INP=LK(NKM)+LI( 2)+NJM
      U(INP)=U(INP-1)
      TE(INP)=TE(INP-1)
      ED(INP)=ED(INP-1)
      T(INP)=T(INP-1)
      VIS(INP)=VIS(INP-1)
!
      INP=LK(NKM)+LI(NIM)+NJM
      U(INP)=U(INP-1)
      TE(INP)=TE(INP-1)
      ED(INP)=ED(INP-1)
      T(INP)=T(INP-1)
      VIS(INP)=VIS(INP-1)
      RETURN
      END
