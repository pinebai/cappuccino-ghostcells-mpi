      SUBROUTINE GRAD_GAUSSc(U,DUDX,DUDY,DUDZ)
!
!***********************************************************************
!
!     Calculates cell centered gradient using Gauss theorem
!     PARAMETERS
!     U - field, the gradient of which we are looking for
!     DUDX,DUDY,DUDZ - arrays where the gradient components are stored
!
!     Gauss gradient rule:
!     ------->                                 ->
!     grad(u) = 1/vol * sum_{i=1}^{i=nf} (u)_f*Sf
!     where:
!     grad(u) - cell centered gradient vector
!     (u)_f   - face interpolated value of scalar u
!     vol     - cell volume
!     Sf      - cell face area vector
!     nf      - number of facec in a cell
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE GEOMETRY

      IMPLICIT NONE

      REAL(PREC), DIMENSION(NXYZA), INTENT(IN) :: U
      REAL(PREC), DIMENSION(NXYZA), INTENT(INOUT) :: DUDX,DUDY,DUDZ

      INTEGER :: I,J,K,INP,LC
      REAL(PREC) :: VOLR
      REAL(PREC), DIMENSION(NXYZA) :: DFXO,DFYO,DFZO
!
!.....INITIALIZE GRADIENT
      ! call GRAD_GAUSS(U,DFXO,DFYO,DFZO)
!
!.....START ITERATIVE CALCULATION OF GRADIENTS
      DO LC=1,NIGRAD
!
!.....INITIALIZE NEW GRADIENT
      dUdx(:)=0.0d0
      dUdy(:)=0.0d0
      dUdz(:)=0.0d0

!.....CALCULATE TERMS INTEGRATED OVER SURFACES
!.....ONLY INNER SURFACES

!.....EAST CELL - FACE
      CALL GRADCO(NIMM,NJMM,NKMM,NJ,1,NIJ, &
                  FX,AR1X,AR1Y,AR1Z,     &
                  U,DUDX,DUDY,DUDZ,DFXO,DFYO,DFZO)
!.....NORTH CELL - FACE
      CALL GRADCO(NIMM,NJMM,NKMM,1,NIJ,NJ, &
                  FY,AR2X,AR2Y,AR2Z,     &
                  U,DUDX,DUDY,DUDZ,DFXO,DFYO,DFZO)
!.....TOP   CELL - FACE
      CALL GRADCO(NIMM,NJMM,NKMM,NIJ,NJ,1, &
                  FZ,AR3X,AR3Y,AR3Z,    &
                  U,DUDX,DUDY,DUDZ,DFXO,DFYO,DFZO)
!
!.......CALCULATE GRADIENT COMPONENTS AT CV-CENTERS
      DO K=3,NKMM
      DO I=3,NIMM
      DO J=3,NJMM
        INP=LK(K)+LI(I)+J
        volr=1./vol(inp)
        dUdx(inp)=dUdx(inp)*volr
        dUdy(inp)=dUdy(inp)*volr
        dUdz(inp)=dUdz(inp)*volr
      ENDDO
      ENDDO
      ENDDO
!
!.......SET OLD GRADIENT = NEW GRADIENT FOR THE NEXT ITERATION
      IF(LC.NE.NIGRAD) THEN
        DFXO=dUdx
        DFYO=dUdy
        DFZO=dUdz
      ENDIF

      ENDDO ! LC-loop

      RETURN
      END

      SUBROUTINE GRADCO(NIE,NJE,NKE,IDEW,IDNS,IDTB, &
                        FIF,SX,SY,SZ, &
                        FI,DFX,DFY,DFZ,DFXO,DFYO,DFZO)
!=======================================================================
!
!=======================================================================
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE GEOMETRY

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NIE,NJE,NKE,IDEW,IDNS,IDTB 
      REAL(PREC), DIMENSION(NXYZA), INTENT(IN) :: FIF
      REAL(PREC), DIMENSION(NXYZA), INTENT(IN) :: FI
      REAL(PREC), DIMENSION(NXYZA), INTENT(IN) :: SX,SY,SZ
      REAL(PREC), DIMENSION(NXYZA)             :: DFX,DFY,DFZ
      REAL(PREC), DIMENSION(NXYZA), INTENT(IN) :: DFXO,DFYO,DFZO

      INTEGER :: I,J,K,INP
      INTEGER :: INE,INS,INB,INBS
      REAL(PREC) :: XI,YI,ZI,DFXI,DFYI,DFZI
      REAL(PREC) :: XF,YF,ZF
      REAL(PREC) :: FIE,DFXE,DFYE,DFZE
      REAL(PREC) :: FXE,FXP

      DO K=2,NKE
      DO I=2,NIE
      DO J=2,NJE
      INP=LK(K)+LI(I)+J

      INE=INP+IDEW
      INS=INP-IDNS
      INB=INP-IDTB 
      INBS=INB-IDNS

!
!.....Coordinates of point on the line connecting center and neighbor,
!     old gradient vector components interpolated for this location.
      FXE=FIF(INP) 
      FXP=1.0d0-FXE
      XI=XC(INP)*FXP+XC(INE)*FXE
      YI=YC(INP)*FXP+YC(INE)*FXE
      ZI=ZC(INP)*FXP+ZC(INE)*FXE
      DFXI=DFXO(INP)*FXP+DFXO(INE)*FXE
      DFYI=DFYO(INP)*FXP+DFYO(INE)*FXE
      DFZI=DFZO(INP)*FXP+DFZO(INE)*FXE

!
!.....Coordinates of the cell-face center
      XF=.25*(X(INP)+X(INS)+X(INB)+X(INBS))
      YF=.25*(Y(INP)+Y(INS)+Y(INB)+Y(INBS))
      ZF=.25*(Z(INP)+Z(INS)+Z(INB)+Z(INBS))

!
!.....Value of the variable at cell-face center
      FIE=FI(INP)*FXP+FI(INE)*FXE+DFXI*(XF-XI)+DFYI*(YF-YI)+DFZI*(ZF-ZI)

!.....(Interpolated mid-face value)x(Area)
      DFXE=FIE*SX(INP)
      DFYE=FIE*SY(INP)
      DFZE=FIE*SZ(INP)
!
!.....ACCUMULATE CONTRIBUTION AT CELL CENTER AND NEIGHBOR
      DFX(INP)=DFX(INP)+DFXE
      DFY(INP)=DFY(INP)+DFYE
      DFZ(INP)=DFZ(INP)+DFZE
!--
      DFX(INE)=DFX(INE)-DFXE
      DFY(INE)=DFY(INE)-DFYE
      DFZ(INE)=DFZ(INE)-DFZE

      END DO !!I-loop
      END DO !!J-loop
      END DO !!K-loop

      RETURN
      END
