!***********************************************************************
!
      SUBROUTINE FLUXMC(INP,IDEW,IDNS,IDTB, &
                        I,J,K,ARVX,ARVY,ARVZ, &
                        FIF,FMCOR)
!
!***********************************************************************
!
!     This routine calculates mass flux correction in the
!     second pressure-correction step which accounts for the
!     effects of non-orthogonality
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE COEF
      USE GEOMETRY
      USE VARIABLES
      USE GRADIENTS

      IMPLICIT NONE
!
!***********************************************************************
!
      INTEGER, INTENT(IN) :: INP,IDEW,IDNS,IDTB
      INTEGER, INTENT(IN) :: I,J,K
      REAL(PREC), DIMENSION(NXYZA) :: ARVX,ARVY,ARVZ
      REAL(PREC), DIMENSION(NXYZA), INTENT(IN) :: FIF
      REAL(PREC), INTENT(INOUT) :: FMCOR
!
!     LOCAL VARIABLES
!
      INTEGER :: INE,INS,INB,INBS
      INTEGER :: indx
      REAL(PREC) :: FXE,FXP
      !REAL(PREC) :: DXS,DYS,DZS,DXT,DYT,DZT
      REAL(PREC) :: ARE,ARX,ARY,ARZ, &
                    XPN,YPN,ZPN
      REAL(PREC) :: RAPR
      REAL(PREC) :: VOLE
      ! REAL(PREC) :: n_x,n_y,n_z, &
      !               X_f,Y_f,Z_f, &
      !               DPpP_x,DPpP_y,DPpP_z, &
      !               DNpN_x,DNpN_y,DNpN_z
      REAL(PREC) :: DPXI, DPYI, DPZI


      INE=INP+IDEW ! Neighbour CV

      INS=INP-IDNS
      INB=INP-IDTB
      INBS=INB-IDNS
!.....
      FXE=FIF(INP)
      FXP=1.0d0-FXE

!.....Precomputed face areas
      ARX=ARVX(INP)
      ARY=ARVY(INP)
      ARZ=ARVZ(INP)

!
!.....DISTANCE VECTOR COMPONENTS
      XPN=XC(INE)-XC(INP)
      YPN=YC(INE)-YC(INP)
      ZPN=ZC(INE)-ZC(INP)

!.....SURFACE VECTOR MAGNITUDE SQUARED
      ARE = sqrt(ARX**2+ARY**2+ARZ**2)

!+++++ANOTHER type correction due to skewness:+++++++++++++++++++++++++++++++++++
      VOLE = ARX*XPN+ARY*YPN+ARZ*ZPN
!.....CELL FACE COEFFICIENTS  1/AP(INP)
      RAPR = -0.5d0*(APU(INP)*DEN(INP)+APU(INE)*DEN(INE))
!.....MASS FLUX = OLD MASS FLUX + MASS FLUX CORRECTION
      DPXI = 0.5d0*(gradP(1,INE)+gradP(1,INP))
      DPYI = 0.5d0*(gradP(2,INE)+gradP(2,INP))
      DPZI = 0.5d0*(gradP(3,INE)+gradP(3,INP))
!.... o-------/-------o correction due to grid skewness
      FMCOR = RAPR*((VOLE*ARX-XPN*ARE)*DPXI &
                   +(VOLE*ARY-YPN*ARE)*DPYI &
                   +(VOLE*ARZ-ZPN*ARE)*DPZI)

! !.....Unit normal vector
!       n_x = ARX/ARE
!       n_y = ARY/ARE
!       n_z = ARZ/ARE

! !.....Coordinates of cell face center 'f'
!       X_f = 0.25*(X(INP)+X(INS)+X(INB)+X(INBS))
!       Y_f = 0.25*(Y(INP)+Y(INS)+Y(INB)+Y(INBS))
!       Z_f = 0.25*(Z(INP)+Z(INS)+Z(INB)+Z(INBS))

! !.....Distances form nodes P and N to auxiliary nodes P' and N'. $\Delta P'P = (r_{P'}-r_{P})$
!       DPpP_x = X_f - XC(INP) - (X_f-XC(INP))*n_x
!       DPpP_y = Y_f - YC(INP) - (Y_f-YC(INP))*n_y
!       DPpP_z = Z_f - ZC(INP) - (Z_f-ZC(INP))*n_z

!       DNpN_x = X_f - XC(INE) - (X_f-XC(INE))*n_x
!       DNpN_y = Y_f - YC(INE) - (Y_f-YC(INE))*n_y
!       DNpN_z = Z_f - ZC(INE) - (Z_f-ZC(INE))*n_z

! !.....APU==1./AP x density - interpolated            
!       RAPR = (APU(INP)*DEN(INP)*vol(inp)*FXP+APU(INE)*DEN(INE)*vol(ine)*FXE)


! !.....MASS FLUX CORRECTION FOR THE SECOND P'-EQUATION (SOURCE TERM)
!       FMCOR = RAPR*ARE*((gradP(1,INE)*DNpN_x-gradP(1,INP)*DPpP_x)   & 
!                        +(gradP(2,INE)*DNpN_y-gradP(2,INP)*DPpP_y)   &
!                        +(gradP(3,INE)*DNpN_z-gradP(3,INP)*DPpP_z))  &
!                        /((XPN*nx)+(YPN*ny)+(ZPN*nz)) 

! !.....Underrelax Mass-flux correction
!       !  fmcor = fmcor*0.6d0



      RETURN
      END
