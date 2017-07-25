!***********************************************************************
!
      SUBROUTINE FLUXSCT(NIE,NJE,NKE,IDEW,IDNS,IDTB,FI,GRADFI,IFI, &
                         ARVX,ARVY,ARVZ, &
                         FIF,ACFE,ACFW,FCF)
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE COEF
      USE COEFB
      USE GEOMETRY
      USE VARIABLES
      USE BUOY
      USE TIME_MOD
      USE OBSTACLE
      USE GRADIENTS

      IMPLICIT NONE
!
!***********************************************************************
!

      INTEGER, INTENT(IN) :: NIE, NJE, NKE, IDEW, IDNS, IDTB, IFI
      REAL(PREC), DIMENSION(NXYZA) :: ARVX,ARVY,ARVZ
      REAL(PREC), DIMENSION(NXYZA) :: FIF
      REAL(PREC), DIMENSION(NXYZA) :: ACFE, ACFW
      REAL(PREC), DIMENSION(NXYZA) :: FCF, FI
      REAL(PREC), DIMENSION(3,NXYZA) :: GRADFI

! 
!     LOCAL VARIABLES
!
      INTEGER :: I, J, K, INP, &
                 INE, INS, INB,INBS
      INTEGER :: IJK
!@      INTEGER :: indx
      REAL(PREC) :: GAM,PRTR,FXE,FXP,ARE,VOLE, &
                    VISTE,VISOBE,GAME,DE, &
                    CE,CP,FII,FM
      REAL(PREC) :: XF,YF,ZF,XI,YI,ZI
      REAL(PREC) :: ARX,ARY,ARZ
      REAL(PREC) :: XPN,YPN,ZPN
      REAL(PREC) :: nxx,nyy,nzz,ixi1,ixi2,ixi3,dpn,costheta,costn
      REAL(PREC) :: fdfie,fdfii,fcfie,fcfii,ffic,suadd
      REAL(PREC) :: d2x,d2y,d2z,d1x,d1y,d1z
      REAL(PREC) :: dfixi,dfiyi,dfizi
      REAL(PREC) :: dfixj,dfiyj,dfizj

      GAM=GDS(IFI)
      PRTR=PRTINV(IFI)
!-------------------------------------------
!.....CALCULATE EAST,TOP,NORTH  CELL FACE
!-------------------------------------------
      DO 100 K=2,NKE
      DO 100 I=2,NIE
      DO 100 J=2,NJE
      INP=LK(K)+LI(I)+J
      INE=INP+IDEW
      INS=INP-IDNS
      INB=INP-IDTB
      INBS=INB-IDNS
!.....INTERPOLATION FACTORS IN FIRST,SECOND AND THIRD WAY
      FXE=FIF(INP)
      FXP=1.-FXE

!.....COMPONENTS OF THREE VECTORS
!.....FIRST
      XPN=XC(INE)-XC(INP)
      YPN=YC(INE)-YC(INP)
      ZPN=ZC(INE)-ZC(INP)

!.....Distance from P to neighbor N
      dpn=sqrt(XPN**2+YPN**2+ZPN**2)     
!.....Components of the unit vector i_ksi
      ixi1=xpn/dpn
      ixi2=ypn/dpn
      ixi3=zpn/dpn


!.....Precomputed face areas
      ARX=ARVX(INP)
      ARY=ARVY(INP)
      ARZ=ARVZ(INP)

!.....CELL FACE AREA
      ARE=SQRT(ARX**2+ARY**2+ARZ**2)

!.....Unit vectors of the normal
      nxx=arx/are
      nyy=ary/are
      nzz=arz/are
!.....Angle between vectorsa n and i_xi - we need cosine
      costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

!.....Relaxation factor for higher-order cell face gradient
      ! Minimal correction: nrelax = +1 :
      !costn = costheta
      ! Orthogonal correction: nrelax =  0 : 
      costn = 1.0d0
      ! Over-relaxed approach: nrelax = -1 :
      !costn = 1./costheta
      ! In general, nrelax can be any signed integer from some 
      ! reasonable interval [-nrelax,nrelax] (or even real number): 
      !costn = costheta**nrelax

!.....FIRST .(SECOND X THIRD) = VOL
      VOLE=XPN*ARX+YPN*ARY+ZPN*ARZ
!.....Overrelaxed correction vector d2, where S=dpn+d2
      d1x = costn
      d1y = costn
      d1z = costn
      !
      d2x = xpn*costn
      d2y = ypn*costn
      d2z = zpn*costn

!.....CELL FACE DIFFUSSION COEFFICIENT
      GAME=PRTR*(VIS(INP)*FXP+VIS(INE)*FXE)/VOLE

      VISTE=(VIS(INP)-VISOB(INP))*FXP+(VIS(INE)-VISOB(INE))*FXE
      VISOBE=VISOB(INP)*FXP+VISOB(INE)*FXE

      IF(IFI.EQ.IEN) GAME=(VISTE*PRT1+VISOBE*PRM1)
      IF(IFI.EQ.IVART)GAME=(VISTE*PRTR+VISOBE*PRM1)
!-----------------------------------------------------------
!     [Pr=nu/alfa; Sc=nu/beta; Le=Sc/Pr ]: Le=1.
!-----------------------------------------------------------
      IF(IFI.EQ.ICON) GAME=(VISTE*PRTR+VISOBE*PRM1)

!++++VALUES OF FI AT CELL FACE CENTER and EXPLICIT DIFFUSION FLUXES+++++++

!.....Coordinates of cell-face center
      XF = 0.25*(X(INP)+X(INS)+X(INB)+X(INBS))
      YF = 0.25*(Y(INP)+Y(INS)+Y(INB)+Y(INBS))
      ZF = 0.25*(Z(INP)+Z(INS)+Z(INB)+Z(INBS))

!.....Coordinates of point e'
      XI=XC(INP)*FXP+XC(INE)*FXE
      YI=YC(INP)*FXP+YC(INE)*FXE
      ZI=ZC(INP)*FXP+ZC(INE)*FXE

!.....Interpolate gradients defined at CV centers to faces
      DFIXI = GRADFI(1,INP)*FXP+GRADFI(1,INE)*FXE
      DFIYI = GRADFI(2,INP)*FXP+GRADFI(2,INE)*FXE
      DFIZI = GRADFI(3,INP)*FXP+GRADFI(3,INE)*FXE

!.....Interpolate variable FI defined at CV centers to face using corrected CDS:
!        |________Ue'_________|_______________Ucorr_________________|
      FII=FI(INP)*FXP+FI(INE)*FXE+DFIXI*(XF-XI)+DFIYI*(YF-YI)+DFIZI*(ZF-ZI)

!.....The cell face interpolated gradient (d phi / dx_i)_j:
      dfixj = dfixi*d1x + arx/vole*( fi(ine)-fi(inp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
      dfiyj = dfiyi*d1y + ary/vole*( fi(ine)-fi(inp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
      dfizj = dfizi*d1z + arz/vole*( fi(ine)-fi(inp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Explicit diffusion
      fdfie = game*(dfixj*arx + dfiyj*ary + dfizj*arz)   
      ! Implicit diffussion 
      fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++END: VALUES OF FI AT CELL FACE CENTER and EXPLICIT DIFFUSION FLUXES+++++++

!.....DIFUSION COEFFICIENT
      DE=GAME*ARE/DPN

!.....CONVECTION FLUXES - UDS
      FM=FCF(INP)
      CE=MIN(FM,ZERO) 
      CP=MAX(FM,ZERO)

!.....SYSTEM MATRIX COEFFOCIENTS
      ACFE(INP)=DE-CE
      ACFW(INE)=DE+CP
!
!---------------------------------------------
!     CONVECTIVE FLUXES [ CENTRAL DIFFERENCING SCHEME (CDS) ]
!---------------------------------------------
      FCFIE=FM*FII
      FCFII=CE*FI(INE)+CP*FI(INP)
      FFIC = GAM*(FCFIE-FCFII)
!-------------------------------------------------------
!.....EXPLICIT PART OF FLUXES
!-------------------------------------------------------
      SUADD=-FFIC+FDFIE-FDFII 
      SV(INP)=SV(INP)+SUADD
      SU(INP)=SU(INP)+SUADD
!-------------------------------------------------------
      BP(INE)=SUADD
  100 CONTINUE
!-------------------------------------------------------
!     [VECTORIZATION ROUTINE:  ]
!-------------------------------------------------------
      DO K=2,NKE
      DO I=2,NIE
      DO J=2,NJE
      INP=LK(K)+LI(I)+J+IDEW
      SU(INP)=SU(INP)-BP(INP)
      SV(INP)=SV(INP)-BP(INP)
      END DO !J-loop
      END DO !I-loop
      END DO !K-loop
!-------------
      DO IJK=ICST,ICEN
      BP(IJK)=0.D0
      END DO ! IJK-loop
!-------------
      RETURN
      END
