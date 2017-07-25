      SUBROUTINE FLUXSCM(IDEW,IDNS,IDTB,FI,GRADFI,IFI, &
                         ARVX,ARVY,ARVZ, &
                         FIF,ACFE,ACFW,FCF)
!##########################################################
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE GEOMETRY
      USE COEF
      USE COEFB
      USE VARIABLES
      USE BUOY
      USE TIME_MOD
      USE GRADIENTS
      USE OMEGA_Turb_Models
      use fieldManipulation

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IDEW, IDNS, IDTB, IFI
      REAL(PREC), DIMENSION(NXYZA) :: ARVX,ARVY,ARVZ
      REAL(PREC), DIMENSION(NXYZA) :: FIF
      REAL(PREC), DIMENSION(NXYZA) :: ACFE, ACFW
      REAL(PREC), DIMENSION(NXYZA) :: FCF, FI
      REAL(PREC), DIMENSION(3,NXYZA) :: GRADFI

!     LOCAL VARIABLES

      INTEGER :: I, J, K, INP, INE, &
                 INS,INB,INBS,IJK
      REAL(PREC) :: GAM,PRTR,FXE,FXP,ARE,VOLE, &
                    VISTE,GAME,DE, &
                    CE,CP,FII,FM
      REAL(PREC) :: XF,YF,ZF,XI,YI,ZI
      REAL(PREC) :: ARX,ARY,ARZ
      REAL(PREC) :: XPN,YPN,ZPN
      REAL(PREC) :: nxx,nyy,nzz,ixi1,ixi2,ixi3,dpn,costheta,costn
      REAL(PREC) :: fdfie,fdfii,fcfie,fcfii,ffic,suadd
      REAL(PREC) :: d1x,d1y,d1z
      REAL(PREC) :: xpp,ypp,zpp,xep,yep,zep,xpnp,ypnp,zpnp,volep
      REAL(PREC) :: nablaFIxdnnp,nablaFIxdppp
      REAL(PREC) :: dfixi,dfiyi,dfizi
      REAL(PREC) :: dfixii,dfiyii,dfizii
      !@REAL(PREC) :: r1,r2,PSIE,PSIW
      !REAL(PREC) :: fxp1,fxe1


      GAM=GDS(IFI)

!.....Usually it is constant:
      PRTR=PRTINV(IFI)

      DFIXI = 0.0d0
      DFIYI = 0.0d0
      DFIZI = 0.0d0
!-------------------------------------------
!.....CALCULATE EAST,TOP,NORTH  CELL FACE
!-------------------------------------------
      DO K=2,NKM
      DO I=2,NIM
      DO J=2,NJM
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
      ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
      !costn = costheta**nrelax

!.....FIRST .(SECOND X THIRD) = VOL
      VOLE=XPN*ARX+YPN*ARY+ZPN*ARZ

!.....For Menter SST model:
      IF (SST.OR.SAS.OR.EARSM_WJ.OR.EARSM_M) THEN
        IF(IFI.EQ.ITE) PRTR=PRTINV_TE(INP)*fxp+PRTINV_TE(INE)*fxe
        IF(IFI.EQ.IED) PRTR=PRTINV_ED(INP)*fxp+PRTINV_ED(INE)*fxe
      END IF

!.....CELL FACE DIFFUSSION COEFFICINT
      VISTE=(VIS(INP)-VISCOS)*FXP+(VIS(INE)-VISCOS)*FXE
      GAME=(VISTE*PRTR+VISCOS)

!++++VALUES OF FI AT CELL FACE CENTER and EXPLICIT DIFFUSION FLUXES+++++++

!.....Coordinates of cell-face center
      XF = 0.25*(X(INP)+X(INS)+X(INB)+X(INBS))
      YF = 0.25*(Y(INP)+Y(INS)+Y(INB)+Y(INBS))
      ZF = 0.25*(Z(INP)+Z(INS)+Z(INB)+Z(INBS))

!.....Coordinates of point e'
      XI=XC(INP)*FXP+XC(INE)*FXE
      YI=YC(INP)*FXP+YC(INE)*FXE
      ZI=ZC(INP)*FXP+ZC(INE)*FXE


!.....Find points P' and E'
      xpp=xf-(xf-xc(inp))*nxx
      ypp=yf-(yf-yc(inp))*nyy
      zpp=zf-(zf-zc(inp))*nzz

      xep=xf-(xf-xc(ine))*nxx
      yep=yf-(yf-yc(ine))*nyy 
      zep=zf-(zf-zc(ine))*nzz     


      xpnp = xep-xpp
      ypnp = yep-ypp
      zpnp = zep-zpp

      volep = arx*xpnp+ary*ypnp+arz*zpnp

!.....Interpolation factors
      !fxp1 = sqrt((xep-xf)**2+(yep-yf)**2+(zep-zf)**2)/sqrt(xpnp*2+ypnp**2+zpnp**2)
      !fxe1 = 1.-fxp1

!.....Overrelaxed correction vector d2, where S=dpn+d2
      d1x = costn
      d1y = costn
      d1z = costn
      !
      xpnp = xpnp*costn
      ypnp = ypnp*costn
      zpnp = zpnp*costn

!.....Interpolate gradients defined at CV centers to faces
      DFIXI = GRADFI(1,INP)*FXP+GRADFI(1,INE)*FXE
      DFIYI = GRADFI(2,INP)*FXP+GRADFI(2,INE)*FXE
      DFIZI = GRADFI(3,INP)*FXP+GRADFI(3,INE)*FXE


!.....The cell face interpolated gradient (d phi / dx_i)_j:
      ! Nonorthogonal corrections:         ___
      ! nablaFIxdnnp =>> dot_product(gradU,dNN')
      ! And:                               ___
      ! nablaFIxdnnp =>> dot_product(gradU,dPP')
      !IF(IFI.EQ.ITE) THEN
      !nablaFIxdnnp = gradTE(1,ine)*(xep-xc(ine))+gradTE(2,ine)*(yep-yc(ine))+gradTE(3,ine)*(zep-zc(ine))
      !nablaFIxdppp = gradTE(1,inp)*(xpp-xc(inp))+gradTE(2,inp)*(ypp-yc(inp))+gradTE(3,inp)*(zpp-zc(inp))
      !ELSEIF(IFI.EQ.IED) THEN
      !nablaFIxdnnp = gradED(1,ine)*(xep-xc(ine))+gradED(2,ine)*(yep-yc(ine))+gradED(3,ine)*(zep-zc(ine))
      !nablaFIxdppp = gradED(1,inp)*(xpp-xc(inp))+gradED(2,inp)*(ypp-yc(inp))+gradED(3,inp)*(zpp-zc(inp))
      !ENDIF
      nablaFIxdnnp = gradFI(1,ine)*(xep-xc(ine))+gradFI(2,ine)*(yep-yc(ine))+gradFI(3,ine)*(zep-zc(ine))
      nablaFIxdppp = gradFI(1,inp)*(xpp-xc(inp))+gradFI(2,inp)*(ypp-yc(inp))+gradFI(3,inp)*(zpp-zc(inp))

      dfixii = dfixi*d1x + arx/volep*( fi(ine)+nablaFIxdnnp-fi(inp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
      dfiyii = dfiyi*d1y + ary/volep*( fi(ine)+nablaFIxdnnp-fi(inp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
      dfizii = dfizi*d1z + arz/volep*( fi(ine)+nablaFIxdnnp-fi(inp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Explicit diffusion
      fdfie = game*(dfixii*arx + dfiyii*ary + dfizii*arz)   
      ! Implicit diffussion 
      fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++END: VALUES OF FI AT CELL FACE CENTER and EXPLICIT DIFFUSION FLUXES+++++++-

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
!     [ CENTRAL DIFFERENCING SCHEME (CDS) ]
!---------------------------------------------
!.....Interpolate variable FI defined at CV centers to face using corrected CDS:
!         |________Ue'___________|_______________Ucorr_____________________|
!      FII=FI(INP)*FXP+FI(INE)*FXE+DFIXI*(XF-XI)+DFIYI*(YF-YI)+DFIZI*(ZF-ZI)
      FII = face_interpolated(FI,gradFI,inp,idew,idns,idtb,fxp,fxe)
!.....Another approach:
      !FII= (fi(inp)+nablaFixdppp)*fxp1 + (fi(ine)+nablaFIxdnnp)*fxe1

!.....Explicit second order convection
      FCFIE=FM*FII
!---------------------------------------------
!     Darwish-Moukalled TVD schemes for unstructured girds, IJHMT, 2003. 
!---------------------------------------------
!     Find r's - the gradient ratio. This is universal for all schemes.
!     If flow goes from P to E
      !@r1 = (2*gradFI(1,inp)*xpn + 2*gradFI(2,inp)*ypn + 2*gradFI(3,inp)*zpn)/(FI(INE)-FI(INP)) - 1.0d0  
!     If flow goes from E to P
      !@r2 = (2*gradFI(1,ine)*xpn + 2*gradFI(2,ine)*ypn + 2*gradFI(3,ine)*zpn)/(FI(INP)-FI(INE)) - 1.0d0 
!     Find Psi for [ MUSCL ] :
      !@PSIW = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
      !@PSIE = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
!     High order flux at cell face
      !@FCFIE =  CE*(FI(INE) + FXE*PSIE*(FI(INP)-FI(INE)))+ &
      !@         CP*(FI(INP) + FXP*PSIW*(FI(INE)-FI(INP)))

!.....Explicit first order convection
      FCFII=CE*FI(INE)+CP*FI(INP)
!.....Deffered correction for convection = gama_blending*(high-low)
      FFIC = GAM*(FCFIE-FCFII)
!-------------------------------------------------------
!.....EXPLICIT PART OF FLUXES
!-------------------------------------------------------
      SUADD=-FFIC+FDFIE-FDFII 
      SV(INP)=SV(INP)+SUADD
      SU(INP)=SU(INP)+SUADD
!-------------------------------------------------------
!      SV(INE)=SV(INE)-SUADD
!      SU(INE)=SU(INE)-SUADD
!-------------------------------------------------------
      BP(INE)=SUADD

      ENDDO
      ENDDO
      ENDDO

      DO K=2,NKM
      DO I=2,NIM
      DO J=2,NJM
      INP=LK(K)+LI(I)+J+IDEW

       SU(INP)=SU(INP)-BP(INP)
       SV(INP)=SV(INP)-BP(INP)

      END DO !J-loop
      END DO !I-loop
      END DO !K-loop

      DO IJK=ICST,ICEN
      BP(IJK)=0.0D0
      END DO

      RETURN
      END
