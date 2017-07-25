!***********************************************************************
!
      SUBROUTINE FLUXUVW(IDEW,IDNS,IDTB, &
                         ARVX,ARVY,ARVZ, &
                         FIF,ACFE,ACFW,FCF)
!
!***********************************************************************
!
!  Fluid flow trough control volumes:
!   _______________________________________________________________
!  |           |            |            |            |            |
!  |           |            |            |            |            |
!  |     o WW ===>   o W   ===>   o P   ===>   o E   ===>   o  EE  |
!  |           |            |            |            |            |
!  |           |            |            |            |            |
!  |___________|____________|____________|____________|____________|
!
!
!***********************************************************************
!
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
      use fieldManipulation

      IMPLICIT NONE
!
!***********************************************************************
!
      INTEGER, INTENT(IN) :: IDEW,IDNS,IDTB 
      REAL(PREC), DIMENSION(NXYZA) :: FIF,ACFE,ACFW,FCF 
      REAL(PREC), DIMENSION(NXYZA) :: ARVX,ARVY,ARVZ
!
!     LOCAL VARIABLES
!
      INTEGER :: I, J, K, IJK, INP,       &
                 INE,INW,INN,INS,INB,INT, &
                 INBS,INEE
      REAL(PREC) :: FXW,FXPW,FXEE,FXEW
      ! REAL(PREC) :: DXU,DYU,DZU,DXU1,DYU1,DZU1,Uinu,Vinu,Winu,Uinu1,Vinu1,Winu1
      REAL(PREC) :: GAM,FXE,FXP,ARE,VOLE,GAME,&
                    UE,VE,WE, & 
                    DE,FLCF,CE,CP, &
                    G11, G12, G21, G22
      REAL(PREC) :: R1,R2,R3,R4,R5,R6,                   &
                    PSIE1,PSIE2,PSIE3,PSIW1,PSIW2,PSIW3, &
                    FUUDS,FVUDS,FWUDS,FUHIGH,FVHIGH,FWHIGH
      REAL(PREC) :: XF,YF,ZF,XI,YI,ZI,ARX,ARY,ARZ, &
                    DUXI,DUYI,DUZI,DVXI,DVYI,DVZI, &
                    DWXI,DWYI,DWZI
      REAL(PREC) :: onethird, twothirds
      REAL(PREC) :: nxx,nyy,nzz,xpp,ypp,zpp,xep,yep,zep,xpnp,ypnp,zpnp,volep
      REAL(PREC) :: ixi1,ixi2,ixi3
      REAL(PREC) :: costheta,d1x,d1y,d1z,costn
      REAL(PREC) :: xpn,ypn,zpn,dpn  
      REAL(PREC) :: nablaFIxdnnp,nablaFIxdppp
      REAL(PREC) :: fdue,fdve,fdwe,fdui,fdvi,fdwi
      REAL(PREC) :: duxii,dvxii,dwxii, &
                    duyii,dvyii,dwyii, &
                    duzii,dvzii,dwzii
      ! REAL(PREC) :: fxp1,fxe1
    
!----------------------------------------------------------------------
      GAM=GDS(IU) 
      onethird = 1./3.
      twothirds = 2./3.

!.....INITIALIZE VALUES FOR HIGH ORDER CONV. FLUXES
      PSIE1=0.0D0;PSIE2=0.0D0;PSIE3=0.0D0
      PSIW1=0.0D0;PSIW2=0.0D0;PSIW3=0.0D0         
      FUHIGH=0.0D0;FVHIGH=0.0D0;FWHIGH=0.0D0
!
!.....CALCULATE EAST,TOP,NORTH  CELL FACE
      DO K=2,NKMm
      DO I=2,NIMm
      DO J=2,NJMm
      INP=LK(K)+LI(I)+J
      INE=INP+IDEW
      INW=INP-IDEW
      INN=INP+IDNS
      INS=INP-IDNS
      INB=INP-IDTB
      INT=INP+IDTB

      INEE=INP+2*IDEW
      INBS=INB-IDNS
!
!.....INTERPOLATION FACTORS IN FIRST,SECOND AND THIRD DIRECTION
      FXE=FIF(INP) 
      FXP=1.-FXE

      FXPW = FIF(INW)
      FXW = 1.-FXPW
      FXEE = FIF(INE)
      FXEW = 1.-FXEE 

      XPN=XC(INE)-XC(INP)
      YPN=YC(INE)-YC(INP)
      ZPN=ZC(INE)-ZC(INP)

      dpn=sqrt(xpn**2+ypn**2+zpn**2)
!.....Components of the unit vector i_ksi which passes trough P and E
      ixi1=xpn/dpn
      ixi2=ypn/dpn
      ixi3=zpn/dpn

!.....Precomputed face areas
      ARX=ARVX(INP)
      ARY=ARVY(INP)
      ARZ=ARVZ(INP)

!.....CELL FACE AREA
      ARE=sqrt(ARX**2+ARY**2+ARZ**2)
!.....Unit vectors of the Surface normal
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

!.....CELL FACE
      GAME=(VIS(INP)*FXP+VIS(INE)*FXE)

!++++VELOCITIES AT CELL FACE CENTER and EXPLICIT DIFFUSION FLUXES+++++++

!.....Coordinates of cell-face center
      XF = 0.25d0*(X(INP)+X(INS)+X(INB)+X(INBS))
      YF = 0.25d0*(Y(INP)+Y(INS)+Y(INB)+Y(INBS))
      ZF = 0.25d0*(Z(INP)+Z(INS)+Z(INB)+Z(INBS))

!.....Coordinates of point e'
      XI=XC(INP)*FXP+XC(INE)*FXE
      YI=YC(INP)*FXP+YC(INE)*FXE
      ZI=ZC(INP)*FXP+ZC(INE)*FXE

!.....Find points P' and E'
      xpp=xf-(xf-xc(inp))*nxx; ypp=yf-(yf-yc(inp))*nyy; zpp=zf-(zf-zc(inp))*nzz
      xep=xf-(xf-xc(ine))*nxx; yep=yf-(yf-yc(ine))*nyy; zep=zf-(zf-zc(ine))*nzz  
   
      xpnp = xep-xpp; ypnp = yep-ypp; zpnp = zep-zpp
      volep = arx*xpnp+ary*ypnp+arz*zpnp

!.....Interpolation factors
      ! fxp1 = sqrt((xep-xf)**2+(yep-yf)**2+(zep-zf)**2)/sqrt(xpnp*2+ypnp**2+zpnp**2)
      ! fxe1 = 1.-fxp1

!.....Overrelaxed correction vector d2, where S=dpn+d2
      d1x = costn
      d1y = costn
      d1z = costn
      !
      xpnp = xpnp*costn
      ypnp = ypnp*costn
      zpnp = zpnp*costn

!..... U - velocity component .........................................................
!.....Interpolate gradients defined at CV centers to faces
      DUXI = GRADU(1,INP)*FXP+GRADU(1,INE)*FXE
      DUYI = GRADU(2,INP)*FXP+GRADU(2,INE)*FXE
      DUZI = GRADU(3,INP)*FXP+GRADU(3,INE)*FXE

!.....du/dx_i interpolated at cell face:
      ! Nonorthogonal corrections:         ___
      ! nablaFIxdnnp =>> dot_product(gradU,dNN')
      ! And:                               ___
      ! nablaFIxdnnp =>> dot_product(gradU,dPP')
      nablaFIxdnnp = gradU(1,ine)*(xep-xc(ine))+gradU(2,ine)*(yep-yc(ine))+gradU(3,ine)*(zep-zc(ine))
      nablaFIxdppp = gradU(1,inp)*(xpp-xc(inp))+gradU(2,inp)*(ypp-yc(inp))+gradU(3,inp)*(zpp-zc(inp))

      duxii = duxi*d1x + arx/volep*( u(ine)+nablaFIxdnnp-u(inp)-nablaFixdppp-duxi*xpnp-duyi*ypnp-duzi*zpnp ) 
      duyii = duyi*d1y + ary/volep*( u(ine)+nablaFIxdnnp-u(inp)-nablaFixdppp-duxi*xpnp-duyi*ypnp-duzi*zpnp ) 
      duzii = duzi*d1z + arz/volep*( u(ine)+nablaFIxdnnp-u(inp)-nablaFixdppp-duxi*xpnp-duyi*ypnp-duzi*zpnp ) 

!.....face interpolated value CDS - nonorthogonality corrected
!        |________Ue'_________|_______________Ucorr___________________|
      UE=U(INP)*FXP+U(INE)*FXE!+(DUXI*(XF-XI)+DUYI*(YF-YI)+DUZI*(ZF-ZI))
      ! UE = face_interpolated(U,gradU,inp,idew,idns,idtb,fxp,fxe)
!.....Another approach:
      ! ue = (u(inp)+nablaFIxdppp)*fxp1 + (u(ine)+nablaFIxdnnp)*fxe1

!..... V - velocity component .........................................................
      DVXI = GRADV(1,INP)*FXP+GRADV(1,INE)*FXE
      DVYI = GRADV(2,INP)*FXP+GRADV(2,INE)*FXE
      DVZI = GRADV(3,INP)*FXP+GRADV(3,INE)*FXE

!.....dv/dx_i interpolated at cell face:
      nablaFIxdnnp = gradV(1,ine)*(xep-xc(ine))+gradV(2,ine)*(yep-yc(ine))+gradV(3,ine)*(zep-zc(ine))
      nablaFIxdppp = gradV(1,inp)*(xpp-xc(inp))+gradV(2,inp)*(ypp-yc(inp))+gradV(3,inp)*(zpp-zc(inp))

      dvxii = dvxi*d1x + arx/volep*( v(ine)+nablaFIxdnnp-v(inp)-nablaFixdppp-dvxi*xpnp-dvyi*ypnp-dvzi*zpnp ) 
      dvyii = dvyi*d1y + ary/volep*( v(ine)+nablaFIxdnnp-v(inp)-nablaFixdppp-dvxi*xpnp-dvyi*ypnp-dvzi*zpnp ) 
      dvzii = dvzi*d1z + arz/volep*( v(ine)+nablaFIxdnnp-v(inp)-nablaFixdppp-dvxi*xpnp-dvyi*ypnp-dvzi*zpnp ) 

!.....face interpolated value CDS - nonorthogonality corrected
!        |________Ve'_________|_______________Vcorr___________________|
      VE=V(INP)*FXP+V(INE)*FXE!+(DVXI*(XF-XI)+DVYI*(YF-YI)+DVZI*(ZF-ZI))
      ! VE = face_interpolated(V,gradV,inp,idew,idns,idtb,fxp,fxe)
!.....Another approach:
      ! ve = (v(inp)+nablaFIxdppp)*fxp1 + (v(ine)+nablaFIxdnnp)*fxe1

!..... W - velocity component .........................................................
      DWXI = GRADW(1,INP)*FXP+GRADW(1,INE)*FXE
      DWYI = GRADW(2,INP)*FXP+GRADW(2,INE)*FXE
      DWZI = GRADW(3,INP)*FXP+GRADW(3,INE)*FXE

!.....dw/dx_i interpolated at cell face:
      nablaFIxdnnp = gradW(1,ine)*(xep-xc(ine))+gradW(2,ine)*(yep-yc(ine))+gradW(3,ine)*(zep-zc(ine))
      nablaFIxdppp = gradW(1,inp)*(xpp-xc(inp))+gradW(2,inp)*(ypp-yc(inp))+gradW(3,inp)*(zpp-zc(inp))

      dwxii = dwxi*d1x + arx/volep*( w(ine)+nablaFIxdnnp-w(inp)-nablaFixdppp-dwxi*xpnp-dwyi*ypnp-dwzi*zpnp ) 
      dwyii = dwyi*d1y + ary/volep*( w(ine)+nablaFIxdnnp-w(inp)-nablaFixdppp-dwxi*xpnp-dwyi*ypnp-dwzi*zpnp ) 
      dwzii = dwzi*d1z + arz/volep*( w(ine)+nablaFIxdnnp-w(inp)-nablaFixdppp-dwxi*xpnp-dwyi*ypnp-dwzi*zpnp ) 

!.....face interpolated value CDS - nonorthogonality corrected
!        |________We'_________|_______________Wcorr___________________|
      WE=W(INP)*FXP+W(INE)*FXE!+(DWXI*(XF-XI)+DWYI*(YF-YI)+DWZI*(ZF-ZI))
      ! WE = face_interpolated(W,gradW,inp,idew,idns,idtb,fxp,fxe)
!.....Another approach:
      ! we = (w(inp)+nablaFIxdppp)*fxp1 + (w(ine)+nablaFIxdnnp)*fxe1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     We calculate explicit and implicit diffsion fde and fdi,
!     later se put their difference (fde-fdi) to RHS vector
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! Explicit diffussion 
       fdue = (dUxii+dUxii)*arx + (dUyii+dVxii)*ary + (dUzii+dWxii)*arz
       fdve = (dUyii+dVxii)*arx + (dVyii+dVyii)*ary + (dVzii+dWyii)*arz
       fdwe = (dUzii+dWxii)*arx + (dWyii+dVzii)*ary + (dWzii+dWzii)*arz

       fdue = game*fdue
       fdve = game*fdve
       fdwe = game*fdwe
       
       ! Implicit diffussion 
       fdui = game*are/dpn*(duxi*xpn+duyi*ypn+duzi*zpn)
       fdvi = game*are/dpn*(dvxi*xpn+dvyi*ypn+dvzi*zpn)
       fdwi = game*are/dpn*(dwxi*xpn+dwyi*ypn+dwzi*zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++END: VELOCITIES AT CELL FACE CENTER and EXPLICIT DIFFUSION FLUXES+++++++

!.....DIFUSION COEFFICIENT
      DE=GAME*ARE/DPN

!.....CONVECTION FLUXES - UDS
!     If flow goes P=>E CP=FLCF, CE=0.
!     If flow goes E=>P CE=FLCF, CP=0.
      FLCF=FCF(INP)
      CE=MIN(FLCF,ZERO) 
      CP=MAX(FLCF,ZERO)
!
!.....COEFFICIENTS AE(P) AND AW(E) DUE TO UDS
!
      ACFE(INP)=DE-CE
      ACFW(INE)=DE+CP
!
!.....EXPLICIT CONVECTIVE FLUXES FOR UDS
!
      FUUDS=CP*U(INP)+CE*U(INE)
      FVUDS=CP*V(INP)+CE*V(INE)
      FWUDS=CP*W(INP)+CE*W(INE)

!=====START CDS SCHEME===============================
      IF(LCDS.EQ.1) THEN
!
!.....EXPLICIT CONVECTIVE FLUXES FOR CDS
!
      FUHIGH=FLCF*Ue
      FVHIGH=FLCF*Ve
      FWHIGH=FLCF*We


!=====END CDS SCHEME================================
      END IF


!=====START LUDS SCHEME=============================
      IF(LLUDS.EQ.1) THEN

!.....EXPLICIT CONVECTIVE FLUXES FOR LUDS SCHEME
!     $Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_e[kg/s] * Phi_e[m/s]$
!     Phi_e is found by extrapolation from upwind nodes, see eq. (3.29) in Sasa's Thesis.
!

      FUHIGH=CE*(U(INE)+(U(INE)-U(INEE))*FIF(INE))+ &
             CP*(U(INP)+(U(INP)-U(INW))*(1-FIF(INW)))
!     mass flux| interpolation of velocity to face |

      FVHIGH=CE*(V(INE)+(V(INE)-V(INEE))*FIF(INE))+ &
             CP*(V(INP)+(V(INP)-V(INW))*(1-FIF(INW)))
       
      FWHIGH=CE*(W(INE)+(W(INE)-W(INEE))*FIF(INE))+ &
             CP*(W(INP)+(W(INP)-W(INW))*(1-FIF(INW)))

!=====END LUDS SCHEME====================================
      END IF

!=====START QUICK SCHEME=================================
      IF(LQUDS.EQ.1) THEN


      IF(K.GT.2.AND.K.LT.NK-1.AND. &
         J.GT.2.AND.J.LT.NJ-1.AND. &
         I.GT.2.AND.I.LT.NI-1) THEN

!.....Find coefficients associated with quadratic function
      g11=((2-FIF(INW))*FIF(INP)**2)/  &
           (1+FIF(INP)-FIF(INW))
      g12=((1-FIF(INP))*(1-FIF(INW))**2)/  &
           (1+FIF(INP)-FIF(INW))
      g21=((1+FIF(INE))*(1-FIF(INP))**2)/  &
           (1+FIF(INE)-FIF(INP))
      g22=(FIF(INE)**2*FIF(INP))/ &
           (1+FIF(INE)-FIF(INP))

!.....EXPLICIT CONVECTIVE FLUXES FOR QUICK SCHEME
!     $Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_e[kg/s] * Phi_e[m/s]$
!     Phi_e is found by extrapolation from upwind nodes, see eq. (3.30) in Sasa's Thesis.

      FUHIGH=CE*(g21*U(INP)+(1-g21+g22)*U(INE)-g22*U(INEE))+ &
             CP*(g11*U(INE)+(1-g11+g12)*U(INP)-g12*U(INW))


      FVHIGH=CE*(g21*V(INP)+(1-g21+g22)*V(INE)-g22*V(INEE))+ &
             CP*(g11*V(INE)+(1-g11+g12)*V(INP)-g12*V(INW))
       
      FWHIGH=CE*(g21*W(INP)+(1-g21+g22)*W(INE)-g22*W(INEE))+ &
             CP*(g11*W(INE)+(1-g11+g12)*W(INP)-g12*W(INW))

      END IF

      END IF
!=============END QUICK SCHEME========================
!--------------------------------------------------------------------------------------------
!     BOUNDED HIGH-ORDER CONVECTIVE SCHEMES (Waterson & Deconinck JCP 224 (2007) pp. 182-207)
!     NOTES:
!     SMART- converges slow, requires stronger underrelaxation (like URF = 0.2), otherwise
!     probably smallest error when converges
!     AVL-SMART - better convergence than SMART, works fine with URF = 0.5
!     MUSCL - good convergence, URF = 0.6, even 0.8, the most robust one.
!     ...
!--------------------------------------------------------------------------------------------
      IF(LSMART.EQ.1.OR.LAVL.EQ.1.OR.LMUSCL.EQ.1.OR.LUMIST.EQ.1.OR.LGAMMA.EQ.1) THEN

! !+++++Find r's. This is universal for all schemes. NON-ORTHOGONAL GRID! +++++++
! !.....Distances of fictitious upwind point U from 'real' upwind cell center.
! !.....If flow goes from P to E, distance from the new node U from closest cell-center
!       DXU = (XC(INP) - (XC(INE)-XC(INP))) - XC(INW)
!       DYU = (YC(INP) - (YC(INE)-YC(INP))) - YC(INW)
!       DZU = (ZC(INP) - (ZC(INE)-ZC(INP))) - ZC(INW)
! !.....If flow goes from P to E, distance from the new node U from closest cell-center
!       DXU1 = (XC(INE) + (XC(INP)-XC(INW))) - XC(INEE)
!       DYU1 = (YC(INE) + (YC(INP)-YC(INW))) - YC(INEE)
!       DZU1 = (ZC(INE) + (ZC(INP)-ZC(INW))) - ZC(INEE)
! 
! !.....Variable value at fictitious upwind point
!       Uinu=U(INW)+DXU*GRADU(1,INW)+DYU*GRADU(2,INW)+DZU*GRADU(3,INW)
!       Vinu=V(INW)+DXU*GRADV(1,INW)+DYU*GRADV(2,INW)+DZU*GRADV(3,INW)
!       Winu=W(INW)+DXU*GRADW(1,INW)+DYU*GRADW(2,INW)+DZU*GRADW(3,INW)
! !
!       Uinu1=U(INEE)+DXU1*GRADU(1,INEE)+DYU1*GRADU(2,INEE)+DZU1*GRADU(3,INEE)
!       Vinu1=V(INEE)+DXU1*GRADV(1,INEE)+DYU1*GRADV(2,INEE)+DZU1*GRADV(3,INEE)
!       Winu1=W(INEE)+DXU1*GRADW(1,INEE)+DYU1*GRADW(2,INEE)+DZU1*GRADW(3,INEE)
! 
! !.....Gradient ratio using value at fictitious upwind point:
! !.....If flow goes from P to E
!       r1 = (U(INE)-U(INP))/(U(INP)-Uinu)
!       r2 = (V(INE)-V(INP))/(V(INP)-Vinu)
!       r3 = (W(INE)-W(INP))/(W(INP)-Winu)
! !.....If flow goes from E to P
!       r4 = (U(INP)-U(INE))/(U(INE)-Uinu1)
!       r5 = (V(INP)-V(INE))/(V(INE)-Vinu1)
!       r6 = (W(INP)-W(INE))/(W(INE)-Winu1)
      

!+++++Find r's. This is universal for all schemes. NON-UNIFORM GRID! ++++++
!.....If flow goes from P to E
!#      r1 = (U(INE)-U(INP))*(XC(INP)-XC(INW))/((U(INP)-U(INW))*(XC(INE)-XC(INP)))
!#      r2 = (V(INE)-V(INP))*(YC(INP)-YC(INW))/((V(INP)-V(INW))*(YC(INE)-YC(INP)))
!#      r3 = (W(INE)-W(INP))*(ZC(INP)-ZC(INW))/((W(INP)-W(INW))*(ZC(INE)-ZC(INP)))
!.....If flow goes from E to P
!#      r4 = (U(INP)-U(INE))*(XC(INE)-XC(INEE))/((U(INE)-U(INEE))*(XC(INP)-XC(INE)))
!#      r5 = (V(INP)-V(INE))*(YC(INE)-YC(INEE))/((V(INE)-V(INEE))*(YC(INP)-YC(INE)))
!#      r6 = (W(INP)-W(INE))*(ZC(INE)-ZC(INEE))/((W(INE)-W(INEE))*(ZC(INP)-ZC(INE)))

!+++++Find r's. This is universal for all schemes. UNIFORM GRID! ++++++
!.....If flow goes from P to E
!UO      r1 = (U(INE)-U(INP))/(U(INP)-U(INW))
!UO      r2 = (V(INE)-V(INP))/(V(INP)-V(INW))
!UO      r3 = (W(INE)-W(INP))/(W(INP)-W(INW))
!.....If flow goes from E to P
!UO      r4 = (U(INP)-U(INE))/(U(INE)-U(INEE))
!UO      r5 = (V(INP)-V(INE))/(V(INE)-V(INEE))
!UO      r6 = (W(INP)-W(INE))/(W(INE)-W(INEE))  

!+++++PATCH-VERSION 3++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!.... Find r's. This is universal for all schemes.
!.....If flow goes from P to E
      r1 = (2*gradU(1,inp)*xpn + 2*gradU(2,inp)*ypn + 2*gradU(3,inp)*zpn)/(U(INE)-U(INP)) - 1.0d0  
      r2 = (2*gradV(1,inp)*xpn + 2*gradV(2,inp)*ypn + 2*gradV(3,inp)*zpn)/(V(INE)-V(INP)) - 1.0d0 
      r3 = (2*gradW(1,inp)*xpn + 2*gradW(2,inp)*ypn + 2*gradW(3,inp)*zpn)/(W(INE)-W(INP)) - 1.0d0 
!.....If flow goes from E to P
      r4 = (2*gradU(1,ine)*xpn + 2*gradU(2,ine)*ypn + 2*gradU(3,ine)*zpn)/(U(INP)-U(INE)) - 1.0d0 
      r5 = (2*gradV(1,ine)*xpn + 2*gradV(2,ine)*ypn + 2*gradV(3,ine)*zpn)/(V(INP)-V(INE)) - 1.0d0 
      r6 = (2*gradW(1,ine)*xpn + 2*gradW(2,ine)*ypn + 2*gradW(3,ine)*zpn)/(W(INP)-W(INE)) - 1.0d0  
!+++++PATCH-VERSION 3++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!=====SMART SCHEME================================
      IF(LSMART.EQ.1) THEN
!.....PSI for SMART scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(2.*r1, 0.75*r1+0.25, 4.))
      PSIW2 = max(0., min(2.*r2, 0.75*r2+0.25, 4.))
      PSIW3 = max(0., min(2.*r3, 0.75*r3+0.25, 4.))
!.....If flow goes from E to P
      PSIE1 = max(0., min(2.*r4, 0.75*r4+0.25, 4.))
      PSIE2 = max(0., min(2.*r5, 0.75*r5+0.25, 4.))
      PSIE3 = max(0., min(2.*r6, 0.75*r6+0.25, 4.))
!=====END SMART SCHEME=============================
      END IF
!
!=====AVL-SMART SCHEME=============================
      IF(LAVL.EQ.1) THEN
!.....PSI for AVL-SMART scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(1.5*r1, 0.75*r1+0.25, 2.5))
      PSIW2 = max(0., min(1.5*r2, 0.75*r2+0.25, 2.5))
      PSIW3 = max(0., min(1.5*r3, 0.75*r3+0.25, 2.5))
!.....If flow goes from E to P
      PSIE1 = max(0., min(1.5*r4, 0.75*r4+0.25, 2.5))
      PSIE2 = max(0., min(1.5*r5, 0.75*r5+0.25, 2.5))
      PSIE3 = max(0., min(1.5*r6, 0.75*r6+0.25, 2.5))
!=====END AVL-SMART SCHEME==========================
      END IF
!
!=====MUSCL SCHEME=================================
      IF(LMUSCL.EQ.1) THEN
!.....PSI for MUSCL scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
      PSIW2 = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
      PSIW3 = max(0., min(2.*r3, 0.5*r3+0.5, 2.))
!.....If flow goes from E to P
      PSIE1 = max(0., min(2.*r4, 0.5*r4+0.5, 2.))
      PSIE2 = max(0., min(2.*r5, 0.5*r5+0.5, 2.))
      PSIE3 = max(0., min(2.*r6, 0.5*r6+0.5, 2.))

!.....PSI for Koren scheme:
!.....If flow goes from P to E
!      PSIW1 = max(0., min(2.*r1, twothirds*r1+onethird, 2.))
!      PSIW2 = max(0., min(2.*r2, twothirds*r2+onethird, 2.))
!      PSIW3 = max(0., min(2.*r3, twothirds*r3+onethird, 2.))
!.....If flow goes from E to P
!      PSIE1 = max(0., min(2.*r4, twothirds*r4+onethird, 2.))
!      PSIE2 = max(0., min(2.*r5, twothirds*r5+onethird, 2.))
!      PSIE3 = max(0., min(2.*r6, twothirds*r6+onethird, 2.))

!.....PSI for GPL-1/3-alpha-3/2 scheme:  !!!!>>>Koren i Ova schema su jako slicne
!.....If flow goes from P to E
!      PSIW1 = max(0., min(1.5*r1, twothirds*r1+onethird, 2.))
!      PSIW2 = max(0., min(1.5*r2, twothirds*r2+onethird, 2.))
!      PSIW3 = max(0., min(1.5*r3, twothirds*r3+onethird, 2.))
!.....If flow goes from E to P
!      PSIE1 = max(0., min(1.5*r4, twothirds*r4+onethird, 2.))
!      PSIE2 = max(0., min(1.5*r5, twothirds*r5+onethird, 2.))
!      PSIE3 = max(0., min(1.5*r6, twothirds*r6+onethird, 2.))

!.....PSI for SMARTER; CHARM NOTABLE; ISNAS
!.....If flow goes from P to E
!      PSIW1 = (r1+abs(r1))*(3*r1+1.)/(2*(r1+1.)**2)
!      PSIW2 = (r2+abs(r2))*(3*r2+1.)/(2*(r2+1.)**2)
!      PSIW3 = (r3+abs(r3))*(3*r3+1.)/(2*(r3+1.)**2)
!.....If flow goes from E to P
!      PSIE1 = (r4+abs(r4))*(3*r4+1.)/(2*(r4+1.)**2)
!      PSIE2 = (r5+abs(r5))*(3*r5+1.)/(2*(r5+1.)**2)
!      PSIE3 = (r6+abs(r6))*(3*r6+1.)/(2*(r6+1.)**2)

!.....PSI for OSPRE
!.....If flow goes from P to E
!      PSIW1 = 1.5*r1*(r1+1.)/(r1**2+r1+1)
!      PSIW2 = 1.5*r2*(r2+1.)/(r2**2+r2+1)
!      PSIW3 = 1.5*r3*(r3+1.)/(r3**2+r3+1)
!.....If flow goes from E to P
!      PSIE1 = 1.5*r4*(r4+1.)/(r4**2+r4+1)
!      PSIE2 = 1.5*r5*(r5+1.)/(r5**2+r5+1)
!      PSIE3 = 1.5*r6*(r6+1.)/(r6**2+r6+1)

!.....PSI for BSOU-BLUI-Chakravarthy-Osher scheme:
!.....If flow goes from P to E
!      PSIW1 = max(0., min(2.*r1,1.))
!      PSIW2 = max(0., min(2.*r2,1.))
!      PSIW3 = max(0., min(2.*r3,1.))
!.....If flow goes from E to P
!      PSIE1 = max(0., min(2.*r4,1.))
!      PSIE2 = max(0., min(2.*r5,1.))
!      PSIE3 = max(0., min(2.*r6,1.))
!=====END MUSCL SCHEME=============================
      END IF
!=====UMIST SCHEME=================================
      IF(LUMIST.EQ.1) THEN
!.....PSI for UMIST scheme:
!.....If flow goes from P to E
!      PSIW1 = max(0., min(2.*r1, 0.75*r1+0.25, 0.25*r1+0.75, 2.))
!      PSIW2 = max(0., min(2.*r2, 0.75*r2+0.25, 0.25*r2+0.75, 2.))
!      PSIW3 = max(0., min(2.*r3, 0.75*r3+0.25, 0.25*r3+0.75, 2.))
!.....If flow goes from E to P
!      PSIE1 = max(0., min(2.*r4, 0.75*r4+0.25, 0.25*r4+0.75, 2.))
!      PSIE2 = max(0., min(2.*r5, 0.75*r5+0.25, 0.25*r5+0.75, 2.))
!      PSIE3 = max(0., min(2.*r6, 0.75*r6+0.25, 0.25*r6+0.75, 2.))
!=====END UMIST SCHEME=============================
      ENDIF
!=====SPL-3/5 SCHEME===============================
     IF(LUMIST.EQ.1) THEN
!.....PSI for SPL-3/5 scheme (NEW SCHEME derived from Waterson&Deconinck's Symmetric piecewise-linear scheme):
!.....If flow goes from P to E
      PSIW1 = max(0., min(2.*r1, 0.8*r1+0.2, 0.2*r1+0.8, 2.))
      PSIW2 = max(0., min(2.*r2, 0.8*r2+0.2, 0.2*r2+0.8, 2.))
      PSIW3 = max(0., min(2.*r3, 0.8*r3+0.2, 0.2*r3+0.8, 2.))
!.....If flow goes from E to P
      PSIE1 = max(0., min(2.*r4, 0.8*r4+0.2, 0.2*r4+0.8, 2.))
      PSIE2 = max(0., min(2.*r5, 0.8*r5+0.2, 0.2*r5+0.8, 2.))
      PSIE3 = max(0., min(2.*r6, 0.8*r6+0.2, 0.2*r6+0.8, 2.))
!=====END UMIST SCHEME=============================
      ENDIF
!=====GAMMA SCHEME================================
      IF(LGAMMA.EQ.1) THEN
!.....PSI for GAMMA scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(r1, 2.*r1/(r1+1.)))
      PSIW2 = max(0., min(r2, 2.*r2/(r2+1.)))
      PSIW3 = max(0., min(r3, 2.*r3/(r3+1.)))
!.....If flow goes from E to P
      PSIE1 = max(0., min(r4, 2.*r4/(r4+1.)))
      PSIE2 = max(0., min(r5, 2.*r5/(r5+1.)))
      PSIE3 = max(0., min(r6, 2.*r6/(r6+1.)))
!=====END GAMMA SCHEME=============================
      END IF

!.....EXPLICIT CONVECTIVE FLUXES FOR HIGH ORDER BOUNDED SCHEMES
!     $Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_e[kg/s] * Phi_e[m/s]$
!     Phi_e is found by extrapolation from upwind nodes, see eq. (3.29) in Sasa's Thesis.
!     Additional multiplication with PSI is application of flux limiters,
!     see eq. (10) in Waterson&Deconinck paper.

!@      fif(ine)=0.5d0  ! A hack for uniform grid 
!@      fif(inw)=0.5d0  !

!      FUHIGH=CE*(U(INE)*(1+FIF(INE)*PSIE1)-Uinu1*FIF(INE)*PSIE1)+ &
!       CP*(U(INP)*(1+(1-FIF(INW))*PSIW1)-Uinu*(1-FIF(INW))*PSIW1)

!      FVHIGH=CE*(V(INE)*(1+FIF(INE)*PSIE2)-Vinu1*FIF(INE)*PSIE2)+ &
!       CP*(V(INP)*(1+(1-FIF(INW))*PSIW2)-Vinu*(1-FIF(INW))*PSIW2)
       
!      FWHIGH=CE*(W(INE)*(1+FIF(INE)*PSIE3)-Winu1*FIF(INE)*PSIE3)+ &
!       CP*(W(INP)*(1+(1-FIF(INW))*PSIW3)-Winu*(1-FIF(INW))*PSIW3)

!.....Ver.3
!......Darwish-Moukalled TVD schemes for unstructured girds, IJHMT, 2003.
      FUHIGH = CE*(U(INE) + FXE*PSIE1*(U(INP)-U(INE)))+ &
               CP*(U(INP) + FXP*PSIW1*(U(INE)-U(INP)))
!       mass flux| bounded interpolation of velocity to face |

      FVHIGH = CE*(V(INE) + FXE*PSIE2*(V(INP)-V(INE)))+ &
               CP*(V(INP) + FXP*PSIW2*(V(INE)-V(INP)))
       
      FWHIGH = CE*(W(INE) + FXE*PSIE3*(W(INP)-W(INE)))+ &
               CP*(W(INP) + FXP*PSIW3*(W(INE)-W(INP)))
!......END: Darwish-Moukalled TVDschemes for unstructured girds, IJHMT, 2003.

!.....END OF BOUNDED HIGH-ORDER SCHEMES
      END IF 

!
!.....EXPLICIT PART OF DIFFUSION FLUXES AND SOURCES DUE TO DEFFERED CORRECTION
!     FOR ALL SCHEMES!
!
      SU(INP)=SU(INP)+GAM*(FUUDS-FUHIGH)+fdue-fdui
      SV(INP)=SV(INP)+GAM*(FVUDS-FVHIGH)+fdve-fdvi
      SW(INP)=SW(INP)+GAM*(FWUDS-FWHIGH)+fdwe-fdwi
!----------------------------------------------------
!......[Vectorization procedure: ]
!----------------------------------------------------
      BP(INE)=-GAM*(FUUDS-FUHIGH)-fdue+fdui
      BT(INE)=-GAM*(FVUDS-FVHIGH)-fdve+fdvi
      BB(INE)=-GAM*(FWUDS-FWHIGH)-fdwe+fdwi
!----------------------------------------------------
      END DO !J-loop
      END DO !I-loop
      END DO !K-loop

      DO K=2,NKM
      DO I=2,NIM
      DO J=2,NJM
      INP=LK(K)+LI(I)+J+IDEW

      SU(INP)=SU(INP)+BP(INP)
      SV(INP)=SV(INP)+BT(INP)
      SW(INP)=SW(INP)+BB(INP)

      END DO !J-loop
      END DO !I-loop
      END DO !K-loop

      DO IJK=ICST,ICEN
      BP(IJK)=0.0D0
      BT(IJK)=0.0D0
      BB(IJK)=0.0D0
      END DO

      RETURN
      END
