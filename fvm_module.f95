!
!
!
! Polymorphism in Fortran is achieved using a generic interface block.
!
! ddt_LLrho_ULL__plus__div_LLphi_ULL__minus__laplacian_LLmu_ULL 
! ddt_LLrho_FiLL__plus__div_LLphi_FiLL__minus__laplacian_LLGame_FiLL 

module finite_volume_method

   implicit none

   interface fvm_div_neglap
      module procedure fvm_div_neglap_velocity
      module procedure fvm_div_neglap_scalar
   end interface fvm_div_neglap

   private  ! hides items not listed on public statement 
   public :: fvm_div_neglap

contains


      SUBROUTINE fvm_div_neglap_velocity(phi,mu,U,V,W)
!  
!******************************************************************************
!
!     Fills matrix with coefficients representing implicit FVM discretization 
!     of Convective term: div(phi*U).
!     Where phi is mass flux surface field, and U is velocity field.
!
!     Inteded for structured 3D grids. Non-orthogonal corrections are included.
!
!     Matrix is stored in diagonal format with seven diagonals.
!     Each diagonal is stored in an array,
!     Off-main diagonnal: ae,aw,an,as,at,ab,
!     Main diagonal: ap.
!     RHS vector is SU.
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!
      Use Types
      Use Parameters
      Use Indexes
      Use Geometry
      Use Coef
      Use Coefb
      Use Variables, only: U,V,W,F1,F2,F3,Vis
      Use Gradients, only: gradU,gradV,gradW

      Implicit None

      !Real(Prec), Dimension(Nxyza), Intent(InOut) :: Phi ! <- mass flux F1,F2,F3 or generic for unstructured mesh.
!
!     Local Variables
!
      Integer :: I, J, K, Inp
      Integer :: Ine, Ins, Inb, Inbs, Idew, Idns, Idtb


!@      Integer :: Indx
      Real(Prec) :: Dxc,Dyc,Dzc,Dxe,Dye,Dze
      Real(Prec) :: Dxu,Dyu,Dzu,Dxu1,Dyu1,Dzu1,Uinu,Vinu,Winu,Uinu1,Vinu1,Winu1
      Real(Prec) :: Dfidx_U,Dfidx_U1, &
                    Dfidy_U,Dfidy_U1, &
                    Dfidz_U,Dfidz_U1
!      Real(Prec) :: Dfidx_F,Dfidy_F,Dfidz_F
!      Real(Prec) :: Dxet,Dxzd,Dyet,Dyzd,Dzet,Dzzd
      Real(Prec) :: Gam,Fxe,Fxp,Are,Vole,Game,&
      Real(Prec) :: Ue,Ve,We
      Real(Prec) :: De,Flcf,Ce,Cp
      Real(Prec) :: G11, G12, G21, G22
      Real(Prec) :: Ae1,Aw1,Shigh1,Shigh2,Shigh3
      Real(Prec) :: R1,R2,R3,R4,R5,R6                   
      Real(Prec) :: Psie1,Psie2,Psie3,Psiw1,Psiw2,Psiw3
      Real(Prec) :: Fuuds,Fvuds,Fwuds,Fuhigh,Fvhigh,Fwhigh
      Real(Prec) :: Xf,Yf,Zf,Xi,Yi,Zi,Arx,Ary,Arz
      Real(Prec) :: Duxi,Duyi,Duzi, &
                    Dvxi,Dvyi,Dvzi, &
                    Dwxi,Dwyi,Dwzi
      Real(Prec) :: Duxii,Dvxii,Dwxii, &
                    Duyii,Dvyii,Dwyii, &
                    Duzii,Dvzii,Dwzii
      Real(Prec) :: Xpn,Ypn,Zpn
      Real(Prec) :: nxx,nyy,nzz,Ixi1,Ixi2,Ixi3,Dpn,Costheta,Costn
      Real(Prec) :: Fdue,Fdve,Fdwe,Fdui,Fdvi,Fdwi
      Real(Prec) :: D2x,D2y,D2z,D1x,D1y,D1z

!----------------------------------------------------------------------
      Gam=Gds(Iu) ! Deffered correction blending coefficient

!.....Initialize Values For High Order Conv. Fluxes
      Psie1=0.0d0;Psie2=0.0d0;Psie3=0.0d0
      Psiw1=0.0d0;Psiw2=0.0d0;Psiw3=0.0d0   
      
      Fuhigh=0.0d0;Fvhigh=0.0d0;Fwhigh=0.0d0

!
!.....Calculate East Cell Face
      idew = nj
      idns = 1
      idtb = nij

      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J
      ine=inp+idew
      inw=inp-idew
      inn=inp+idns
      ins=inp-idns
      inb=inp-idtb
      int=inp+idtb
      inbs=inb-idns
!
!.....Interpolation Factor
      fxe=fx(inp) 
      fxp=1.-fxe

!.....Components Of Distance (From P To Neighbor N) Vector
!.....First
      Xpn=Xc(Ine)-Xc(Inp)
      Ypn=Yc(Ine)-Yc(Inp)
      Zpn=Zc(Ine)-Zc(Inp)

!.....Distance From P To Neighbor N
      Dpn=Sqrt(Xpn**2+Ypn**2+Zpn**2)  
   
!.....Components Of The Unit Vector I_Ksi
      Ixi1=Xpn/Dpn
      Ixi2=Ypn/Dpn
      Ixi3=Zpn/Dpn

!.....Precomputed Face Areas
      Arx=Ar1x(Inp)
      Ary=Ar1y(Inp)
      Arz=Ar1z(Inp)

!.....Cell Face Area
      Are=Sqrt(Arx**2+Ary**2+Arz**2)

!.....Unit Vectors Of The Normal
      nxx=Arx/Are
      nyy=Ary/Are
      nzz=Arz/Are

!.....Angle Between Vectors n And i_xi - We Need Cosine
      costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

!.....Relaxation Factor For Higher-Order Cell Face Gradient
      ! Minimal Correction: Nrelax = +1 :
      !Costn = Costheta
      ! Orthogonal Correction: Nrelax =  0 : 
      Costn = 1.0d0
      ! Over-Relaxed Approach: Nrelax = -1 :
      !Costn = 1./Costheta
      ! In General, Nrelax Can Be Any Signed Integer From Some 
      ! Reasonable Interval [-Nrelax,Nrelax] (Or Maybe Even Real Number): 
      !Costn = Costheta**Nrelax

!.....d_pn . S_f
      Vole=Xpn*Arx+Ypn*Ary+Zpn*Arz

!.....Overrelaxed Correction Vector D2, Where S=Dpn+D2
      D1x = Costn
      D1y = Costn
      D1z = Costn
      !
      D2x = Xpn*Costn
      D2y = Ypn*Costn
      D2z = Zpn*Costn

!.....Cell Face Viscosity
      Game=(Vis(Inp)*Fxp+Vis(Ine)*Fxe)

!++++Velocities At Cell Face Center And Explicit Diffusion Fluxes+++++++

!.....Coordinates Of Cell-Face Center
      Xf = 0.25*(X(Inp)+X(Ins)+X(Inb)+X(Inbs))
      Yf = 0.25*(Y(Inp)+Y(Ins)+Y(Inb)+Y(Inbs))
      Zf = 0.25*(Z(Inp)+Z(Ins)+Z(Inb)+Z(Inbs))

!.....Coordinates Of Point E'
      Xi=Xc(Inp)*Fxp+Xc(Ine)*Fxe
      Yi=Yc(Inp)*Fxp+Yc(Ine)*Fxe
      Zi=Zc(Inp)*Fxp+Zc(Ine)*Fxe

!.....Interpolate Gradients Defined At Cv Centers To Faces
      Duxi = Gradu(1,Inp)*Fxp+Gradu(1,Ine)*Fxe
      Duyi = Gradu(2,Inp)*Fxp+Gradu(2,Ine)*Fxe
      Duzi = Gradu(3,Inp)*Fxp+Gradu(3,Ine)*Fxe

!        |________Ue'_________|_______________Ucorr___________________|
      Ue=U(Inp)*Fxp+U(Ine)*Fxe+(Duxi*(Xf-Xi)+Duyi*(Yf-Yi)+Duzi*(Zf-Zi))

!.....Du/Dx_I Interpolated At Cell Face:
      Duxii = Duxi*D1x + Arx/Vole*( U(Ine)-U(Inp)-Duxi*D2x-Duyi*D2y-Duzi*D2z ) 
      Duyii = Duyi*D1y + Ary/Vole*( U(Ine)-U(Inp)-Duxi*D2x-Duyi*D2y-Duzi*D2z ) 
      Duzii = Duzi*D1z + Arz/Vole*( U(Ine)-U(Inp)-Duxi*D2x-Duyi*D2y-Duzi*D2z ) 

      Dvxi = Gradv(1,Inp)*Fxp+Gradv(1,Ine)*Fxe
      Dvyi = Gradv(2,Inp)*Fxp+Gradv(2,Ine)*Fxe
      Dvzi = Gradv(3,Inp)*Fxp+Gradv(3,Ine)*Fxe

!        |________Ve'_________|_______________Vcorr___________________|
      Ve=V(Inp)*Fxp+V(Ine)*Fxe+(Dvxi*(Xf-Xi)+Dvyi*(Yf-Yi)+Dvzi*(Zf-Zi))

!.....Dv/Dx_I Interpolated At Cell Face:
      Dvxii = Dvxi*D1x + Arx/Vole*( V(Ine)-V(Inp)-Dvxi*D2x-Dvyi*D2y-Dvzi*D2z ) 
      Dvyii = Dvyi*D1y + Ary/Vole*( V(Ine)-V(Inp)-Dvxi*D2x-Dvyi*D2y-Dvzi*D2z ) 
      Dvzii = Dvzi*D1z + Arz/Vole*( V(Ine)-V(Inp)-Dvxi*D2x-Dvyi*D2y-Dvzi*D2z ) 

      Dwxi = Gradw(1,Inp)*Fxp+Gradw(1,Ine)*Fxe
      Dwyi = Gradw(2,Inp)*Fxp+Gradw(2,Ine)*Fxe
      Dwzi = Gradw(3,Inp)*Fxp+Gradw(3,Ine)*Fxe

!        |________We'_________|_______________Wcorr___________________|
      We=W(Inp)*Fxp+W(Ine)*Fxe+(Dwxi*(Xf-Xi)+Dwyi*(Yf-Yi)+Dwzi*(Zf-Zi))

!.....Dw/Dx_I Interpolated At Cell Face:
      Dwxii = Dwxi*D1x + Arx/Vole*( W(Ine)-W(Inp)-Dwxi*D2x-Dwyi*D2y-Dwzi*D2z ) 
      Dwyii = Dwyi*D1y + Ary/Vole*( W(Ine)-W(Inp)-Dwxi*D2x-Dwyi*D2y-Dwzi*D2z ) 
      Dwzii = Dwzi*D1z + Arz/Vole*( W(Ine)-W(Inp)-Dwxi*D2x-Dwyi*D2y-Dwzi*D2z ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     We calculate explicit and implicit diffusion fde and fdi,
!     later we put their difference (fde-fdi) to rhs vector.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ! Explicit Diffussion 
       Fdue = Game*(Duxii+Duxii)*Arx + (Duyii+Dvxii)*Ary + (Duzii+Dwxii)*Arz
       Fdve = Game*(Duyii+Dvxii)*Arx + (Dvyii+Dvyii)*Ary + (Dvzii+Dwyii)*Arz
       Fdwe = Game*(Duzii+Dwxii)*Arx + (Dwyii+Dvzii)*Ary + (Dwzii+Dwzii)*Arz
       
       ! Implicit Diffussion 
       Fdui = Game*Are/Dpn*(Duxi*Xpn+Duyi*Ypn+Duzi*Zpn)
       Fdvi = Game*Are/Dpn*(Dvxi*Xpn+Dvyi*Ypn+Dvzi*Zpn)
       Fdwi = Game*Are/Dpn*(Dwxi*Xpn+Dwyi*Ypn+Dwzi*Zpn)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++End: Velocities At Cell Face Center And Explicit Diffusion Fluxes+++++++


!.....Difusion Coefficient
      De=Game*Are/Dpn

!.....Convection Fluxes - Uds
!     If Flow Goes P=>E Cp=Flcf, Ce=0.
!     If Flow Goes E=>P Ce=Flcf, Cp=0.
      Flcf = F1(Inp)
      Ce = Min(Flcf,Zero) 
      Cp = Max(Flcf,Zero)
!
!.....Coefficients Ae(P) And Aw(E) Due To Uds
!
      Ae(Inp) = De-Ce
      Aw(Ine) = De+Cp

!
!.....EXPLICIT CONVECTIVE FLUXES FOR UDS
!
      Fuuds=Cp*U(Inp)+Ce*U(Ine)
      Fvuds=Cp*V(Inp)+Ce*V(Ine)
      Fwuds=Cp*W(Inp)+Ce*W(Ine)

!=====START CDS SCHEME===============================
      IF(LCDS.EQ.1) THEN
!
!.....EXPLICIT CONVECTIVE FLUXES FOR CDS
!
      Fuhigh=Flcf*Ue
      Fvhigh=Flcf*Ve
      Fwhigh=Flcf*We

!=====END CDS SCHEME================================
      END IF


!=====START MUSCL SCHEME===============================
      IF(LMUSCL.EQ.1) THEN
!
!.....EXPLICIT CONVECTIVE FLUXES FOR MUSCL
!
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
!.....PSI for MUSCL scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
      PSIW2 = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
      PSIW3 = max(0., min(2.*r3, 0.5*r3+0.5, 2.))
!.....If flow goes from E to P
      PSIE1 = max(0., min(2.*r4, 0.5*r4+0.5, 2.))
      PSIE2 = max(0., min(2.*r5, 0.5*r5+0.5, 2.))
      PSIE3 = max(0., min(2.*r6, 0.5*r6+0.5, 2.))

!.....Ver.3
!......Darwish-Moukalled TVD schemes for unstructured girds, IJHMT, 2003.
      Fuhigh = Ce*(U(Ine) + Fxe*Psie1*(U(Inp)-U(Ine)))+ &
               Cp*(U(Inp) + Fxp*Psiw1*(U(Ine)-U(Inp)))
!       mass flux| bounded interpolation of velocity to face |

      Fvhigh = Ce*(V(Ine) + Fxe*Psie2*(W(Inp)-W(Ine)))+ &
               Cp*(V(Inp) + Fxp*Psiw2*(V(Ine)-V(Inp)))
       
      Fwhigh = Ce*(W(Ine) + Fxe*Psie3*(W(Inp)-W(Ine)))+ &
               Cp*(W(Inp) + Fxp*Psiw3*(W(Ine)-W(Inp)))
!......END: Darwish-Moukalled TVDschemes for unstructured girds, IJHMT, 2003.

!.....END OF BOUNDED HIGH-ORDER SCHEMES
      END IF 

!
!.....EXPLICIT PART AND SOURCES DUE TO DEFFERED CORRECTION
!
      Su(Inp) = Su(Inp)+Gam*(Fuuds-Fuhigh)+Fdue-Fdui
      Sv(Inp) = Sv(Inp)+Gam*(Fvuds-Fvhigh)+Fdve-Fdvi
      Sw(Inp) = Sw(Inp)+Gam*(Fwuds-Fwhigh)+Fdwe-Fdwi
!----------------------------------------------------
!......[Vectorization procedure: ]
!----------------------------------------------------
      Bp(Ine) = -Gam*(Fuuds-Fuhigh)-Fdue+Fdui
      Bt(Ine) = -Gam*(Fvuds-Fvhigh)-Fdve+Fdvi
      Bb(Ine) = -Gam*(Fwuds-Fwhigh)-Fdwe+Fdwi

      END DO
      END DO
      END DO

      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J+Idew

      Su(Inp)=Su(Inp)+Bp(Inp)
      Sv(Inp)=Sv(Inp)+Bt(Inp)
      Sw(Inp)=Sw(Inp)+Bb(Inp)

      End Do !J-Loop
      End Do !I-Loop
      End Do !K-Loop

      do ijk=icst,icen
      bp(ijk)=0.0d0
      bt(ijk)=0.0d0
      bb(ijk)=0.0d0
      end do

!/////Stop here for unstructured grid. we weill have only one face loop!
!/////But for the present structured grid we have to make north and top faces loop...

!
!.....Calculate North Cell Face
      idew = 1
      idns = nij
      idtb = nj

      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J
      ine=inp+idew
      inw=inp-idew
      inn=inp+idns
      ins=inp-idns
      inb=inp-idtb
      int=inp+idtb
      inbs=inb-idns
!
!.....Interpolation Factor
      fxe=fy(inp) 
      fxp=1.-fxe

!.....Components Of Distance (From P To Neighbor N) Vector
!.....First
      Xpn=Xc(Ine)-Xc(Inp)
      Ypn=Yc(Ine)-Yc(Inp)
      Zpn=Zc(Ine)-Zc(Inp)

!.....Distance From P To Neighbor N
      Dpn=Sqrt(Xpn**2+Ypn**2+Zpn**2)  
   
!.....Components Of The Unit Vector I_Ksi
      Ixi1=Xpn/Dpn
      Ixi2=Ypn/Dpn
      Ixi3=Zpn/Dpn

!.....Precomputed Face Areas
      Arx=Ar2x(Inp)
      Ary=Ar2y(Inp)
      Arz=Ar2z(Inp)

!.....Cell Face Area
      Are=Sqrt(Arx**2+Ary**2+Arz**2)

!.....Unit Vectors Of The Normal
      nxx=Arx/Are
      nyy=Ary/Are
      nzz=Arz/Are

!.....Angle Between Vectors n And i_xi - We Need Cosine
      costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

!.....Relaxation Factor For Higher-Order Cell Face Gradient
      ! Minimal Correction: Nrelax = +1 :
      !Costn = Costheta
      ! Orthogonal Correction: Nrelax =  0 : 
      Costn = 1.0d0
      ! Over-Relaxed Approach: Nrelax = -1 :
      !Costn = 1./Costheta
      ! In General, Nrelax Can Be Any Signed Integer From Some 
      ! Reasonable Interval [-Nrelax,Nrelax] (Or Maybe Even Real Number): 
      !Costn = Costheta**Nrelax

!.....d_pn . S_f
      Vole=Xpn*Arx+Ypn*Ary+Zpn*Arz

!.....Overrelaxed Correction Vector D2, Where S=Dpn+D2
      D1x = Costn
      D1y = Costn
      D1z = Costn
      !
      D2x = Xpn*Costn
      D2y = Ypn*Costn
      D2z = Zpn*Costn

!.....Cell Face Viscosity
      Game=(Vis(Inp)*Fxp+Vis(Ine)*Fxe)

!++++Velocities At Cell Face Center And Explicit Diffusion Fluxes+++++++

!.....Coordinates Of Cell-Face Center
      Xf = 0.25*(X(Inp)+X(Ins)+X(Inb)+X(Inbs))
      Yf = 0.25*(Y(Inp)+Y(Ins)+Y(Inb)+Y(Inbs))
      Zf = 0.25*(Z(Inp)+Z(Ins)+Z(Inb)+Z(Inbs))

!.....Coordinates Of Point E'
      Xi=Xc(Inp)*Fxp+Xc(Ine)*Fxe
      Yi=Yc(Inp)*Fxp+Yc(Ine)*Fxe
      Zi=Zc(Inp)*Fxp+Zc(Ine)*Fxe

!.....Interpolate Gradients Defined At Cv Centers To Faces
      Duxi = Gradu(1,Inp)*Fxp+Gradu(1,Ine)*Fxe
      Duyi = Gradu(2,Inp)*Fxp+Gradu(2,Ine)*Fxe
      Duzi = Gradu(3,Inp)*Fxp+Gradu(3,Ine)*Fxe

!        |________Ue'_________|_______________Ucorr___________________|
      Ue=U(Inp)*Fxp+U(Ine)*Fxe+(Duxi*(Xf-Xi)+Duyi*(Yf-Yi)+Duzi*(Zf-Zi))

!.....Du/Dx_I Interpolated At Cell Face:
      Duxii = Duxi*D1x + Arx/Vole*( U(Ine)-U(Inp)-Duxi*D2x-Duyi*D2y-Duzi*D2z ) 
      Duyii = Duyi*D1y + Ary/Vole*( U(Ine)-U(Inp)-Duxi*D2x-Duyi*D2y-Duzi*D2z ) 
      Duzii = Duzi*D1z + Arz/Vole*( U(Ine)-U(Inp)-Duxi*D2x-Duyi*D2y-Duzi*D2z ) 

      Dvxi = Gradv(1,Inp)*Fxp+Gradv(1,Ine)*Fxe
      Dvyi = Gradv(2,Inp)*Fxp+Gradv(2,Ine)*Fxe
      Dvzi = Gradv(3,Inp)*Fxp+Gradv(3,Ine)*Fxe

!        |________Ve'_________|_______________Vcorr___________________|
      Ve=V(Inp)*Fxp+V(Ine)*Fxe+(Dvxi*(Xf-Xi)+Dvyi*(Yf-Yi)+Dvzi*(Zf-Zi))

!.....Dv/Dx_I Interpolated At Cell Face:
      Dvxii = Dvxi*D1x + Arx/Vole*( V(Ine)-V(Inp)-Dvxi*D2x-Dvyi*D2y-Dvzi*D2z ) 
      Dvyii = Dvyi*D1y + Ary/Vole*( V(Ine)-V(Inp)-Dvxi*D2x-Dvyi*D2y-Dvzi*D2z ) 
      Dvzii = Dvzi*D1z + Arz/Vole*( V(Ine)-V(Inp)-Dvxi*D2x-Dvyi*D2y-Dvzi*D2z ) 

      Dwxi = Gradw(1,Inp)*Fxp+Gradw(1,Ine)*Fxe
      Dwyi = Gradw(2,Inp)*Fxp+Gradw(2,Ine)*Fxe
      Dwzi = Gradw(3,Inp)*Fxp+Gradw(3,Ine)*Fxe

!        |________We'_________|_______________Wcorr___________________|
      We=W(Inp)*Fxp+W(Ine)*Fxe+(Dwxi*(Xf-Xi)+Dwyi*(Yf-Yi)+Dwzi*(Zf-Zi))

!.....Dw/Dx_I Interpolated At Cell Face:
      Dwxii = Dwxi*D1x + Arx/Vole*( W(Ine)-W(Inp)-Dwxi*D2x-Dwyi*D2y-Dwzi*D2z ) 
      Dwyii = Dwyi*D1y + Ary/Vole*( W(Ine)-W(Inp)-Dwxi*D2x-Dwyi*D2y-Dwzi*D2z ) 
      Dwzii = Dwzi*D1z + Arz/Vole*( W(Ine)-W(Inp)-Dwxi*D2x-Dwyi*D2y-Dwzi*D2z ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     We calculate explicit and implicit diffusion fde and fdi,
!     later we put their difference (fde-fdi) to rhs vector.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ! Explicit Diffussion 
       Fdue = Game*(Duxii+Duxii)*Arx + (Duyii+Dvxii)*Ary + (Duzii+Dwxii)*Arz
       Fdve = Game*(Duyii+Dvxii)*Arx + (Dvyii+Dvyii)*Ary + (Dvzii+Dwyii)*Arz
       Fdwe = Game*(Duzii+Dwxii)*Arx + (Dwyii+Dvzii)*Ary + (Dwzii+Dwzii)*Arz
       
       ! Implicit Diffussion 
       Fdui = Game*Are/Dpn*(Duxi*Xpn+Duyi*Ypn+Duzi*Zpn)
       Fdvi = Game*Are/Dpn*(Dvxi*Xpn+Dvyi*Ypn+Dvzi*Zpn)
       Fdwi = Game*Are/Dpn*(Dwxi*Xpn+Dwyi*Ypn+Dwzi*Zpn)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++End: Velocities At Cell Face Center And Explicit Diffusion Fluxes+++++++


!.....Difusion Coefficient
      De=Game*Are/Dpn

!.....Convection Fluxes - Uds
!     If Flow Goes P=>E Cp=Flcf, Ce=0.
!     If Flow Goes E=>P Ce=Flcf, Cp=0.
      Flcf = F2(Inp)
      Ce = Min(Flcf,Zero) 
      Cp = Max(Flcf,Zero)
!
!.....Coefficients Ae(P) And Aw(E) Due To Uds
!
      An(Inp) = De-Ce
      As(Ine) = De+cp

!
!.....EXPLICIT CONVECTIVE FLUXES FOR UDS
!
      Fuuds=Cp*U(Inp)+Ce*U(Ine)
      Fvuds=Cp*V(Inp)+Ce*V(Ine)
      Fwuds=Cp*W(Inp)+Ce*W(Ine)

!=====START CDS SCHEME===============================
      IF(LCDS.EQ.1) THEN
!
!.....EXPLICIT CONVECTIVE FLUXES FOR CDS
!
      Fuhigh=Flcf*Ue
      Fvhigh=Flcf*Ve
      Fwhigh=Flcf*We

!=====END CDS SCHEME================================
      END IF


!=====START MUSCL SCHEME===============================
      IF(LMUSCL.EQ.1) THEN
!
!.....EXPLICIT CONVECTIVE FLUXES FOR MUSCL
!
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
!.....PSI for MUSCL scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
      PSIW2 = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
      PSIW3 = max(0., min(2.*r3, 0.5*r3+0.5, 2.))
!.....If flow goes from E to P
      PSIE1 = max(0., min(2.*r4, 0.5*r4+0.5, 2.))
      PSIE2 = max(0., min(2.*r5, 0.5*r5+0.5, 2.))
      PSIE3 = max(0., min(2.*r6, 0.5*r6+0.5, 2.))

!.....Ver.3
!......Darwish-Moukalled TVD schemes for unstructured girds, IJHMT, 2003.
      Fuhigh = Ce*(U(Ine) + Fxe*Psie1*(U(Inp)-U(Ine)))+ &
               Cp*(U(Inp) + Fxp*Psiw1*(U(Ine)-U(Inp)))
!       mass flux| bounded interpolation of velocity to face |

      Fvhigh = Ce*(V(Ine) + Fxe*Psie2*(W(Inp)-W(Ine)))+ &
               Cp*(V(Inp) + Fxp*Psiw2*(V(Ine)-V(Inp)))
       
      Fwhigh = Ce*(W(Ine) + Fxe*Psie3*(W(Inp)-W(Ine)))+ &
               Cp*(W(Inp) + Fxp*Psiw3*(W(Ine)-W(Inp)))
!......END: Darwish-Moukalled TVDschemes for unstructured girds, IJHMT, 2003.

!.....END OF BOUNDED HIGH-ORDER SCHEMES
      END IF 

!
!.....EXPLICIT PART AND SOURCES DUE TO DEFFERED CORRECTION
!
      Su(Inp) = Su(Inp)+Gam*(Fuuds-Fuhigh)+Fdue-Fdui
      Sv(Inp) = Sv(Inp)+Gam*(Fvuds-Fvhigh)+Fdve-Fdvi
      Sw(Inp) = Sw(Inp)+Gam*(Fwuds-Fwhigh)+Fdwe-Fdwi
!----------------------------------------------------
!......[Vectorization procedure: ]
!----------------------------------------------------
      Bp(Ine) = -Gam*(Fuuds-Fuhigh)-Fdue+Fdui
      Bt(Ine) = -Gam*(Fvuds-Fvhigh)-Fdve+Fdvi
      Bb(Ine) = -Gam*(Fwuds-Fwhigh)-Fdwe+Fdwi

      END DO
      END DO
      END DO

      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J+Idew

      Su(Inp)=Su(Inp)+Bp(Inp)
      Sv(Inp)=Sv(Inp)+Bt(Inp)
      Sw(Inp)=Sw(Inp)+Bb(Inp)

      End Do !J-Loop
      End Do !I-Loop
      End Do !K-Loop

      do ijk=icst,icen
      bp(ijk)=0.0d0
      bt(ijk)=0.0d0
      bb(ijk)=0.0d0
      end do


!
!.....Calculate Top Cell Face
      idew = nij
      idns = nj
      idtb = 1

      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J
      ine=inp+idew
      inw=inp-idew
      inn=inp+idns
      ins=inp-idns
      inb=inp-idtb
      int=inp+idtb
      inbs=inb-idns
!
!.....Interpolation Factor
      fxe=fz(inp) 
      fxp=1.-fxe

!.....Components Of Distance (From P To Neighbor N) Vector
!.....First
      Xpn=Xc(Ine)-Xc(Inp)
      Ypn=Yc(Ine)-Yc(Inp)
      Zpn=Zc(Ine)-Zc(Inp)

!.....Distance From P To Neighbor N
      Dpn=Sqrt(Xpn**2+Ypn**2+Zpn**2)  
   
!.....Components Of The Unit Vector I_Ksi
      Ixi1=Xpn/Dpn
      Ixi2=Ypn/Dpn
      Ixi3=Zpn/Dpn

!.....Precomputed Face Areas
      Arx=Ar3x(Inp)
      Ary=Ar3y(Inp)
      Arz=Ar3z(Inp)

!.....Cell Face Area
      Are=Sqrt(Arx**2+Ary**2+Arz**2)

!.....Unit Vectors Of The Normal
      nxx=Arx/Are
      nyy=Ary/Are
      nzz=Arz/Are

!.....Angle Between Vectors n And i_xi - We Need Cosine
      costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

!.....Relaxation Factor For Higher-Order Cell Face Gradient
      ! Minimal Correction: Nrelax = +1 :
      !Costn = Costheta
      ! Orthogonal Correction: Nrelax =  0 : 
      Costn = 1.0d0
      ! Over-Relaxed Approach: Nrelax = -1 :
      !Costn = 1./Costheta
      ! In General, Nrelax Can Be Any Signed Integer From Some 
      ! Reasonable Interval [-Nrelax,Nrelax] (Or Maybe Even Real Number): 
      !Costn = Costheta**Nrelax

!.....d_pn . S_f
      Vole=Xpn*Arx+Ypn*Ary+Zpn*Arz

!.....Overrelaxed Correction Vector D2, Where S=Dpn+D2
      D1x = Costn
      D1y = Costn
      D1z = Costn
      !
      D2x = Xpn*Costn
      D2y = Ypn*Costn
      D2z = Zpn*Costn

!.....Cell Face Viscosity
      Game=(Vis(Inp)*Fxp+Vis(Ine)*Fxe)

!++++Velocities At Cell Face Center And Explicit Diffusion Fluxes+++++++

!.....Coordinates Of Cell-Face Center
      Xf = 0.25*(X(Inp)+X(Ins)+X(Inb)+X(Inbs))
      Yf = 0.25*(Y(Inp)+Y(Ins)+Y(Inb)+Y(Inbs))
      Zf = 0.25*(Z(Inp)+Z(Ins)+Z(Inb)+Z(Inbs))

!.....Coordinates Of Point E'
      Xi=Xc(Inp)*Fxp+Xc(Ine)*Fxe
      Yi=Yc(Inp)*Fxp+Yc(Ine)*Fxe
      Zi=Zc(Inp)*Fxp+Zc(Ine)*Fxe

!.....Interpolate Gradients Defined At Cv Centers To Faces
      Duxi = Gradu(1,Inp)*Fxp+Gradu(1,Ine)*Fxe
      Duyi = Gradu(2,Inp)*Fxp+Gradu(2,Ine)*Fxe
      Duzi = Gradu(3,Inp)*Fxp+Gradu(3,Ine)*Fxe

!        |________Ue'_________|_______________Ucorr___________________|
      Ue=U(Inp)*Fxp+U(Ine)*Fxe+(Duxi*(Xf-Xi)+Duyi*(Yf-Yi)+Duzi*(Zf-Zi))

!.....Du/Dx_I Interpolated At Cell Face:
      Duxii = Duxi*D1x + Arx/Vole*( U(Ine)-U(Inp)-Duxi*D2x-Duyi*D2y-Duzi*D2z ) 
      Duyii = Duyi*D1y + Ary/Vole*( U(Ine)-U(Inp)-Duxi*D2x-Duyi*D2y-Duzi*D2z ) 
      Duzii = Duzi*D1z + Arz/Vole*( U(Ine)-U(Inp)-Duxi*D2x-Duyi*D2y-Duzi*D2z ) 

      Dvxi = Gradv(1,Inp)*Fxp+Gradv(1,Ine)*Fxe
      Dvyi = Gradv(2,Inp)*Fxp+Gradv(2,Ine)*Fxe
      Dvzi = Gradv(3,Inp)*Fxp+Gradv(3,Ine)*Fxe

!        |________Ve'_________|_______________Vcorr___________________|
      Ve=V(Inp)*Fxp+V(Ine)*Fxe+(Dvxi*(Xf-Xi)+Dvyi*(Yf-Yi)+Dvzi*(Zf-Zi))

!.....Dv/Dx_I Interpolated At Cell Face:
      Dvxii = Dvxi*D1x + Arx/Vole*( V(Ine)-V(Inp)-Dvxi*D2x-Dvyi*D2y-Dvzi*D2z ) 
      Dvyii = Dvyi*D1y + Ary/Vole*( V(Ine)-V(Inp)-Dvxi*D2x-Dvyi*D2y-Dvzi*D2z ) 
      Dvzii = Dvzi*D1z + Arz/Vole*( V(Ine)-V(Inp)-Dvxi*D2x-Dvyi*D2y-Dvzi*D2z ) 

      Dwxi = Gradw(1,Inp)*Fxp+Gradw(1,Ine)*Fxe
      Dwyi = Gradw(2,Inp)*Fxp+Gradw(2,Ine)*Fxe
      Dwzi = Gradw(3,Inp)*Fxp+Gradw(3,Ine)*Fxe

!        |________We'_________|_______________Wcorr___________________|
      We=W(Inp)*Fxp+W(Ine)*Fxe+(Dwxi*(Xf-Xi)+Dwyi*(Yf-Yi)+Dwzi*(Zf-Zi))

!.....Dw/Dx_I Interpolated At Cell Face:
      Dwxii = Dwxi*D1x + Arx/Vole*( W(Ine)-W(Inp)-Dwxi*D2x-Dwyi*D2y-Dwzi*D2z ) 
      Dwyii = Dwyi*D1y + Ary/Vole*( W(Ine)-W(Inp)-Dwxi*D2x-Dwyi*D2y-Dwzi*D2z ) 
      Dwzii = Dwzi*D1z + Arz/Vole*( W(Ine)-W(Inp)-Dwxi*D2x-Dwyi*D2y-Dwzi*D2z ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     We calculate explicit and implicit diffusion fde and fdi,
!     later we put their difference (fde-fdi) to rhs vector.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ! Explicit Diffussion 
       Fdue = Game*(Duxii+Duxii)*Arx + (Duyii+Dvxii)*Ary + (Duzii+Dwxii)*Arz
       Fdve = Game*(Duyii+Dvxii)*Arx + (Dvyii+Dvyii)*Ary + (Dvzii+Dwyii)*Arz
       Fdwe = Game*(Duzii+Dwxii)*Arx + (Dwyii+Dvzii)*Ary + (Dwzii+Dwzii)*Arz
       
       ! Implicit Diffussion 
       Fdui = Game*Are/Dpn*(Duxi*Xpn+Duyi*Ypn+Duzi*Zpn)
       Fdvi = Game*Are/Dpn*(Dvxi*Xpn+Dvyi*Ypn+Dvzi*Zpn)
       Fdwi = Game*Are/Dpn*(Dwxi*Xpn+Dwyi*Ypn+Dwzi*Zpn)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++End: Velocities At Cell Face Center And Explicit Diffusion Fluxes+++++++


!.....Difusion Coefficient
      De=Game*Are/Dpn

!.....Convection Fluxes - Uds
!     If Flow Goes P=>E Cp=Flcf, Ce=0.
!     If Flow Goes E=>P Ce=Flcf, Cp=0.
      Flcf = F3(Inp)
      Ce = Min(Flcf,Zero) 
      Cp = Max(Flcf,Zero)
!
!.....Coefficients Ae(P) And Aw(E) Due To Uds
!
      At(Inp) = De-Ce
      Ab(Ine) = De+cp

!
!.....EXPLICIT CONVECTIVE FLUXES FOR UDS
!
      Fuuds=Cp*U(Inp)+Ce*U(Ine)
      Fvuds=Cp*V(Inp)+Ce*V(Ine)
      Fwuds=Cp*W(Inp)+Ce*W(Ine)

!=====START CDS SCHEME===============================
      IF(LCDS.EQ.1) THEN
!
!.....EXPLICIT CONVECTIVE FLUXES FOR CDS
!
      Fuhigh=Flcf*Ue
      Fvhigh=Flcf*Ve
      Fwhigh=Flcf*We

!=====END CDS SCHEME================================
      END IF


!=====START MUSCL SCHEME===============================
      IF(LMUSCL.EQ.1) THEN
!
!.....EXPLICIT CONVECTIVE FLUXES FOR MUSCL
!
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
!.....PSI for MUSCL scheme:
!.....If flow goes from P to E
      PSIW1 = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
      PSIW2 = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
      PSIW3 = max(0., min(2.*r3, 0.5*r3+0.5, 2.))
!.....If flow goes from E to P
      PSIE1 = max(0., min(2.*r4, 0.5*r4+0.5, 2.))
      PSIE2 = max(0., min(2.*r5, 0.5*r5+0.5, 2.))
      PSIE3 = max(0., min(2.*r6, 0.5*r6+0.5, 2.))

!.....Ver.3
!......Darwish-Moukalled TVD schemes for unstructured girds, IJHMT, 2003.
      Fuhigh = Ce*(U(Ine) + Fxe*Psie1*(U(Inp)-U(Ine)))+ &
               Cp*(U(Inp) + Fxp*Psiw1*(U(Ine)-U(Inp)))
!       mass flux| bounded interpolation of velocity to face |

      Fvhigh = Ce*(V(Ine) + Fxe*Psie2*(W(Inp)-W(Ine)))+ &
               Cp*(V(Inp) + Fxp*Psiw2*(V(Ine)-V(Inp)))
       
      Fwhigh = Ce*(W(Ine) + Fxe*Psie3*(W(Inp)-W(Ine)))+ &
               Cp*(W(Inp) + Fxp*Psiw3*(W(Ine)-W(Inp)))
!......END: Darwish-Moukalled TVDschemes for unstructured girds, IJHMT, 2003.

!.....END OF BOUNDED HIGH-ORDER SCHEMES
      END IF 

!
!.....EXPLICIT PART AND SOURCES DUE TO DEFFERED CORRECTION
!
      Su(Inp) = Su(Inp)+Gam*(Fuuds-Fuhigh)+Fdue-Fdui
      Sv(Inp) = Sv(Inp)+Gam*(Fvuds-Fvhigh)+Fdve-Fdvi
      Sw(Inp) = Sw(Inp)+Gam*(Fwuds-Fwhigh)+Fdwe-Fdwi
!----------------------------------------------------
!......[Vectorization procedure: ]
!----------------------------------------------------
      Bp(Ine) = -Gam*(Fuuds-Fuhigh)-Fdue+Fdui
      Bt(Ine) = -Gam*(Fvuds-Fvhigh)-Fdve+Fdvi
      Bb(Ine) = -Gam*(Fwuds-Fwhigh)-Fdwe+Fdwi

      END DO
      END DO
      END DO

      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J+Idew

      Su(Inp)=Su(Inp)+Bp(Inp)
      Sv(Inp)=Sv(Inp)+Bt(Inp)
      Sw(Inp)=Sw(Inp)+Bb(Inp)

      End Do !J-Loop
      End Do !I-Loop
      End Do !K-Loop

      do ijk=icst,icen
      bp(ijk)=0.0d0
      bt(ijk)=0.0d0
      bb(ijk)=0.0d0
      end do

      Return
      End



      SUBROUTINE fvm_div_neglap_scalar(fi,gradfi,Gama)
!  
!******************************************************************************
!
!     Fills matrix with coefficients representing implicit FVM discretization 
!     of Convective term: div(psi*Fi).
!     Where psi is mass flux surface field, and Fi is scalar field.
!
!     Inteded for structured 3D grids. Non-orthogonal corrections are included.
!
!     Matrix is stored in diagonal format with seven diagonals.
!     Each diagonal is stored in an array,
!     Off-main diagonnal: ae,aw,an,as,at,ab,
!     Main diagonal: ap.
!     RHS vector is SU.
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!
      Use Types
      Use Parameters
      Use Indexes
      Use Geometry
      Use Coef
      Use Coefb
      Use Variables, only: F1,F2,F3,Vis

      Implicit None
!
!***********************************************************************
!

      Real(Prec), Dimension(Nxyza),   Intent(in) :: Fi     ! Scalar field
      Real(Prec), Dimension(3,Nxyza), Intent(in) :: Gradfi ! Scalar gradient field
      Real(Prec), Dimension(Nxyza),   Intent(in) :: Gama   ! Diffussion coefficient field
!
!     Local Variables
!
      Integer :: Idew, Idns, Idtb, Ifi
      Integer :: I, J, K, Inp
      Integer :: Ine,Ins,Inb,Inbs
      Integer :: Ijk,Lik,Lkk
!@      Integer :: Indx
      !Real(Prec) :: Dxet,Dxzd,Dyet,Dyzd,Dzet,Dzzd
      Real(Prec) :: Gam,Prtr,Fxe,Fxp,Are,Vole
      Real(Prec) :: Viste,Game,De
      Real(Prec) :: Ce,Cp,Fii,Fm
      Real(Prec) :: Xf,Yf,Zf,Xi,Yi,Zi
      Real(Prec) :: Arx,Ary,Arz
      Real(Prec) :: Xpn,Ypn,Zpn
      Real(Prec) :: Nxx,Nyy,Nzz,Ixi1,Ixi2,Ixi3,Dpn,Costheta,Costn
      Real(Prec) :: Fdfie,Fdfii,Fcfie,Fcfii,Ffic,Suadd
      Real(Prec) :: D2x,D2y,D2z,D1x,D1y,D1z
      Real(Prec) :: Dfixi,Dfiyi,Dfizi
      Real(Prec) :: Dfixii,Dfiyii,Dfizii
      Real(Prec) :: R1,R2,Psie,Psiw

!.....Blending (Deffered correction) Coefficient For Convection 
      Gam=Gds(Ifi)

!.....Usually It Is Constant:
      !Prtr=Prtinv(Ifi) !<-passed as input
    
      Dfixi = 0.0d0
      Dfiyi = 0.0d0
      Dfizi = 0.0d0
!-------------------------------------------
!.....Calculate East,Top,North  Cell Face
!-------------------------------------------

!
!.....Calculate East Cell Face
      idew = nj
      idns = 1
      idtb = nij

      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J
      Ine=Inp+Idew
      Ins=Inp-Idns
      Inb=Inp-Idtb
      Inbs=Inb-Idns

!.....Interpolation Factors In First,Second And Third Way
      Fxe=Fx(Inp)
      Fxp=1.-Fxe
!/////Correction For Triple Periodicity://///////////////
!@      Indx=0
!@      If(I.Eq.Nie.And.Idew.Eq.Nj)  Indx=Lk(K)+Li(2)+J
!@      If(J.Eq.Nje.And.Idew.Eq.1)   Indx=Lk(K)+Li(I)+2
!@      If(K.Eq.Nke.And.Idew.Eq.Nij) Indx=Lk(2)+Li(I)+J
!@      If(Indx/=0) Then !................................
!@        Ine=Indx
!@        Fxe=0.5d0 
!@        Fxp=0.5d0
!@      End If  !.........................................
!/////End: Correction For Triple Periodicity://///////////

!.....Components Of Three Vectors
!.....First
      Xpn=Xc(Ine)-Xc(Inp)
      Ypn=Yc(Ine)-Yc(Inp)
      Zpn=Zc(Ine)-Zc(Inp)
!/////Correction For Triple Periodicity://///////////
!@      If(I.Eq.Nie.And.Idew.Eq.Nj)  Xpn=(X(Inp)-Xc(Inp))+(X(Ine)-Xc(Ine))
!@      If(J.Eq.Nje.And.Idew.Eq.1)   Ypn=(Y(Inp)-Yc(Inp))+(Y(Ine)-Yc(Ine))
!@      If(K.Eq.Nke.And.Idew.Eq.Nij) Zpn=(Z(Inp)-Zc(Inp))+(Z(Ine)-Zc(Ine))
!/////End: Correction For Triple Periodicity://///////////

!.....Distance From P To Neighbor N
      Dpn=Sqrt(Xpn**2+Ypn**2+Zpn**2)     
!.....Components Of The Unit Vector I_Ksi
      Ixi1=Xpn/Dpn
      Ixi2=Ypn/Dpn
      Ixi3=Zpn/Dpn

!.....Second
!      Dxet=.5*(X(Inp)-X(Ins)+X(Inb)-X(Inbs))
!      Dyet=.5*(Y(Inp)-Y(Ins)+Y(Inb)-Y(Inbs))
!      Dzet=.5*(Z(Inp)-Z(Ins)+Z(Inb)-Z(Inbs))
!.....Third
!      Dxzd=.5*(X(Inp)-X(Inb)+X(Ins)-X(Inbs))
!      Dyzd=.5*(Y(Inp)-Y(Inb)+Y(Ins)-Y(Inbs))
!      Dzzd=.5*(Z(Inp)-Z(Inb)+Z(Ins)-Z(Inbs))
!.....Second X Third Define  Always Cf Area
!      Arx=Dyet*Dzzd-Dyzd*Dzet
!      Ary=Dxzd*Dzet-Dxet*Dzzd
!      Arz=Dxet*Dyzd-Dyet*Dxzd
!.....Precomputed Face Areas
      Arx=Ar1x(Inp)
      Ary=Ar1y(Inp)
      Arz=Ar1z(Inp)

!.....Cell Face Area
      Are=Sqrt(Arx**2+Ary**2+Arz**2)

!.....Unit Vectors Of The Normal
      nxx=Arx/Are
      nyy=Ary/Are
      nzz=Arz/Are
!.....Angle Between Vectorsa N And I_Xi - We Need Cosine
      Costheta=nxx*Ixi1+nyy*Ixi2+nzz*Ixi3

!.....Relaxation Factor For Higher-Order Cell Face Gradient
      ! Minimal Correction: Nrelax = +1 :
      !Costn = Costheta
      ! Orthogonal Correction: Nrelax =  0 : 
      Costn = 1.0d0
      ! Over-Relaxed Approach: Nrelax = -1 :
      !Costn = 1./Costheta
      ! In General, Nrelax Can Be Any Signed Integer From Some 
      ! Reasonable Interval [-Nrelax,Nrelax] (Or Even Real Number): 
      !Costn = Costheta**Nrelax

!.....First .(Second X Third) = Vol
      Vole=Xpn*Arx+Ypn*Ary+Zpn*Arz
!.....Overrelaxed Correction Vector D2, Where S=Dpn+D2
      D1x = Costn
      D1y = Costn
      D1z = Costn
      !
      D2x = Xpn*Costn
      D2y = Ypn*Costn
      D2z = Zpn*Costn

!.....For Menter Sst Model:
      !If (Sst.Or.Sas.Or.Earsm_Wj.Or.Earsm_M) Then
      !  If(Ifi.Eq.Ite) Prtr=Prtinv_Te(Inp)
      !  If(Ifi.Eq.Ied) Prtr=Prtinv_Ed(Inp)
      !End If
!.....Cell Face Diffussion Coefficient
      !Viste=(Vis(Inp)-Viscos)*Fxp+(Vis(Ine)-Viscos)*Fxe
      !Game=(Viste*Prtr+Viscos)
      Game = Gama(Inp)*Fxp+Gama(Ine)*Fxe

!++++Values Of Fi At Cell Face Center And Explicit Diffusion Fluxes+++++++

!.....Coordinates Of Cell-Face Center
      Xf = 0.25*(X(Inp)+X(Ins)+X(Inb)+X(Inbs))
      Yf = 0.25*(Y(Inp)+Y(Ins)+Y(Inb)+Y(Inbs))
      Zf = 0.25*(Z(Inp)+Z(Ins)+Z(Inb)+Z(Inbs))

!.....Coordinates Of Point E'
      Xi=Xc(Inp)*Fxp+Xc(Ine)*Fxe
      Yi=Yc(Inp)*Fxp+Yc(Ine)*Fxe
      Zi=Zc(Inp)*Fxp+Zc(Ine)*Fxe
!/////Correction For Triple Periodicity://///////////
!@      If(I.Eq.Nie.And.Idew.Eq.Nj)  Xi=Xc(Inp)*Fxp+(Xc(Inp)+Xpn)*Fxe
!@      If(J.Eq.Nje.And.Idew.Eq.1)   Yi=Yc(Inp)*Fxp+(Yc(Inp)+Ypn)*Fxe
!@      If(K.Eq.Nke.And.Idew.Eq.Nij) Zi=Zc(Inp)*Fxp+(Zc(Inp)+Zpn)*Fxe
!/////End: Correction For Triple Periodicity://///////////

!.....Interpolate Gradients Defined At Cv Centers To Faces
      Dfixi = Gradfi(1,Inp)*Fxp+Gradfi(1,Ine)*Fxe
      Dfiyi = Gradfi(2,Inp)*Fxp+Gradfi(2,Ine)*Fxe
      Dfizi = Gradfi(3,Inp)*Fxp+Gradfi(3,Ine)*Fxe


!.....The Cell Face Interpolated Gradient (D Phi / Dx_I)_J:
      Dfixii = Dfixi*D1x + Arx/Vole*( Fi(Ine)-Fi(Inp)-Dfixi*D2x-Dfiyi*D2y-Dfizi*D2z ) 
      Dfiyii = Dfiyi*D1y + Ary/Vole*( Fi(Ine)-Fi(Inp)-Dfixi*D2x-Dfiyi*D2y-Dfizi*D2z ) 
      Dfizii = Dfizi*D1z + Arz/Vole*( Fi(Ine)-Fi(Inp)-Dfixi*D2x-Dfiyi*D2y-Dfizi*D2z ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Explicit Diffusion
      Fdfie = Game*(Dfixii*Arx + Dfiyii*Ary + Dfizii*Arz)   
      ! Implicit Diffussion 
      Fdfii = Game*Are/Dpn*(Dfixi*Xpn+Dfiyi*Ypn+Dfizi*Zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++End: Values Of Fi At Cell Face Center And Explicit Diffusion Fluxes+++++++

!.....Difusion Coefficient
      De=Game*Are/Dpn

!.....Convection Fluxes - Uds
      Fm=F1(Inp)
      Ce=Min(Fm,Zero) 
      Cp=Max(Fm,Zero)

!.....System Matrix Coeffocients
      Ae(Inp)=De-Ce
      Aw(Ine)=De+Cp
!
!---------------------------------------------
!     Convective Fluxes [ Central Differencing Scheme (Cds) ]
!---------------------------------------------
!.....Interpolate Variable Fi Defined At Cv Centers To Face Using Corrected Cds:
!         |__________Ue'_________|_______________Ucorr______________________|
      !Fii=Fi(Inp)*Fxp+Fi(Ine)*Fxe+(Dfixi*(Xf-Xi)+Dfiyi*(Yf-Yi)+Dfizi*(Zf-Zi))
      !Fcfie=Fm*Fii
!---------------------------------------------
!     Darwish-Moukalled Tvd Schemes For Unstructured Girds, Ijhmt, 2003. 
!---------------------------------------------
!     Find R'S - The Gradient Ratio. This Is Universal For All Schemes.
!     If Flow Goes From P To E
      r1 = (2*Gradfi(1,Inp)*Xpn + 2*Gradfi(2,Inp)*Ypn + 2*Gradfi(3,Inp)*Zpn)/(Fi(Ine)-Fi(Inp)) - 1.0d0  
!     if Flow Goes From E To P
      r2 = (2*Gradfi(1,Ine)*Xpn + 2*Gradfi(2,Ine)*Ypn + 2*Gradfi(3,Ine)*Zpn)/(Fi(Inp)-Fi(Ine)) - 1.0d0 
!     Find Psi For [ Muscl ] :
      Psiw = Max(0., Min(2.*r1, 0.5*r1+0.5, 2.))
      Psie = Max(0., Min(2.*r2, 0.5*r2+0.5, 2.))
!     High Order Flux At Cell Face
      Fcfie =  Ce*(Fi(Ine) + Fxe*Psie*(Fi(Inp)-Fi(Ine)))+ &
               Cp*(Fi(Inp) + Fxp*Psiw*(Fi(Ine)-Fi(Inp)))

!.....First Order Upwind Part
      Fcfii=Ce*Fi(Ine)+Cp*Fi(Inp)
!.....Blend High Order Explicit And First Order Upwind 
      Ffic = Gam*(Fcfie-Fcfii)
!-------------------------------------------------------
!.....Explicit Part Of Fluxes
!-------------------------------------------------------
      Suadd = -Ffic+Fdfie-Fdfii 
      Sv(Inp) = Sv(Inp)+Suadd
      Su(Inp) = Su(Inp)+Suadd
!-------------------------------------------------------
      Bp(Ine)=Suadd

      End Do !J-Loop
      End Do !I-Loop
      End Do !K-Loop
!-------------------------------------------------------
!     [Vectorization Routine:  ]
!-------------------------------------------------------
      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J+Idew

       Su(Inp) = Su(Inp)-Bp(Inp)
       Sv(Inp) = Sv(Inp)-Bp(Inp)

      End Do !J-Loop
      End Do !I-Loop
      End Do !K-Loop
!-------------------------------------------------------
      Do Ijk=Icst,Icen
      Bp(Ijk)=0.0d0
      End Do
!-------------------------------------------------------


!
!.....Calculate North Cell Face
      idew = 1
      idns = nij
      idtb = nj

      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J
      Ine=Inp+Idew
      Ins=Inp-Idns
      Inb=Inp-Idtb
      Inbs=Inb-Idns

!.....Interpolation Factors In First,Second And Third Way
      Fxe=Fy(Inp)
      Fxp=1.-Fxe
!/////Correction For Triple Periodicity://///////////////
!@      Indx=0
!@      If(I.Eq.Nie.And.Idew.Eq.Nj)  Indx=Lk(K)+Li(2)+J
!@      If(J.Eq.Nje.And.Idew.Eq.1)   Indx=Lk(K)+Li(I)+2
!@      If(K.Eq.Nke.And.Idew.Eq.Nij) Indx=Lk(2)+Li(I)+J
!@      If(Indx/=0) Then !................................
!@        Ine=Indx
!@        Fxe=0.5d0 
!@        Fxp=0.5d0
!@      End If  !.........................................
!/////End: Correction For Triple Periodicity://///////////

!.....Components Of Three Vectors
!.....First
      Xpn=Xc(Ine)-Xc(Inp)
      Ypn=Yc(Ine)-Yc(Inp)
      Zpn=Zc(Ine)-Zc(Inp)
!/////Correction For Triple Periodicity://///////////
!@      If(I.Eq.Nie.And.Idew.Eq.Nj)  Xpn=(X(Inp)-Xc(Inp))+(X(Ine)-Xc(Ine))
!@      If(J.Eq.Nje.And.Idew.Eq.1)   Ypn=(Y(Inp)-Yc(Inp))+(Y(Ine)-Yc(Ine))
!@      If(K.Eq.Nke.And.Idew.Eq.Nij) Zpn=(Z(Inp)-Zc(Inp))+(Z(Ine)-Zc(Ine))
!/////End: Correction For Triple Periodicity://///////////

!.....Distance From P To Neighbor N
      Dpn=Sqrt(Xpn**2+Ypn**2+Zpn**2)     
!.....Components Of The Unit Vector I_Ksi
      Ixi1=Xpn/Dpn
      Ixi2=Ypn/Dpn
      Ixi3=Zpn/Dpn

!.....Second
!      Dxet=.5*(X(Inp)-X(Ins)+X(Inb)-X(Inbs))
!      Dyet=.5*(Y(Inp)-Y(Ins)+Y(Inb)-Y(Inbs))
!      Dzet=.5*(Z(Inp)-Z(Ins)+Z(Inb)-Z(Inbs))
!.....Third
!      Dxzd=.5*(X(Inp)-X(Inb)+X(Ins)-X(Inbs))
!      Dyzd=.5*(Y(Inp)-Y(Inb)+Y(Ins)-Y(Inbs))
!      Dzzd=.5*(Z(Inp)-Z(Inb)+Z(Ins)-Z(Inbs))
!.....Second X Third Define  Always Cf Area
!      Arx=Dyet*Dzzd-Dyzd*Dzet
!      Ary=Dxzd*Dzet-Dxet*Dzzd
!      Arz=Dxet*Dyzd-Dyet*Dxzd
!.....Precomputed Face Areas
      Arx=Ar2x(Inp)
      Ary=Ar2y(Inp)
      Arz=Ar2z(Inp)

!.....Cell Face Area
      Are=Sqrt(Arx**2+Ary**2+Arz**2)

!.....Unit Vectors Of The Normal
      nxx=Arx/Are
      nyy=Ary/Are
      nzz=Arz/Are
!.....Angle Between Vectorsa N And I_Xi - We Need Cosine
      Costheta=nxx*Ixi1+nyy*Ixi2+nzz*Ixi3

!.....Relaxation Factor For Higher-Order Cell Face Gradient
      ! Minimal Correction: Nrelax = +1 :
      !Costn = Costheta
      ! Orthogonal Correction: Nrelax =  0 : 
      Costn = 1.0d0
      ! Over-Relaxed Approach: Nrelax = -1 :
      !Costn = 1./Costheta
      ! In General, Nrelax Can Be Any Signed Integer From Some 
      ! Reasonable Interval [-Nrelax,Nrelax] (Or Even Real Number): 
      !Costn = Costheta**Nrelax

!.....First .(Second X Third) = Vol
      Vole=Xpn*Arx+Ypn*Ary+Zpn*Arz
!.....Overrelaxed Correction Vector D2, Where S=Dpn+D2
      D1x = Costn
      D1y = Costn
      D1z = Costn
      !
      D2x = Xpn*Costn
      D2y = Ypn*Costn
      D2z = Zpn*Costn

!.....For Menter Sst Model:
      !If (Sst.Or.Sas.Or.Earsm_Wj.Or.Earsm_M) Then
      !  If(Ifi.Eq.Ite) Prtr=Prtinv_Te(Inp)
      !  If(Ifi.Eq.Ied) Prtr=Prtinv_Ed(Inp)
      !End If
!.....Cell Face Diffussion Coefficient
      Viste=(Vis(Inp)-Viscos)*Fxp+(Vis(Ine)-Viscos)*Fxe
      Game=(Viste*Prtr+Viscos)

!++++Values Of Fi At Cell Face Center And Explicit Diffusion Fluxes+++++++

!.....Coordinates Of Cell-Face Center
      Xf = 0.25*(X(Inp)+X(Ins)+X(Inb)+X(Inbs))
      Yf = 0.25*(Y(Inp)+Y(Ins)+Y(Inb)+Y(Inbs))
      Zf = 0.25*(Z(Inp)+Z(Ins)+Z(Inb)+Z(Inbs))

!.....Coordinates Of Point E'
      Xi=Xc(Inp)*Fxp+Xc(Ine)*Fxe
      Yi=Yc(Inp)*Fxp+Yc(Ine)*Fxe
      Zi=Zc(Inp)*Fxp+Zc(Ine)*Fxe
!/////Correction For Triple Periodicity://///////////
!@      If(I.Eq.Nie.And.Idew.Eq.Nj)  Xi=Xc(Inp)*Fxp+(Xc(Inp)+Xpn)*Fxe
!@      If(J.Eq.Nje.And.Idew.Eq.1)   Yi=Yc(Inp)*Fxp+(Yc(Inp)+Ypn)*Fxe
!@      If(K.Eq.Nke.And.Idew.Eq.Nij) Zi=Zc(Inp)*Fxp+(Zc(Inp)+Zpn)*Fxe
!/////End: Correction For Triple Periodicity://///////////

!.....Interpolate Gradients Defined At Cv Centers To Faces
      Dfixi = Gradfi(1,Inp)*Fxp+Gradfi(1,Ine)*Fxe
      Dfiyi = Gradfi(2,Inp)*Fxp+Gradfi(2,Ine)*Fxe
      Dfizi = Gradfi(3,Inp)*Fxp+Gradfi(3,Ine)*Fxe


!.....The Cell Face Interpolated Gradient (D Phi / Dx_I)_J:
      Dfixii = Dfixi*D1x + Arx/Vole*( Fi(Ine)-Fi(Inp)-Dfixi*D2x-Dfiyi*D2y-Dfizi*D2z ) 
      Dfiyii = Dfiyi*D1y + Ary/Vole*( Fi(Ine)-Fi(Inp)-Dfixi*D2x-Dfiyi*D2y-Dfizi*D2z ) 
      Dfizii = Dfizi*D1z + Arz/Vole*( Fi(Ine)-Fi(Inp)-Dfixi*D2x-Dfiyi*D2y-Dfizi*D2z ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Explicit Diffusion
      Fdfie = Game*(Dfixii*Arx + Dfiyii*Ary + Dfizii*Arz)   
      ! Implicit Diffussion 
      Fdfii = Game*Are/Dpn*(Dfixi*Xpn+Dfiyi*Ypn+Dfizi*Zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++End: Values Of Fi At Cell Face Center And Explicit Diffusion Fluxes+++++++

!.....Difusion Coefficient
      De=Game*Are/Dpn

!.....Convection Fluxes - Uds
      Fm=F2(Inp)
      Ce=Min(Fm,Zero) 
      Cp=Max(Fm,Zero)

!.....System Matrix Coeffocients
      An(Inp)=De-Ce
      As(Ine)=De+Cp
!
!---------------------------------------------
!     Convective Fluxes [ Central Differencing Scheme (Cds) ]
!---------------------------------------------
!.....Interpolate Variable Fi Defined At Cv Centers To Face Using Corrected Cds:
!         |__________Ue'_________|_______________Ucorr______________________|
      !Fii=Fi(Inp)*Fxp+Fi(Ine)*Fxe+(Dfixi*(Xf-Xi)+Dfiyi*(Yf-Yi)+Dfizi*(Zf-Zi))
      !Fcfie=Fm*Fii
!---------------------------------------------
!     Darwish-Moukalled Tvd Schemes For Unstructured Girds, Ijhmt, 2003. 
!---------------------------------------------
!     Find R'S - The Gradient Ratio. This Is Universal For All Schemes.
!     If Flow Goes From P To E
      r1 = (2*Gradfi(1,Inp)*Xpn + 2*Gradfi(2,Inp)*Ypn + 2*Gradfi(3,Inp)*Zpn)/(Fi(Ine)-Fi(Inp)) - 1.0d0  
!     if Flow Goes From E To P
      r2 = (2*Gradfi(1,Ine)*Xpn + 2*Gradfi(2,Ine)*Ypn + 2*Gradfi(3,Ine)*Zpn)/(Fi(Inp)-Fi(Ine)) - 1.0d0 
!     Find Psi For [ Muscl ] :
      Psiw = Max(0., Min(2.*r1, 0.5*r1+0.5, 2.))
      Psie = Max(0., Min(2.*r2, 0.5*r2+0.5, 2.))
!     High Order Flux At Cell Face
      Fcfie =  Ce*(Fi(Ine) + Fxe*Psie*(Fi(Inp)-Fi(Ine)))+ &
               Cp*(Fi(Inp) + Fxp*Psiw*(Fi(Ine)-Fi(Inp)))

!.....First Order Upwind Part
      Fcfii=Ce*Fi(Ine)+Cp*Fi(Inp)
!.....Blend High Order Explicit And First Order Upwind 
      Ffic = Gam*(Fcfie-Fcfii)
!-------------------------------------------------------
!.....Explicit Part Of Fluxes
!-------------------------------------------------------
      Suadd = -Ffic+Fdfie-Fdfii 
      Sv(Inp) = Sv(Inp)+Suadd
      Su(Inp) = Su(Inp)+Suadd
!-------------------------------------------------------
      Bp(Ine)=Suadd

      End Do !J-Loop
      End Do !I-Loop
      End Do !K-Loop
!-------------------------------------------------------
!     [Vectorization Routine:  ]
!-------------------------------------------------------
      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J+Idew

       Su(Inp) = Su(Inp)-Bp(Inp)
       Sv(Inp) = Sv(Inp)-Bp(Inp)

      End Do !J-Loop
      End Do !I-Loop
      End Do !K-Loop
!-------------------------------------------------------
      Do Ijk=Icst,Icen
      Bp(Ijk)=0.0d0
      End Do
!-------------------------------------------------------

!
!.....Calculate Top Cell Face
      idew = nij
      idns = nj
      idtb = 1

      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J
      Ine=Inp+Idew
      Ins=Inp-Idns
      Inb=Inp-Idtb
      Inbs=Inb-Idns

!.....Interpolation Factors In First,Second And Third Way
      Fxe=Fz(Inp)
      Fxp=1.-Fxe
!/////Correction For Triple Periodicity://///////////////
!@      Indx=0
!@      If(I.Eq.Nie.And.Idew.Eq.Nj)  Indx=Lk(K)+Li(2)+J
!@      If(J.Eq.Nje.And.Idew.Eq.1)   Indx=Lk(K)+Li(I)+2
!@      If(K.Eq.Nke.And.Idew.Eq.Nij) Indx=Lk(2)+Li(I)+J
!@      If(Indx/=0) Then !................................
!@        Ine=Indx
!@        Fxe=0.5d0 
!@        Fxp=0.5d0
!@      End If  !.........................................
!/////End: Correction For Triple Periodicity://///////////

!.....Components Of Three Vectors
!.....First
      Xpn=Xc(Ine)-Xc(Inp)
      Ypn=Yc(Ine)-Yc(Inp)
      Zpn=Zc(Ine)-Zc(Inp)
!/////Correction For Triple Periodicity://///////////
!@      If(I.Eq.Nie.And.Idew.Eq.Nj)  Xpn=(X(Inp)-Xc(Inp))+(X(Ine)-Xc(Ine))
!@      If(J.Eq.Nje.And.Idew.Eq.1)   Ypn=(Y(Inp)-Yc(Inp))+(Y(Ine)-Yc(Ine))
!@      If(K.Eq.Nke.And.Idew.Eq.Nij) Zpn=(Z(Inp)-Zc(Inp))+(Z(Ine)-Zc(Ine))
!/////End: Correction For Triple Periodicity://///////////

!.....Distance From P To Neighbor N
      Dpn=Sqrt(Xpn**2+Ypn**2+Zpn**2)     
!.....Components Of The Unit Vector I_Ksi
      Ixi1=Xpn/Dpn
      Ixi2=Ypn/Dpn
      Ixi3=Zpn/Dpn

!.....Second
!      Dxet=.5*(X(Inp)-X(Ins)+X(Inb)-X(Inbs))
!      Dyet=.5*(Y(Inp)-Y(Ins)+Y(Inb)-Y(Inbs))
!      Dzet=.5*(Z(Inp)-Z(Ins)+Z(Inb)-Z(Inbs))
!.....Third
!      Dxzd=.5*(X(Inp)-X(Inb)+X(Ins)-X(Inbs))
!      Dyzd=.5*(Y(Inp)-Y(Inb)+Y(Ins)-Y(Inbs))
!      Dzzd=.5*(Z(Inp)-Z(Inb)+Z(Ins)-Z(Inbs))
!.....Second X Third Define  Always Cf Area
!      Arx=Dyet*Dzzd-Dyzd*Dzet
!      Ary=Dxzd*Dzet-Dxet*Dzzd
!      Arz=Dxet*Dyzd-Dyet*Dxzd
!.....Precomputed Face Areas
      Arx=Ar3x(Inp)
      Ary=Ar3y(Inp)
      Arz=Ar3z(Inp)

!.....Cell Face Area
      Are=Sqrt(Arx**2+Ary**2+Arz**2)

!.....Unit Vectors Of The Normal
      nxx=Arx/Are
      nyy=Ary/Are
      nzz=Arz/Are
!.....Angle Between Vectorsa N And I_Xi - We Need Cosine
      Costheta=nxx*Ixi1+nyy*Ixi2+nzz*Ixi3

!.....Relaxation Factor For Higher-Order Cell Face Gradient
      ! Minimal Correction: Nrelax = +1 :
      !Costn = Costheta
      ! Orthogonal Correction: Nrelax =  0 : 
      Costn = 1.0d0
      ! Over-Relaxed Approach: Nrelax = -1 :
      !Costn = 1./Costheta
      ! In General, Nrelax Can Be Any Signed Integer From Some 
      ! Reasonable Interval [-Nrelax,Nrelax] (Or Even Real Number): 
      !Costn = Costheta**Nrelax

!.....First .(Second X Third) = Vol
      Vole=Xpn*Arx+Ypn*Ary+Zpn*Arz
!.....Overrelaxed Correction Vector D2, Where S=Dpn+D2
      D1x = Costn
      D1y = Costn
      D1z = Costn
      !
      D2x = Xpn*Costn
      D2y = Ypn*Costn
      D2z = Zpn*Costn

!.....For Menter Sst Model:
      !If (Sst.Or.Sas.Or.Earsm_Wj.Or.Earsm_M) Then
      !  If(Ifi.Eq.Ite) Prtr=Prtinv_Te(Inp)
      !  If(Ifi.Eq.Ied) Prtr=Prtinv_Ed(Inp)
      !End If
!.....Cell Face Diffussion Coefficient
      Viste=(Vis(Inp)-Viscos)*Fxp+(Vis(Ine)-Viscos)*Fxe
      Game=(Viste*Prtr+Viscos)

!++++Values Of Fi At Cell Face Center And Explicit Diffusion Fluxes+++++++

!.....Coordinates Of Cell-Face Center
      Xf = 0.25*(X(Inp)+X(Ins)+X(Inb)+X(Inbs))
      Yf = 0.25*(Y(Inp)+Y(Ins)+Y(Inb)+Y(Inbs))
      Zf = 0.25*(Z(Inp)+Z(Ins)+Z(Inb)+Z(Inbs))

!.....Coordinates Of Point E'
      Xi=Xc(Inp)*Fxp+Xc(Ine)*Fxe
      Yi=Yc(Inp)*Fxp+Yc(Ine)*Fxe
      Zi=Zc(Inp)*Fxp+Zc(Ine)*Fxe
!/////Correction For Triple Periodicity://///////////
!@      If(I.Eq.Nie.And.Idew.Eq.Nj)  Xi=Xc(Inp)*Fxp+(Xc(Inp)+Xpn)*Fxe
!@      If(J.Eq.Nje.And.Idew.Eq.1)   Yi=Yc(Inp)*Fxp+(Yc(Inp)+Ypn)*Fxe
!@      If(K.Eq.Nke.And.Idew.Eq.Nij) Zi=Zc(Inp)*Fxp+(Zc(Inp)+Zpn)*Fxe
!/////End: Correction For Triple Periodicity://///////////

!.....Interpolate Gradients Defined At Cv Centers To Faces
      Dfixi = Gradfi(1,Inp)*Fxp+Gradfi(1,Ine)*Fxe
      Dfiyi = Gradfi(2,Inp)*Fxp+Gradfi(2,Ine)*Fxe
      Dfizi = Gradfi(3,Inp)*Fxp+Gradfi(3,Ine)*Fxe


!.....The Cell Face Interpolated Gradient (D Phi / Dx_I)_J:
      Dfixii = Dfixi*D1x + Arx/Vole*( Fi(Ine)-Fi(Inp)-Dfixi*D2x-Dfiyi*D2y-Dfizi*D2z ) 
      Dfiyii = Dfiyi*D1y + Ary/Vole*( Fi(Ine)-Fi(Inp)-Dfixi*D2x-Dfiyi*D2y-Dfizi*D2z ) 
      Dfizii = Dfizi*D1z + Arz/Vole*( Fi(Ine)-Fi(Inp)-Dfixi*D2x-Dfiyi*D2y-Dfizi*D2z ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Explicit Diffusion
      Fdfie = Game*(Dfixii*Arx + Dfiyii*Ary + Dfizii*Arz)   
      ! Implicit Diffussion 
      Fdfii = Game*Are/Dpn*(Dfixi*Xpn+Dfiyi*Ypn+Dfizi*Zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++End: Values Of Fi At Cell Face Center And Explicit Diffusion Fluxes+++++++

!.....Difusion Coefficient
      De=Game*Are/Dpn

!.....Convection Fluxes - Uds
      Fm=F3(Inp)
      Ce=Min(Fm,Zero) 
      Cp=Max(Fm,Zero)

!.....System Matrix Coeffocients
      At(Inp)=De-Ce
      Ab(Ine)=De+Cp
!
!---------------------------------------------
!     Convective Fluxes [ Central Differencing Scheme (Cds) ]
!---------------------------------------------
!.....Interpolate Variable Fi Defined At Cv Centers To Face Using Corrected Cds:
!         |__________Ue'_________|_______________Ucorr______________________|
      !Fii=Fi(Inp)*Fxp+Fi(Ine)*Fxe+(Dfixi*(Xf-Xi)+Dfiyi*(Yf-Yi)+Dfizi*(Zf-Zi))
      !Fcfie=Fm*Fii
!---------------------------------------------
!     Darwish-Moukalled Tvd Schemes For Unstructured Girds, Ijhmt, 2003. 
!---------------------------------------------
!     Find R'S - The Gradient Ratio. This Is Universal For All Schemes.
!     If Flow Goes From P To E
      r1 = (2*Gradfi(1,Inp)*Xpn + 2*Gradfi(2,Inp)*Ypn + 2*Gradfi(3,Inp)*Zpn)/(Fi(Ine)-Fi(Inp)) - 1.0d0  
!     if Flow Goes From E To P
      r2 = (2*Gradfi(1,Ine)*Xpn + 2*Gradfi(2,Ine)*Ypn + 2*Gradfi(3,Ine)*Zpn)/(Fi(Inp)-Fi(Ine)) - 1.0d0 
!     Find Psi For [ Muscl ] :
      Psiw = Max(0., Min(2.*r1, 0.5*r1+0.5, 2.))
      Psie = Max(0., Min(2.*r2, 0.5*r2+0.5, 2.))
!     High Order Flux At Cell Face
      Fcfie =  Ce*(Fi(Ine) + Fxe*Psie*(Fi(Inp)-Fi(Ine)))+ &
               Cp*(Fi(Inp) + Fxp*Psiw*(Fi(Ine)-Fi(Inp)))

!.....First Order Upwind Part
      Fcfii=Ce*Fi(Ine)+Cp*Fi(Inp)
!.....Blend High Order Explicit And First Order Upwind 
      Ffic = Gam*(Fcfie-Fcfii)
!-------------------------------------------------------
!.....Explicit Part Of Fluxes
!-------------------------------------------------------
      Suadd = -Ffic+Fdfie-Fdfii 
      Sv(Inp) = Sv(Inp)+Suadd
      Su(Inp) = Su(Inp)+Suadd
!-------------------------------------------------------
      Bp(Ine)=Suadd

      End Do !J-Loop
      End Do !I-Loop
      End Do !K-Loop
!-------------------------------------------------------
!     [Vectorization Routine:  ]
!-------------------------------------------------------
      Do K=2,Nkm
      Do I=2,Nim
      Do J=2,Njm
      Inp=Lk(K)+Li(I)+J+Idew

       Su(Inp) = Su(Inp)-Bp(Inp)
       Sv(Inp) = Sv(Inp)-Bp(Inp)

      End Do !J-Loop
      End Do !I-Loop
      End Do !K-Loop
!-------------------------------------------------------
      Do Ijk=Icst,Icen
      Bp(Ijk)=0.0d0
      End Do
!-------------------------------------------------------

      Return
      End Subroutine fvm_div_neglap_scalar

end module finite_volume_method
