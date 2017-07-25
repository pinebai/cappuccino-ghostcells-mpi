      MODULE INTERPOLATION
  
      PUBLIC

      CONTAINS

!***********************************************************************
!
      pure function face_interpolated(u,gradu,inp,idew,idns,idtb,fxp,fxe) result(ue)
!
!***********************************************************************
!
!     Variable interpolated at cell face center with non-orthogonality 
!     correction.
!     This is broken down version of two sided interpolation with cell
!     values at two sides and corresponding gradients.
!
!***********************************************************************
!
      use types
      use types
      use parameters
      use indexes
      use geometry

      implicit none
!
!***********************************************************************
!

!.....Result
      real(prec) :: ue
!.....Arguments
      real(prec), dimension(nxyza), intent(in) :: u
      real(prec), dimension(3,nxyza), intent(in) :: gradu
      integer, intent(in) :: inp,idew,idns,idtb
      real(prec), intent(in) :: fxe, fxp
!.....Locals
      integer :: ine,ins,inb,inbs
      real(prec) :: xpn,ypn,zpn,xf,yf,zf, xi,yi,zi

      INE = INP+IDEW
      INS = INP-IDNS
      INB = INP-IDTB
      INBS = INB-IDNS

!.....Distance vector between cell centers
      XPN=XC(INE)-XC(INP)
      YPN=YC(INE)-YC(INP)
      ZPN=ZC(INP)-ZC(INP)

!.....Coordinates of cell-face center - j
      XF = 0.25*(X(INP)+X(INS)+X(INB)+X(INBS))
      YF = 0.25*(Y(INP)+Y(INS)+Y(INB)+Y(INBS))
      ZF = 0.25*(Z(INP)+Z(INS)+Z(INB)+Z(INBS))

!.....Coordinates of intersection point - j'
      XI=XC(INP)*FXP+XC(INE)*FXE
      YI=YC(INP)*FXP+YC(INE)*FXE
      ZI=ZC(INP)*FXP+ZC(INE)*FXE

      UE = 0.5*( (U(INP)+U(INE)) +                        &
                 (                                        &
                   (gradU(1,INP)+gradU(1,INE))*(XF-XI) +  &
                   (gradU(2,INP)+gradU(2,INE))*(YF-YI) +  &
                   (gradU(3,INP)+gradU(3,INE))*(ZF-ZI)    &
                 ) +                                      &
                 (                                        &
                    gradU(1,INP)*XPN*FXP +                &
                    gradU(2,INP)*YPN*FXP +                &
                    gradU(3,INP)*ZPN*FXP                  &
                  ) +                                     &
                  (                                       &
                    gradU(1,INE)*XPN*FXE +                &
                    gradU(2,INE)*YPN*FXE +                &
                    gradU(3,INE)*ZPN*FXE                  &
                  )                                       &
               )

      end function


      SUBROUTINE face_velocities(Ue,Un,Us,Ut,Ub,    &
                                 Ve,Vn,Vs,Vt,Vb,    &
                                 We,Wn,Ws,Wt,Wb,    &
                                 inp,idew,idns,idtb)
!=======================================================================
!     Calculates face value using values of variables and their gradients
!     at multiple neighbours cell-centers,
!     in least-squares sense.
!=======================================================================
      USE TYPES
      USE PARAMETERS
      USE INDEXES     ! only NJ, NIJ
      USE GEOMETRY    ! xc,yc,zc
      USE VARIABLES
      USE GRADIENTS

      IMPLICIT NONE

      REAL(PREC), INTENT(OUT) ::   Ue,Un,Us,Ut,Ub,  &
                                   Ve,Vn,Vs,Vt,Vb,  &
                                   We,Wn,Ws,Wt,Wb
      !REAL(PREC),DIMENSION(NXYZA), INTENT(IN) :: U,V,W
      !REAL(PREC),DIMENSION(3,NXYZA), INTENT(IN) :: gradU,gradV,gradW
      INTEGER, INTENT(IN) :: inp,idew,idns,idtb

!
!     Locals
!
      INTEGER v1, v2, v3, v4
      REAL(PREC) :: xf, yf, zf
      REAL(PREC) :: xcp,ycp,zcp
      REAL(PREC) :: gradfidr

!.....Cell center coordinates
      xcp = xc(inp)
      ycp = yc(inp)
      zcp = zc(inp)

!=====First let's calculate values at mid-face point 'e'==========================

!.....Vertices used in calculation
      v1 = inp
      v2 = inp-idtb
      v3 = inp-idns
      v4 = v2-idns
!.....Point coordinates
      xf = 0.25*(x(v1)+x(v2)+x(v3)+x(v4))
      yf = 0.25*(y(v1)+y(v2)+y(v3)+y(v4))
      zf = 0.25*(z(v1)+z(v2)+z(v3)+z(v4))
!.....Inner product (gradPhi,dr)
      gradfidr=gradU(1,inp)*(xf-xcp)+gradU(2,inp)*(yf-ycp)+gradU(3,inp)*(zf-zcp) 
      Ue = U(inp) + gradfidr 
      gradfidr=gradV(1,inp)*(xf-xcp)+gradV(2,inp)*(yf-ycp)+gradV(3,inp)*(zf-zcp) 
      Ve = V(inp) + gradfidr 
      gradfidr=gradW(1,inp)*(xf-xcp)+gradW(2,inp)*(yf-ycp)+gradW(3,inp)*(zf-zcp) 
      We = W(inp) + gradfidr
!==================================================================================  

!=====Calculate values at mid-egde point 'n'=======================================

!.....Vertices used in calculation
      v1 = inp
      v2 = inp-idtb
!.....Point coordinates
      xf = 0.5*(x(v1)+x(v2))
      yf = 0.5*(y(v1)+y(v2))
      zf = 0.5*(z(v1)+z(v2))
!.....Inner product (gradPhi,dr)
      gradfidr=gradU(1,inp)*(xf-xcp)+gradU(2,inp)*(yf-ycp)+gradU(3,inp)*(zf-zcp) 
      Un = U(inp) + gradfidr 
      gradfidr=gradV(1,inp)*(xf-xcp)+gradV(2,inp)*(yf-ycp)+gradV(3,inp)*(zf-zcp) 
      Vn = V(inp) + gradfidr 
      gradfidr=gradW(1,inp)*(xf-xcp)+gradW(2,inp)*(yf-ycp)+gradW(3,inp)*(zf-zcp) 
      Wn = W(inp) + gradfidr
!==================================================================================

!=====Calculate values at mid-egde point 's'=======================================

!.....Vertices used in calculation
      v1 = inp-idns
      v2 = inp-idtb-idns
!.....Point coordinates
      xf = 0.5*(x(v1)+x(v2))
      yf = 0.5*(y(v1)+y(v2))
      zf = 0.5*(z(v1)+z(v2))
!.....Inner product (gradPhi,dr)
      gradfidr=gradU(1,inp)*(xf-xcp)+gradU(2,inp)*(yf-ycp)+gradU(3,inp)*(zf-zcp) 
      Us = U(inp) + gradfidr 
      gradfidr=gradV(1,inp)*(xf-xcp)+gradV(2,inp)*(yf-ycp)+gradV(3,inp)*(zf-zcp) 
      Vs = V(inp) + gradfidr 
      gradfidr=gradW(1,inp)*(xf-xcp)+gradW(2,inp)*(yf-ycp)+gradW(3,inp)*(zf-zcp) 
      Ws = W(inp) + gradfidr
!==================================================================================

!=====Calculate values at mid-egde point 't'=======================================

!.....Vertices used in calculation
      v1 = inp
      v2 = inp-idns
!.....Point coordinates
      xf = 0.5*(x(v1)+x(v2))
      yf = 0.5*(y(v1)+y(v2))
      zf = 0.5*(z(v1)+z(v2))
!.....Inner product (gradPhi,dr)
      gradfidr=gradU(1,inp)*(xf-xcp)+gradU(2,inp)*(yf-ycp)+gradU(3,inp)*(zf-zcp) 
      Ut = U(inp) + gradfidr 
      gradfidr=gradV(1,inp)*(xf-xcp)+gradV(2,inp)*(yf-ycp)+gradV(3,inp)*(zf-zcp) 
      Vt = V(inp) + gradfidr 
      gradfidr=gradW(1,inp)*(xf-xcp)+gradW(2,inp)*(yf-ycp)+gradW(3,inp)*(zf-zcp) 
      Wt = W(inp) + gradfidr
!==================================================================================

!=====Calculate values at mid-egde point 'b'=======================================

!.....Vertices used in calculation
      v1 = inp-idtb
      v2 = inp-idtb-idns
!.....Point coordinates
      xf = 0.5*(x(v1)+x(v2))
      yf = 0.5*(y(v1)+y(v2))
      zf = 0.5*(z(v1)+z(v2))
!.....Inner product (gradPhi,dr)
      gradfidr=gradU(1,inp)*(xf-xcp)+gradU(2,inp)*(yf-ycp)+gradU(3,inp)*(zf-zcp) 
      Ub = U(inp) + gradfidr 
      gradfidr=gradV(1,inp)*(xf-xcp)+gradV(2,inp)*(yf-ycp)+gradV(3,inp)*(zf-zcp) 
      Vb = V(inp) + gradfidr 
      gradfidr=gradW(1,inp)*(xf-xcp)+gradW(2,inp)*(yf-ycp)+gradW(3,inp)*(zf-zcp) 
      Wb = W(inp) + gradfidr
!==================================================================================   

      RETURN
      END SUBROUTINE
!
!
!
      DOUBLE PRECISION FUNCTION face_value(fi,inp,face,gradfi)
!=======================================================================
!     Calculates face value using values of variables and their gradients
!     at multiple neighbours cell-centers,
!     in least-squares sense.
!=======================================================================
      USE TYPES
      USE PARAMETERS
      USE INDEXES     ! only NJ, NIJ
      USE GEOMETRY    ! xc,yc,zc

      IMPLICIT NONE

!     Input
      CHARACTER face
      INTEGER inp
      REAL(PREC),DIMENSION(NXYZA) :: fi
      REAL(PREC),DIMENSION(3,NXYZA) :: gradfi

!     Locals
      INTEGER innb
      INTEGER v1, v2, v3, v4
      REAL(PREC) :: xf, yf, zf
      REAL(PREC) ::  phi_p, phi_n
      REAL(PREC) :: xcp,ycp,zcp
      REAL(PREC) :: xcn,ycn,zcn
      REAL(PREC) :: gradfi_p_x,gradfi_p_y,gradfi_p_z
      REAL(PREC) :: gradfi_n_x,gradfi_n_y,gradfi_n_z
      REAL(PREC) :: nr
      REAL(PREC) :: gradfidr

      IF (face.eq.'e') THEN
        innb = inp + nj
        v1 = inp
        v2 = inp-nij
        v3 = inp-1
        v4 = v2-1
      ELSEIF (face.eq.'w') THEN
        innb = inp - nj
        v1 = inp-nj
        v2 = inp-nj-nij
        v3 = inp-nj-1
        v4 = v2-1
      ELSEIF (face.eq.'n') THEN
        innb = inp+1
        v1 = inp
        v2 = inp - nij
        v3 = inp - nj
        v4 = v3 - nij
      ELSEIF (face.eq.'s') THEN
        innb = inp-1
        v1 = inp-1
        v2 = inp-1-nij
        v3 = inp-1-nj
        v4 = v3 - nij
      ELSEIF (face.eq.'t') THEN
        innb = inp + nij
        v1 = inp
        v2 = inp-1
        v3 = inp-nj
        v4 = inp-nj-1
      ELSEIF (face.eq.'b') THEN
        innb = inp - nij
        v1 = inp-nij
        v2 = inp-nij-1
        v3 = inp-nij-nj
        v4 = inp-nij-nj-1
      ENDIF

!.....find face-center coordinates......................................
      xf = 0.25*(x(v1)+x(v2)+x(v3)+x(v4))
      yf = 0.25*(y(v1)+y(v2)+y(v3)+y(v4))
      zf = 0.25*(z(v1)+z(v2)+z(v3)+z(v4))

!.....Values at cell center's of neighbouring cells:
      phi_p = fi(inp)

!.....for higher precision..............................................
      phi_n = fi(innb)

      xcp = xc(inp)
      ycp = yc(inp)
      zcp = zc(inp)

      xcn = xc(innb)
      ycn = yc(innb)
      zcn = zc(innb)

      gradfi_p_x = gradfi(1,inp)
      gradfi_p_y = gradfi(2,inp)
      gradfi_p_z = gradfi(3,inp)

      gradfi_n_x = gradfi(1,innb)
      gradfi_n_y = gradfi(2,innb)
      gradfi_n_z = gradfi(3,innb)

      nr = 0.50D0

!.....gradfixdr = (sum(gradphi_nb(i,:)*r_nb2f(i,:)), i=1,n)
!@      gradfidr=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) 
      gradfidr=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) &
              +gradfi_n_x*(xf-xcn)+gradfi_n_y*(yf-ycn)+gradfi_n_z*(zf-zcn)

!@      face_value = ( phi_p + gradfidr )
      face_value = nr*( phi_p + phi_n + gradfidr)

      END FUNCTION

      DOUBLE PRECISION FUNCTION slope_limited_face_value(fi,inp,face,gradfi,fimax,fimin) result(face_value)
!=======================================================================
!     Calculates face value using values of variables and their gradients
!     at cell-center
!     Cell-centered gradient limited using slope limiter:
!     Wang modified Venkataktirshnan slope limiter
!     Ref.: Z. J. Wang. "A Fast Nested Multi-grid Viscous Flow Solver for Adaptive Cartesian/Quad Grids",
!     International Journal for Numerical Methods in Fluids. 33. 657–680. 2000.
!     The same slope limiter is used in Fluent.
!=======================================================================
      USE TYPES
      USE PARAMETERS
      USE INDEXES     ! only NJ, NIJ
      USE GEOMETRY    ! xc,yc,zc,vol

      IMPLICIT NONE

!     Input
      CHARACTER face
      INTEGER inp
      REAL(PREC),DIMENSION(NXYZA) :: fi
      REAL(PREC),DIMENSION(3,NXYZA) :: gradfi

!     Locals
      INTEGER innb
      INTEGER v1, v2, v3, v4
      REAL(PREC) :: xf, yf, zf
      REAL(PREC) :: phi_p
      REAL(PREC) :: gradfidr,slopelimit
      REAL(PREC) :: deltam,deltap,epsi
      REAL(PREC) :: phi_max,phi_min
      REAL(PREC) :: fimax,fimin

      IF (face.eq.'e') THEN
        innb = inp + nj
        v1 = inp
        v2 = inp-nij
        v3 = inp-1
        v4 = v2-1
      ELSEIF (face.eq.'w') THEN
        innb = inp - nj
        v1 = inp-nj
        v2 = inp-nj-nij
        v3 = inp-nj-1
        v4 = v2-1
      ELSEIF (face.eq.'n') THEN
        innb = inp+1
        v1 = inp
        v2 = inp - nij
        v3 = inp - nj
        v4 = v3 - nij
      ELSEIF (face.eq.'s') THEN
        innb = inp-1
        v1 = inp-1
        v2 = inp-1-nij
        v3 = inp-1-nj
        v4 = v3 - nij
      ELSEIF (face.eq.'t') THEN
        innb = inp + nij
        v1 = inp
        v2 = inp-1
        v3 = inp-nj
        v4 = inp-nj-1
      ELSEIF (face.eq.'b') THEN
        innb = inp - nij
        v1 = inp-nij
        v2 = inp-nij-1
        v3 = inp-nij-nj
        v4 = inp-nij-nj-1
      ENDIF

!.....find face-center coordinates......................................
      xf = 0.25*(x(v1)+x(v2)+x(v3)+x(v4))
      yf = 0.25*(y(v1)+y(v2)+y(v3)+y(v4))
      zf = 0.25*(z(v1)+z(v2)+z(v3)+z(v4))

!.....Values at cell center:
      phi_p = fi(inp)

!.....gradfixdr = (sum(gradphi_nb(i,:)*r_nb2f(i,:)), i=1,n)
      gradfidr=gradfi(1,inp)*(xf-xc(inp))+gradfi(2,inp)*(yf-yc(inp))+gradfi(3,inp)*(zf-zc(inp)) 

!.....Find unlimited value:
      face_value =  phi_p + gradfidr 

!:::::Define slope limiter:

!.....max and min values over current cell and neighbors
      phi_max = max(fi(inp),fi(inp-1),fi(inp+1),fi(inp+nj),fi(inp-nj),fi(inp+nij),fi(inp-nij))
      phi_min = min(fi(inp),fi(inp-1),fi(inp+1),fi(inp+nj),fi(inp-nj),fi(inp+nij),fi(inp-nij))

      deltam = face_value - phi_p
      if (deltam .gt. 0.0d0) then
          deltap = phi_max-phi_p
      else
          deltap = phi_min-phi_p
      endif

!.....Original Venkatakrishnan K=[0,?], we take fixed K=0.05
!      epsi = (0.05*vol(inp))**3 
  
!.....Wang proposition for epsilon
      epsi = 0.05*( fimax-fimin )
      epsi = epsi*epsi

      slopelimit = 1./(deltam+small) *((deltap+epsi)*deltam+2*deltam**2*deltap) &
                                     /(deltap**2+2*deltam**2+deltap*deltam+epsi+small)


      face_value =  phi_p + slopelimit*gradfidr 

      END FUNCTION



      END MODULE INTERPOLATION

