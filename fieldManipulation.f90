    module fieldManipulation
!
!***********************************************************************
!
 
    public

    contains

    pure function volumeWeightedAverage(U) result(wAvgU)
    use types
    use parameters
    use indexes
    use geometry

    implicit none
!
!***********************************************************************
!
!...Output
    real(prec) :: wAvgU

!...Input
    real(prec), dimension(nxyza), intent(in) :: U

!...Locals
    integer :: i,j,k,inp
    real(prec) :: sumvol
  
    sumvol = 0.0_dp
    wAvgU = 0.0_dp  
    do k=3,nkmm
    do i=3,nimm
    do j=3,njmm
    inp = lk(k)+li(i)+j
    sumvol = sumvol + vol(inp)
    wAvgU = wAvgU + (Vol(inp)*U(inp))
    enddo
    enddo
    enddo
    
    wAvgU = wAvgU / sumvol

    end function



!***********************************************************************
!
      subroutine add_random_noise_to_field(phi,percent)
!
!***********************************************************************
!
!     Level is Max perturbation aplitude in percent (%) from given field value.
!   
!     Example usage:
!       call add_random_noise_to_field(U,10)
!
      use types
      use parameters
      use indexes

      implicit none
!
!***********************************************************************
!
      real(prec), dimension(nxyza), intent(inout) :: phi
      integer, intent(in) :: percent
      
      integer :: i, j, k, inp
      real(prec) :: level
      real(prec) :: perturb

      level = dble(percent)

      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      ! Random number based fluctuation of mean profile            
      call init_random_seed()
      call random_number(perturb)
      
      ! perturb is now between 0. and 1., we want it to be from 0 to 2*amplitude
      ! e.g. perturb = 0.9+perturb/5. when Max perturbation is +/- 10% of mean profile
      perturb = ( 1.0_dp - level/100.0_dp ) + perturb * (2*level/100.0_dp)

      phi(inp) = perturb*phi(inp)

      enddo
      enddo
      enddo

      end subroutine


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
!     np: non-periodic
!
!***********************************************************************
!
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
      xpn=xc(ine)-xc(inp)
      ypn=yc(ine)-yc(inp)
      zpn=zc(inp)-zc(inp)

!.....Coordinates of cell-face center - j
      Xf = 0.25*(X(Inp)+X(Ins)+X(Inb)+X(Inbs))
      Yf = 0.25*(Y(Inp)+Y(Ins)+Y(Inb)+Y(Inbs))
      Zf = 0.25*(Z(Inp)+Z(Ins)+Z(Inb)+Z(Inbs))

!.....Coordinates of intersection point - j'
      Xi=Xc(Inp)*Fxp+Xc(Ine)*Fxe
      Yi=Yc(Inp)*Fxp+Yc(Ine)*Fxe
      Zi=Zc(Inp)*Fxp+Zc(Ine)*Fxe

      Ue = 0.5*( (U(Inp)+U(Ine)) +                        &
                 (                                        &
                   (Gradu(1,Inp)+Gradu(1,Ine))*(Xf-Xi) +  &
                   (Gradu(2,Inp)+Gradu(2,Ine))*(Yf-Yi) +  &
                   (Gradu(3,Inp)+Gradu(3,Ine))*(Zf-Zi)    &
                 ) +                                      &
                 (                                        &
                    Gradu(1,Inp)*Xpn*Fxp +                &
                    Gradu(2,Inp)*Ypn*Fxp +                &
                    Gradu(3,Inp)*Zpn*Fxp                  &
                  ) +                                     &
                  (                                       &
                    Gradu(1,Ine)*Xpn*Fxe +                &
                    Gradu(2,Ine)*Ypn*Fxe +                &
                    Gradu(3,Ine)*Zpn*Fxe                  &
                  )                                       &
               )

      end function



      pure function face_value_central(fi,gradfi,inp,idew,idns,idtb) result(ue)
!=======================================================================
!     Calculates face value using values of variables and their gradients
!     at multiple neighbours cell-centers,
!     in least-squares sense.
!=======================================================================
      use types
      use parameters
      use indexes     ! only nj, nij
      use geometry    ! xc,yc,zc

      implicit none

!.....Result
      real(prec) :: ue
!.....Arguments
      real(prec), dimension(nxyza), intent(in) :: fi
      real(prec), dimension(3,nxyza), intent(in) :: gradfi
      integer, intent(in) :: inp,idew,idns,idtb

!.....Locals
      integer :: ins,inb,inbs

!     Locals
      integer innb
      real(prec) :: xf, yf, zf
      real(prec) ::  phi_p, phi_n
      real(prec) :: xcp,ycp,zcp
      real(prec) :: xcn,ycn,zcn
      real(prec) :: gradfi_p_x,gradfi_p_y,gradfi_p_z
      real(prec) :: gradfi_n_x,gradfi_n_y,gradfi_n_z
      real(prec) :: gradfidr

      ins = inp-idns
      inb = inp-idtb
      inbs = inb-idns

      innb = inp + idew


!.....Coordinates of cell-face center - j
      xf = 0.25d0*(x(inp)+x(ins)+x(inb)+x(inbs))
      yf = 0.25d0*(y(inp)+y(ins)+y(inb)+y(inbs))
      zf = 0.25d0*(z(inp)+z(ins)+z(inb)+z(inbs))

!.....Values at cell center:
      phi_p = fi(inp)
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


!.....gradfixdr = (sum(gradphi_nb(i,:)*r_nb2f(i,:)), i=1,n)
      gradfidr=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) &
              +gradfi_n_x*(xf-xcn)+gradfi_n_y*(yf-ycn)+gradfi_n_z*(zf-zcn)

      ue = 0.5d0*( phi_p + phi_n + gradfidr)

      end function


      pure function face_value_2nd_upwind(fi,gradfi,inp,inpv,idew,idns,idtb) result(ue)
!=======================================================================
!     Calculates face value using values of variables and their gradients
!     at multiple neighbours cell-centers,
!     in least-squares sense.
!=======================================================================
      use types
      use parameters
      use indexes     ! only nj, nij
      use geometry    ! xc,yc,zc

      implicit none

!.....Result
      real(prec) :: ue
!.....Arguments
      real(prec), dimension(nxyza), intent(in) :: fi
      real(prec), dimension(3,nxyza), intent(in) :: gradfi
      integer, intent(in) :: inp,inpv,idew,idns,idtb

!.....Locals
      integer :: ine,ins,inb,inbs
      real(prec) :: xf, yf, zf
      real(prec) ::  phi_p
      real(prec) :: xcp,ycp,zcp
      real(prec) :: gradfi_p_x,gradfi_p_y,gradfi_p_z
      real(prec) :: gradfidr

      ine = inpv+idew
      ins = inpv-idns
      inb = inpv-idtb
      inbs = inb-idns

!.....Coordinates of cell-face center - j
      xf = 0.25d0*(x(inpv)+x(ins)+x(inb)+x(inbs))
      yf = 0.25d0*(y(inpv)+y(ins)+y(inb)+y(inbs))
      zf = 0.25d0*(z(inpv)+z(ins)+z(inb)+z(inbs))

!.....Values at cell center:
      phi_p = fi(inp)

      xcp = xc(inp)
      ycp = yc(inp)
      zcp = zc(inp)

      gradfi_p_x = gradfi(1,inp)
      gradfi_p_y = gradfi(2,inp)
      gradfi_p_z = gradfi(3,inp)


!.....gradfixdr = (sum(gradphi_nb(i,:)*r_nb2f(i,:)), i=1,n)
      gradfidr=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp)

      ue = phi_p + gradfidr

      end function



      pure function face_value_muscl(fi,gradfi,inp,inpv,idew,idns,idtb) result(ue)
!=======================================================================
!     Calculates face value using values of variables and their gradients
!     at multiple neighbours cell-centers,
!     in least-squares sense.
!=======================================================================
      use types
      use parameters
      use indexes     ! only nj, nij
      use geometry    ! xc,yc,zc

      implicit none

!.....Result
      real(prec) :: ue
!.....Arguments
      real(prec), dimension(nxyza), intent(in) :: fi
      real(prec), dimension(3,nxyza), intent(in) :: gradfi
      integer, intent(in) :: inp,inpv,idew,idns,idtb

!.....Locals
      integer :: ins,inb,inbs
      integer :: innb
      real(prec) :: xf, yf, zf
      real(prec) ::  phi_p, phi_n,theta
      real(prec) :: xcp,ycp,zcp
      real(prec) :: xcn,ycn,zcn
      real(prec) :: gradfi_p_x,gradfi_p_y,gradfi_p_z
      real(prec) :: gradfi_n_x,gradfi_n_y,gradfi_n_z
      real(prec) :: gradfidr_2nd_upwind,gradfidr_central,face_value_2nd_upwind,face_value_central

      ins = inpv-idns
      inb = inpv-idtb
      inbs = inb-idns

      innb = inp + idew

      ! theta = 1/8
      theta = 0.125d0

!.....Coordinates of cell-face center - j
      xf = 0.25d0*(x(inpv)+x(ins)+x(inb)+x(inbs))
      yf = 0.25d0*(y(inpv)+y(ins)+y(inb)+y(inbs))
      zf = 0.25d0*(z(inpv)+z(ins)+z(inb)+z(inbs))

!.....Values at cell center:
      phi_p = fi(inp)
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


!.....gradfixdr = (sum(gradphi_nb(i,:)*r_nb2f(i,:)), i=1,n)
      gradfidr_2nd_upwind=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) 
      gradfidr_central=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) &
                      +gradfi_n_x*(xf-xcn)+gradfi_n_y*(yf-ycn)+gradfi_n_z*(zf-zcn)

      face_value_2nd_upwind = ( phi_p + gradfidr_2nd_upwind )
      face_value_central = 0.5_dp*( phi_p + phi_n + gradfidr_central)

      ue = theta*face_value_central + (1.0_dp-theta)*face_value_2nd_upwind

      end function


!==== Slope limited interpolation:============================================================


      pure function limited_face_value_2nd_upwind(fi,gradfi,inp,idew,idns,idtb,fimax,fimin) result(ue)
!=======================================================================
!     Calculates face value using values of variables and their gradients
!     at cell-center
!     Cell-centered gradient limited using slope limiter:
!     Wang modified Venkataktirshnan slope limiter
!     Ref.: Z. J. Wang. "A Fast Nested Multi-grid Viscous Flow Solver for Adaptive Cartesian/Quad Grids",
!     International Journal for Numerical Methods in Fluids. 33. 657–680. 2000.
!     The same slope limiter is used in Fluent.
!=======================================================================
      use types
      use parameters
      use indexes     ! only nj, nij
      use geometry    ! xc,yc,zc

      implicit none

!.....Result
      real(prec) :: ue
!.....Arguments
      real(prec), dimension(nxyza), intent(in) :: fi
      real(prec), dimension(3,nxyza), intent(in) :: gradfi
      integer, intent(in) :: inp,idew,idns,idtb
      real(prec), intent(in) :: fimax,fimin

!.....Locals
      integer :: ine,ins,inb,inbs

!     Locals
      real(prec) :: xf, yf, zf
      real(prec) :: phi_p,phi_max,phi_min,deltam,deltap,epsi,slopelimit
      real(prec) :: xcp,ycp,zcp
      real(prec) :: gradfi_p_x,gradfi_p_y,gradfi_p_z
      real(prec) :: gradfidr

      ine = inp+idew
      ins = inp-idns
      inb = inp-idtb
      inbs = inb-idns

!.....Coordinates of cell-face center - j
      xf = 0.25d0*(x(inp)+x(ins)+x(inb)+x(inbs))
      yf = 0.25d0*(y(inp)+y(ins)+y(inb)+y(inbs))
      zf = 0.25d0*(z(inp)+z(ins)+z(inb)+z(inbs))

!.....Values at cell center:
      phi_p = fi(inp)

      xcp = xc(inp)
      ycp = yc(inp)
      zcp = zc(inp)

      gradfi_p_x = gradfi(1,inp)
      gradfi_p_y = gradfi(2,inp)
      gradfi_p_z = gradfi(3,inp)


!.....gradfixdr = (sum(gradphi_nb(i,:)*r_nb2f(i,:)), i=1,n)
      gradfidr=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp)

      ue = phi_p + gradfidr



!:::::Define slope limiter:

!.....max and min values over current cell and neighbors
      phi_max = max(fi(inp),fi(inp-1),fi(inp+1),fi(inp+nj),fi(inp-nj),fi(inp+nij),fi(inp-nij))
      phi_min = min(fi(inp),fi(inp-1),fi(inp+1),fi(inp+nj),fi(inp-nj),fi(inp+nij),fi(inp-nij))

      deltam = ue - phi_p
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


      ue =  phi_p + slopelimit*gradfidr 

      end function

!==== END: Slope limited interpolation.============================================================


!***********************************************************************
!
      pure function von_karman_lengthscale() result(lvk)
!
!***********************************************************************
!
      use types
      use indexes
      use geometry
      use variables
      use boundc
      use gradients
      
      implicit none
!
!***********************************************************************
!
!....Output
     real(prec), dimension(nxyza) :: lvk

!....Input
     !(None)

!....Locals
     integer :: i,j,k,inp
     real(dp) :: uscnd,volr
     real(dp) :: dudxe,dudxw,dudxn,dudxs,dudxt,dudxb
     real(dp) :: dvdxe,dvdxw,dvdxn,dvdxs,dvdxt,dvdxb
     real(dp) :: dwdxe,dwdxw,dwdxn,dwdxs,dwdxt,dwdxb
     real(dp) :: dudye,dudyw,dudyn,dudys,dudyt,dudyb
     real(dp) :: dvdye,dvdyw,dvdyn,dvdys,dvdyt,dvdyb
     real(dp) :: dwdye,dwdyw,dwdyn,dwdys,dwdyt,dwdyb
     real(dp) :: dudze,dudzw,dudzn,dudzs,dudzt,dudzb
     real(dp) :: dvdze,dvdzw,dvdzn,dvdzs,dvdzt,dvdzb
     real(dp) :: dwdze,dwdzw,dwdzn,dwdzs,dwdzt,dwdzb
     real(dp) :: d2udx2,d2udy2,d2udz2
     real(dp) :: d2vdx2,d2vdy2,d2vdz2
     real(dp) :: d2wdx2,d2wdy2,d2wdz2
     
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j
     
      volr = 1./vol(inp)
      

!.....derivatives in  x- direction: 
!      
      dudxe=gradu(1,inp+nj)*fx(inp)+gradu(1,inp)*(1.0d0-fx(inp))
      dudxw=gradu(1,inp)*fx(inp-nj)+gradu(1,inp-nj)*(1.0d0-fx(inp-nj))
      dudxn=gradu(1,inp+1)*fy(inp)+gradu(1,inp)*(1.0d0-fy(inp))
      dudxs=gradu(1,inp)*fy(inp-1)+gradu(1,inp-1)*(1.0d0-fy(inp-1))
      dudxt=gradu(1,inp+nij)*fz(inp)+gradu(1,inp)*(1.0d0-fz(inp))
      dudxb=gradu(1,inp)*fz(inp-nij)+gradu(1,inp-nij)*(1.0d0-fz(inp-nij))
      
      dvdxe=gradv(1,inp+nj)*fx(inp)+gradv(1,inp)*(1.0d0-fx(inp))
      dvdxw=gradv(1,inp)*fx(inp-nj)+gradv(1,inp-nj)*(1.0d0-fx(inp-nj))
      dvdxn=gradv(1,inp+1)*fy(inp)+gradv(1,inp)*(1.0d0-fy(inp))
      dvdxs=gradv(1,inp)*fy(inp-1)+gradv(1,inp-1)*(1.0d0-fy(inp-1))
      dvdxt=gradv(1,inp+nij)*fz(inp)+gradv(1,inp)*(1.0d0-fz(inp))
      dvdxb=gradv(1,inp)*fz(inp-nij)+gradv(1,inp-nij)*(1.0d0-fz(inp-nij))

      dwdxe=gradw(1,inp+nj)*fx(inp)+gradw(1,inp)*(1.0d0-fx(inp))
      dwdxw=gradw(1,inp)*fx(inp-nj)+gradw(1,inp-nj)*(1.0d0-fx(inp-nj))
      dwdxn=gradw(1,inp+1)*fy(inp)+gradw(1,inp)*(1.0d0-fy(inp))
      dwdxs=gradw(1,inp)*fy(inp-1)+gradw(1,inp-1)*(1.0d0-fy(inp-1))
      dwdxt=gradw(1,inp+nij)*fz(inp)+gradw(1,inp)*(1.0d0-fz(inp))
      dwdxb=gradw(1,inp)*fz(inp-nij)+gradw(1,inp-nij)*(1.0d0-fz(inp-nij))
!
!.....derivatives in y- direction:    
      dudye=gradu(2,inp+nj)*fx(inp)+gradu(2,inp)*(1.0d0-fx(inp))
      dudyw=gradu(2,inp)*fx(inp-nj)+gradu(2,inp-nj)*(1.0d0-fx(inp-nj))
      dudyn=gradu(2,inp+1)*fy(inp)+gradu(2,inp)*(1.0d0-fy(inp))
      dudys=gradu(2,inp)*fy(inp-1)+gradu(2,inp-1)*(1.0d0-fy(inp-1))
      dudyt=gradu(2,inp+nij)*fz(inp)+gradu(2,inp)*(1.0d0-fz(inp))
      dudyb=gradu(2,inp)*fz(inp-nij)+gradu(2,inp-nij)*(1.0d0-fz(inp-nij))

      dvdye=gradv(2,inp+nj)*fx(inp)+gradv(2,inp)*(1.0d0-fx(inp))
      dvdyw=gradv(2,inp)*fx(inp-nj)+gradv(2,inp-nj)*(1.0d0-fx(inp-nj))
      dvdyn=gradv(2,inp+1)*fy(inp)+gradv(2,inp)*(1.0d0-fy(inp))
      dvdys=gradv(2,inp)*fy(inp-1)+gradv(2,inp-1)*(1.0d0-fy(inp-1))
      dvdyt=gradv(2,inp+nij)*fz(inp)+gradv(2,inp)*(1.0d0-fz(inp))
      dvdyb=gradv(2,inp)*fz(inp-nij)+gradv(2,inp-nij)*(1.0d0-fz(inp-nij))

      dwdye=gradw(2,inp+nj)*fx(inp)+gradw(2,inp)*(1.0d0-fx(inp))
      dwdyw=gradw(2,inp)*fx(inp-nj)+gradw(2,inp-nj)*(1.0d0-fx(inp-nj))
      dwdyn=gradw(2,inp+1)*fy(inp)+gradw(2,inp)*(1.0d0-fy(inp))
      dwdys=gradw(2,inp)*fy(inp-1)+gradw(2,inp-1)*(1.0d0-fy(inp-1))
      dwdyt=gradw(2,inp+nij)*fz(inp)+gradw(2,inp)*(1.0d0-fz(inp))
      dwdyb=gradw(2,inp)*fz(inp-nij)+gradw(2,inp-nij)*(1.0d0-fz(inp-nij))
!
!.....derivatives in z- direction:      
      dudze=gradu(3,inp+nj)*fx(inp)+gradu(3,inp)*(1.0d0-fx(inp))
      dudzw=gradu(3,inp)*fx(inp-nj)+gradu(3,inp-nj)*(1.0d0-fx(inp-nj))
      dudzn=gradu(3,inp+1)*fy(inp)+gradu(3,inp)*(1.0d0-fy(inp))
      dudzs=gradu(3,inp)*fy(inp-1)+gradu(3,inp-1)*(1.0d0-fy(inp-1))
      dudzt=gradu(3,inp+nij)*fz(inp)+gradu(3,inp)*(1.0d0-fz(inp))
      dudzb=gradu(3,inp)*fz(inp-nij)+gradu(3,inp-nij)*(1.0d0-fz(inp-nij))

      dvdze=gradv(3,inp+nj)*fx(inp)+gradv(3,inp)*(1.0d0-fx(inp))
      dvdzw=gradv(3,inp)*fx(inp-nj)+gradv(3,inp-nj)*(1.0d0-fx(inp-nj))
      dvdzn=gradv(3,inp+1)*fy(inp)+gradv(3,inp)*(1.0d0-fy(inp))
      dvdzs=gradv(3,inp)*fy(inp-1)+gradv(3,inp-1)*(1.0d0-fy(inp-1))
      dvdzt=gradv(3,inp+nij)*fz(inp)+gradv(3,inp)*(1.0d0-fz(inp))
      dvdzb=gradv(3,inp)*fz(inp-nij)+gradv(3,inp-nij)*(1.0d0-fz(inp-nij))

      dwdze=gradw(3,inp+nj)*fx(inp)+gradw(3,inp)*(1.0d0-fx(inp))
      dwdzw=gradw(3,inp)*fx(inp-nj)+gradw(3,inp-nj)*(1.0d0-fx(inp-nj))
      dwdzn=gradw(3,inp+1)*fy(inp)+gradw(3,inp)*(1.0d0-fy(inp))
      dwdzs=gradw(3,inp)*fy(inp-1)+gradw(3,inp-1)*(1.0d0-fy(inp-1))
      dwdzt=gradw(3,inp+nij)*fz(inp)+gradw(3,inp)*(1.0d0-fz(inp))
      dwdzb=gradw(3,inp)*fz(inp-nij)+gradw(3,inp-nij)*(1.0d0-fz(inp-nij))
!
!.....second derivatives:
      d2udx2 = ((dudxe-dudxw)*ar1x(inp)+(dudxn-dudxs)*ar2x(inp)+ &
                (dudxt-dudxb)*ar3x(inp))*volr 
      d2udy2 = ((dudye-dudyw)*ar1y(inp)+(dudyn-dudys)*ar2y(inp)+ &
                (dudyt-dudyb)*ar3y(inp))*volr
      d2udz2 = ((dudze-dudzw)*ar1z(inp)+(dudzn-dudzs)*ar2z(inp)+ &
                (dudzt-dudzb)*ar3z(inp))*volr
!---------------
      d2vdx2 = ((dvdxe-dvdxw)*ar1x(inp)+(dvdxn-dvdxs)*ar2x(inp)+ &
                (dvdxt-dvdxb)*ar3x(inp))*volr
      d2vdy2 = ((dvdye-dvdyw)*ar1y(inp)+(dvdyn-dvdys)*ar2y(inp)+ &
                (dvdyt-dvdyb)*ar3y(inp))*volr
      d2vdz2 = ((dvdze-dvdzw)*ar1z(inp)+(dvdzn-dvdzs)*ar2z(inp)+ &
                (dvdzt-dvdzb)*ar3z(inp))*volr
!---------------
      d2wdx2 = ((dwdxe-dwdxw)*ar1x(inp)+(dwdxn-dwdxs)*ar2x(inp)+ &
                (dwdxt-dwdxb)*ar3x(inp))*volr
      d2wdy2 = ((dwdye-dwdyw)*ar1y(inp)+(dwdyn-dwdys)*ar2y(inp)+ &
                (dwdyt-dwdyb)*ar3y(inp))*volr
      d2wdz2 = ((dwdze-dwdzw)*ar1z(inp)+(dwdzn-dwdzs)*ar2z(inp)+ &
                (dwdzt-dwdzb)*ar3z(inp))*volr
!---------------


!.....2nd velocity derivative generalized to 3d using the magnitude of
!     velocity laplacian
      uscnd = sqrt((d2udx2+d2udy2+d2udz2)**2+ &
                   (d2vdx2+d2vdy2+d2vdz2)**2+ &
                   (d2wdx2+d2wdy2+d2wdz2)**2) 
                 
!.....von karman length scale
      lvk(inp) = cappa*strain(inp)/uscnd
                 
      end do
      end do
      end do

      end function

    end module fieldManipulation
