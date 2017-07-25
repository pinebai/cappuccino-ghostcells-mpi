!***********************************************************************
!
      subroutine modvis
!
!***********************************************************************
!
!     Update effective viscosity = molecular dynamic visc. + eddy visc.
!     we need fresh vel. gradients for some models - we got them in calcscm
!     when they are evaluated after vel. correction in calcp - these 
!     are the freshest they can be in this outer (simple) interation.
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use geometry
      use coef
      use variables
      use bc
      use boundc
      use title_mod
      use buoy
      use time_mod
      use coefb
      use gradients
      use omega_turb_models
      use fieldmanipulation ! volume weighted average function

      implicit none
!
!***********************************************************************
!
      integer :: i, j, k, inp
      real(prec) :: visold
      real(prec) :: a0,stild, wrlzb,ffi, ass, ust, cmur          ! for realizable k-eps
      real(prec) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz ! gradients
      real(prec) :: s11,s12,s13,s21,s22,s23,s31,s32,s33, &       ! symmetric and assymetric part of vel. grad. tensor - 
                    w12,w13,w23                                  ! strain and rotation tensor components                  
      real(prec) :: fmu
      real(prec) :: wldist,etha,f2_sst,alphast


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(stdkeps.or.rng) then
!====================================================
!.....for standard k-epsilon model
!====================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm
      inp=lk(k)+li(i)+j
      visold=vis(inp)

      vis(inp)=viscos + den(inp)*cmu*te(inp)**2/(ed(inp)+small)
      vis(inp)=urf(ivis)*vis(inp)+(1.-urf(ivis))*visold

      enddo
      enddo
      enddo
!====================================================
      endif


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(durbin) then
!====================================================
!.....for k-epsilon model + durbin time-scale limiter
!====================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm
      inp=lk(k)+li(i)+j
      visold=vis(inp)

      timelimit(inp)=min(te(inp)/(ed(inp)+small), &
       .6/(cmu*sqrt(6.)*strain(inp)))
      vis(inp)=viscos + densit*te(inp)*cmu*timelimit(inp)
      vis(inp)=urf(ivis)*vis(inp)+(1.-urf(ivis))*visold

      enddo
      enddo
      enddo
!====================================================
      endif



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(lowre_lb) then
!====================================================
!.....for lam-bremhorst low re k-epsilon model
!====================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm
      inp=lk(k)+li(i)+j
      visold=vis(inp)
      ret(inp)=den(inp)*te(inp)**2/(viscos*ed(inp)+small) 
      !
      ! launder-sharma :
      !
      !fmu = exp(-3.4d0/(1.+0.02d0*ret(inp))**2)
      !
      ! lam-bremhorst :
      !
      fmu = ret(inp) !<<< we used ret array to write f_mu damping function values in calcscm 

      vis(inp)=viscos + den(inp)*cmu*fmu*te(inp)**2/(ed(inp)+small)
      vis(inp)=urf(ivis)*vis(inp)+(1.-urf(ivis))*visold

      enddo
      enddo
      enddo
!====================================================
      endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(realizable) then
!====================================================
!.....for realizable re k-epsilon model by shih et.al.
!====================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      a0 = 4.0d0   ! original ref. a0 = 4.0; fluent puts a0=4.04!

      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      visold=vis(inp)

!     for easier manipulation
      dudx = gradu(1,inp)
      dudy = gradu(2,inp)
      dudz = gradu(3,inp)
      
      dvdx = gradv(1,inp)
      dvdy = gradv(2,inp)
      dvdz = gradv(3,inp)

      dwdx = gradw(1,inp)
      dwdy = gradw(2,inp)
      dwdz = gradw(3,inp)

!.....find strain rate tensor
!     [s_ij]: |s_ij|=sqrt[2s_ij s_ij]
      s11=dudx
      s12=0.5*(dudy+dvdx)
      s13=0.5*(dudz+dwdx)
      s22=dvdy
      s23=0.5*(dvdz+dwdy) 
      s33=dwdz
      s21=s12  ! symmetry
      s31=s13
      s32=s23
!.....find antisymmetric part of velocity gradient tensor
!     [om_ij]: |om_ij|=sqrt[2 om_ij om_ij] 
      w12=0.5*(dudy - dvdx)
      w13=0.5*(dudz - dwdx)
      w23=0.5*(dvdz - dwdy)



!.....find stilda = sqrt (sij*sij) note there's no factor of 2!!! like in s=sqrt (2*sij*sij)
      stild=sqrt(s11**2+s22**2+s33**2 + 2*(s12**2+s13**2+s23**2))

!     w = sij*sjk*ski/s
      wrlzb = (s11*s11*s11+s11*s12*s21+s11*s13*s31+ &
               s12*s21*s11+s12*s22*s21+s12*s23*s31+ &
               s13*s31*s11+s13*s32*s21+s13*s33*s31+ &
               s21*s11*s12+s21*s12*s22+s21*s13*s32+ &
               s22*s21*s12+s22*s22*s22+s22*s23*s32+ &
               s23*s31*s12+s23*s32*s22+s23*s33*s32+ &
               s31*s11*s13+s31*s12*s23+s31*s13*s33+ &
               s32*s21*s13+s32*s22*s23+s32*s23*s33+ &
               s33*s31*s13+s33*s32*s23+s33*s33*s33)/(stild**3)

      ffi = 1./3. * acos(max(-1.0d0, min( sqrt(6.0d0)*wrlzb ,1.0d0))) ! we make sure the argument of arccos lies in interval [-1,1]
      ass = dsqrt(6.0d0)*cos(ffi)
      ust = dsqrt(s11**2+s22**2+s33**2 + 2*(s12**2+s13**2+s23**2 + w12**2+w13**2+w23**2))

      cmur = 1.0d0/(a0 + ass*ust*te(inp)/(ed(inp) + small))

      vis(inp)=viscos + den(inp)*cmur*te(inp)**2/(ed(inp)+small)

      vis(inp)=urf(ivis)*vis(inp)+(1.-urf(ivis))*visold

      end do !!i-loop
      end do !!j-loop
      end do !!k-loop
!====================================================
      endif


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(sst.or.sas) then
!====================================================
!.....for menter k-omega shear-stress transport (sst)
!====================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm
      inp=lk(k)+li(i)+j
      visold=vis(inp)

      wldist = walldistance(inp)

!=====find etha=========================================================
      etha=max(2*sqrt(te(inp))/(cmu*wldist*ed(inp)), &
               (500*viscos/den(inp))/(wldist**2*ed(inp))) 
!=====find f2===========================================================
      f2_sst = tanh(etha*etha)
      vis(inp)=viscos  &
              +den(inp)*a1*te(inp)/(max(a1,ed(inp), strain(inp)*f2_sst))
!.....low-re version......................................................
      if (lowre) then                                                    !
!.....let's find alpha*                                                  !                                         
      alphast=(0.024+(densit*te(inp))/(6.*viscos*ed(inp)))   &           !           
             /(1.+(densit*te(inp))/(6.*viscos*ed(inp)))                  !
         vis(inp)=viscos  &                                              !
              +den(inp)*te(inp)/(ed(inp)+small)               &          !  
              *1./max(1./alphast, strain(inp)*f2_sst/(a1*ed(inp)))       !                                                   
!.....end of low-re version..............................................!
      end if

!...
      vis(inp)=urf(ivis)*vis(inp)+(1.-urf(ivis))*visold

      enddo
      enddo
      enddo
!====================================================
      endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(wilcox) then
!====================================================
!.....for wilcox k-omega model
!====================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm
      inp=lk(k)+li(i)+j
      visold=vis(inp)

      vis(inp)=viscos  &
              +den(inp)*te(inp)/(ed(inp)+small)

!     if (lowre) then
!.....low-re version......................................................
!.....let's find alpha*                                                  !
      alphast=(0.024+(densit*te(inp))/(6.*viscos*ed(inp)))   &           !
             /(1.+(densit*te(inp))/(6.*viscos*ed(inp)))                  !
      vis(inp)=viscos                                        &           !
              +alphast*den(inp)*te(inp)/(ed(inp)+small)                  !
!........................................................................!
!     endif
!...
      vis(inp)=urf(ivis)*vis(inp)+(1.-urf(ivis))*visold

      enddo
      enddo
      enddo
!====================================================
      endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(earsm_wj) then
!====================================================
!.....for earsm wj model 
!====================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm
      inp=lk(k)+li(i)+j
      visold=vis(inp)

!.....low-re version......................................................
      if (lowre) then                                                    !
!.....let's find alpha*                                                  !                                         
      alphast=(0.0249+(densit*te(inp))/(6.*viscos*ed(inp)))   &          !           
             /(1.+(densit*te(inp))/(6.*viscos*ed(inp)))                  !
      vis(inp)=viscos                                         &          !
            +alphast*cmueff(inp)/cmu*den(inp)*te(inp)/ed(inp)            !                                                      
!.....end of low-re version..............................................!
      end if
!.....high-re version.....................................................
      vis(inp) = viscos + (cmueff(inp)/cmu) * den(inp)*te(inp)/ed(inp)
!.....under-relaxation
      vis(inp)=urf(ivis)*vis(inp)+(1.-urf(ivis))*visold

      enddo
      enddo
      enddo
!====================================================
      endif


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(earsm_m) then
!====================================================
!.....for earsm wj model by menter et.al. (2009)
!====================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm
      inp=lk(k)+li(i)+j
      visold=vis(inp)

!.....low-re version......................................................
      if (lowre) then                                                    !
!.....let's find alpha*                                                  !                                         
      alphast=(0.0249+(densit*te(inp))/(6.*viscos*ed(inp)))   &          !           
             /(1.+(densit*te(inp))/(6.*viscos*ed(inp)))                  !
      vis(inp)=viscos                                         &          !
            +alphast                *den(inp)*te(inp)/ed(inp)            !                                                      
!.....end of low-re version..............................................!
      end if
!.....high-re version.....................................................
      vis(inp) = viscos +                     den(inp)*te(inp)/ed(inp)
!.....under-relaxation
      vis(inp)=urf(ivis)*vis(inp)+(1.-urf(ivis))*visold

      enddo
      enddo
      enddo
!====================================================
      endif

!      magvis =  volumeweightedaverage(vis)
!      write(66,*)'volume weighted average of effective visc.',magvis



! !#############################################
! !     [ obstacle modification: ]
! !#############################################
!       if(iobst.eq.0) return
!       do nsa=1,nobst ! obstacles---
!       ioss=ios(nsa)
!       ioes=ioe(nsa)
!       joss=jos(nsa)
!       joes=joe(nsa)
!       koss=kos(nsa)
!       koes=koe(nsa)

!       do k=koss,koes
!       do i=ioss,ioes
!       do j=joss,joes
!       inp=lk(k)+li(i)+j
!       vis(inp)=0.0d0
!       viscos=0.0d0
!       den(inp)=0.0d0
!       uu(inp)=0.0d0
!       vv(inp)=0.0d0
!       ww(inp)=0.0d0
!       uv(inp)=0.0d0
!       uw(inp)=0.0d0
!       vw(inp)=0.0d0
!       utt(inp)=0.0d0
!       vtt(inp)=0.0d0
!       wtt(inp)=0.0d0
!       end do ! j-loop
!       end do ! i-loop
!       end do ! k-loop

!       end do  ! obstacles----


      return
      end
