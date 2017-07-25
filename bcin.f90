!***********************************************************************
!
      subroutine bcin
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use geometry
      use variables
      use bc
      use inlet

      implicit none
!
!***********************************************************************
!
      integer :: i, j, k, inp, inbc, inh
      real(prec) :: tenom, ednom

      ! triangular hill:
      !real(prec), parameter :: blthick = 0.15 ! bl thickness for flow case 'b'
      !real(prec), parameter :: expp = 0.17    ! exponent p for flow case 'b'
      !real(prec), parameter :: umag = 5.774641_dp    ! exponent p for flow case 'b'

      ! ishihara hill:
      ! real(prec), parameter :: blthick = 0.2 
      ! real(prec), parameter :: expp = 0.135    
      ! real(prec), parameter :: umag = 5.5   

      ! dobric hill
      !real(prec) :: ust  

      ! za channel395
      real(prec) :: zz,perturb

!----------------------------------------------
!.....read inlet values and calc. mass fluxes
!----------------------------------------------
      flowin=0.
      xmonin=0.
!-------------------
!.....West
!-------------------
      do k=2,nkm
      do j=2,njm

      inp=lk(k)+li(2)+j

      inbc=(j-1)*nk+k
      lbhlp=lbw(inbc+jks)

      if(lbhlp.eq.1) then
!.....ako zelis da ti za inflow procita iz fajla inlet:
      ! read(7,*) u(inp),v(inp),w(inp),p(inp),te(inp),ed(inp),t(inp)

      !ed(inp)=ed(inp)/(cmu*te(inp)+small) ! inlet napravljen za k-epsilon

!.....ne treba da cita, vec je procitao
      !u(inp)=u_inl(k)
      !v(inp)=v_inl(k)
      !w(inp)=w_inl(k)
      !te(inp)=te_inl(k)
      !ed(inp)=ed_inl(k)

!.....ishihara hill case:
      !u(inp) = umag*(zc(inp)/blthick)**expp
      !v(inp) = small
      !w(inp) = small
      !te(inp) = tein
      !ed(inp) = edin!/(cmu*te(inp)) ! bez '/(cmu*te(inp))' faktora za k-eps
      !write(7,*) u(inp),v(inp),w(inp),p(inp),te(inp),ed(inp),t(inp)

!.....channel flow with random noise:
      zz=zc(inp)

      ! random number based fluctuation of mean profile            
      call init_random_seed()
      call random_number(perturb)
      
      perturb = 0.9+perturb/5. ! max perturbation is 10% of mean profile

      u(inp) = perturb*magubar*(1.25*(1.-zz**4))
      v(inp) = small
      w(inp) = small
      p(inp) = small

!.....Bolund hill case:
!      u(inp)=log((zc(inp)*120-zc(inbc)*120)/0.0003)
!      v(inp)=0.
!      w(inp)=0.
!      p(inp)=0.
!      te(inp)=0.928
      !lscale=0.4*0.1 ! l=0.4*d; d sam uzeo visinu brda 12 m skalirano
            !faktorom 1/120 jer je bila fora da visina bude 1 m a ne 120
      !ed(inp)=cmu75*te(inp)**1.5/lscale
!      ed(inp)=(0.4)**3/(cappa*zc(inp)*120) !u*=0.4; ed=u*^3/cappa*(z+z0))
     
!.....Triangular hill case:
      !te(inp)=max(small, &
      !        -3.71354692995576e-5*zc(inp)**7+0.0008508042*zc(inp)**6-0.0074501233*zc(inp)**5 &
      !        +0.0287202493*zc(inp)**4-0.0279210481*zc(inp)**3-0.1060775951*zc(inp)**2 &
      !        +0.1756108394*zc(inp)+0.2104104465)
      !ed(inp)=max(small, cmu75*te(inp)**1.5/(0.07*blthick)) !/(cmu*te(inp))

      vis(inp)=viscos
      if(lturb) then
        if(wilcox.or.sst.or.sas.or.earsm_wj.or.earsm_m) vis(inp)=vis(inp)+den(inp)*te(inp)/(ed(inp)+small)
        if(stdkeps.or.durbin.or.rng.or.realizable.or.lowre_lb) vis(inp)=vis(inp)+den(inp)*te(inp)**2*cmu/(ed(inp)+small)  
      endif

      inh=inp
      f1(inh)=den(inp)*(ar1x(inp)*u(inp)+ar1y(inp)*v(inp) &
                       +ar1z(inp)*w(inp))
      flowin=flowin+f1(inh)
      xmonin=xmonin+f1(inh)*dsqrt(u(inp)**2+v(inp)**2+w(inp)**2)
      endif
      enddo
      enddo
!-------------------
!.....East
!-------------------
      do k=2,nkm
      do j=2,njm

      inp=lk(k)+li(ni)+j

      inbc=(j-1)*nk+k
      lbhlp=lbe(inbc+jks)

      if(lbhlp.eq.1) then 
      read(7,*) u(inp),v(inp),w(inp),p(inp),te(inp),ed(inp),t(inp)

!.....Dobric hill case:
!.....log profile for dobric case;
!      inbc = li(ni)+j 
!      if ( (zc(inp)-zc(inbc)) .lt. 500. ) then
!        ust=0.345843
!        u(inp) = min(-(ust/cappa)*log((zc(inp)-zc(inbc))/zzero),-small)
!        te(inp) = ust**2/sqrt(cappa)*(1.-((zc(inp)-zc(inbc))/500.))**2
!      else
!        u(inp) = u(inp)
!        te(inp) = small
!      endif
!      ed(inp) = max(te(inp)**1.5/10.,small) !<-- epsilon eq.
!      !ed(inp) = ed(inp)/0.02973/te(inp) !<-- omega eq.
!      !ed(inp)=ed(inp)/(cmu*te(inp)+small) ! k-eps --> k-omega
!      if (j.eq.80) print*,k,u(inp),te(inp),ed(inp),(zc(inp)-zc(inbc))


      vis(inp)=viscos

      if(lturb) then
        if(wilcox.or.sst.or.sas.or.earsm_wj.or.earsm_m) vis(inp)=vis(inp)+den(inp)*te(inp)/(ed(inp)+small)
        if(stdkeps.or.durbin.or.rng.or.realizable.or.lowre_lb) vis(inp)=vis(inp)+den(inp)*te(inp)**2*cmu/(ed(inp)+small)  
      endif

      inh=inp-nj
      f1(inh)=den(inp)*(ar1x(inp)*u(inp)+ar1y(inp)*v(inp) &
                       +ar1z(inp)*w(inp))
      flowin=flowin-f1(inh)
      xmonin=xmonin-f1(inh)*dsqrt(u(inp)**2+v(inp)**2+w(inp)**2)
      endif
      enddo
      enddo
!-------------------
!.....south
!-------------------
      do k=2,nkm
      do i=2,nim
      inp=lk(k)+li(i)+1
      inbc=(i-1)*nk+k
      lbhlp=lbs(inbc+iks)
      if(lbhlp.eq.1) then 
      read(7,*) u(inp),v(inp),w(inp),p(inp),te(inp),ed(inp),t(inp)
      vis(inp)=viscos
      if(lturb)vis(inp)=vis(inp)+den(inp)*te(inp)**2*cmu/(ed(inp)+small)
      inh=inp
      f2(inh)=den(inp)*(ar1x(inp)*u(inp)+ar1y(inp)*v(inp) &
                       +ar1z(inp)*w(inp))
      flowin=flowin+f2(inh)
      xmonin=xmonin+f2(inh)*dsqrt(u(inp)**2+v(inp)**2+w(inp)**2)
      end if
      enddo
      enddo
!-------------------
!.....north
!-------------------
      do k=2,nkm
      do i=2,nim
      inp=lk(k)+li(i)+nj
      inbc=(i-1)*nk+k
      lbhlp=lbn(inbc+iks)
      if(lbhlp.eq.1) then 
      read(7,*) u(inp),v(inp),w(inp),p(inp),te(inp),ed(inp),t(inp)
      vis(inp)=viscos
      if(lturb)vis(inp)=vis(inp)+den(inp)*te(inp)**2*cmu/(ed(inp)+small)
      inh=inp-1
      f2(inh)=den(inp)*(ar1x(inp)*u(inp)+ar1y(inp)*v(inp) &
                       +ar1z(inp)*w(inp))
      flowin=flowin-f2(inh)
      xmonin=xmonin-f2(inh)*dsqrt(u(inp)**2+v(inp)**2+w(inp)**2)
      endif
      enddo
      enddo
!-------------------
!.....bottom
!-------------------
      do i=2,nim
      do j=2,njm
      inp=lk(1)+li(i)+j
      inbc=(i-1)*nj+j
      lbhlp=lbb(inbc+ijs)
      if(lbhlp.eq.1) then
      read(7,*) u(inp),v(inp),w(inp),p(inp),te(inp),ed(inp),t(inp)
      vis(inp)=viscos
      if(lturb)vis(inp)=vis(inp)+den(inp)*te(inp)**2*cmu/(ed(inp)+small)
      inh=inp
      f3(inh)=den(inp)*(ar1x(inp)*u(inp)+ar1y(inp)*v(inp) &
                       +ar1z(inp)*w(inp))
      flowin=flowin+f3(inh)
      xmonin=xmonin+f3(inh)*dsqrt(u(inp)**2+v(inp)**2+w(inp)**2)
      endif
      enddo
      enddo
!-------------------
!.....top
!-------------------
      do i=2,nim
      do j=2,njm
            
      inp=lk(nk)+li(i)+j

      inbc=(i-1)*nj+j
      lbhlp=lbt(inbc+ijs)

      if(lbhlp.eq.1) then
      read(7,*) u(inp),v(inp),w(inp),p(inp),te(inp),ed(inp),t(inp)
      vis(inp)=viscos
      if(lturb)vis(inp)=vis(inp)+den(inp)*te(inp)**2*cmu/(ed(inp)+small)
      inh=inp-nij
      f3(inh)=den(inp)*(ar1x(inp)*u(inp)+ar1y(inp)*v(inp) &
                       +ar1z(inp)*w(inp))
      flowin=flowin-f3(inh)
      xmonin=xmonin-f3(inh)*dsqrt(u(inp)**2+v(inp)**2+w(inp)**2)
      endif
      enddo
      enddo
!
!.....normalization factors for residual normalization
      if(.not.lout) then
      flowin = 1.
      xmonin = 1.
      tenom  = 1.
      ednom  = 1.
      end if

      ! write(66,'(a,es11.4)')'flowin: ',flowin
      ! write(66,'(a,es11.4)')'xmonin: ',xmonin

      do i=1,nphi
      resor(i) = 0.
      enddo

      snorin(iu)   = 1./(xmonin+small)
      snorin(iv)   = snorin(iu)
      snorin(iw)   = snorin(iu)
      snorin(ip)   = 1./(flowin+small)
      snorin(ite)  = 1./(flowin*tenom+small)
      snorin(ied)  = 1./(flowin*ednom+small)
      snorin(ien)  = 1.
      snorin(ivart)= 1.
      snorin(icon) = 1.
! !
! !###########################################
! !     [obstacle modifications: ]
! !###########################################
!       if(iobst.eq.0) return

!       do nsa=1,nobst
!       ioss=ios(nsa)
!       ioes=ioe(nsa)
!       joss=jos(nsa)
!       joes=joe(nsa)
!       koss=kos(nsa)
!       koes=koe(nsa)

!       do k=koss,koes
!       do i=ioss,ioes
!       do j=joss,joes
!       ijk=lk(k)+li(i)+j
!       u(ijk)=0.
!       v(ijk)=0.
!       w(ijk)=0.
!       t(ijk)=0.
!       te(ijk)=0.
!       ed(ijk)=0.
!       vis(ijk)=0.
!       visob(ijk)=0.
!       den(ijk)=0.
!       uu(ijk)=0.
!       vv(ijk)=0.
!       ww(ijk)=0.
!       uv(ijk)=0.
!       uw(ijk)=0.
!       vw(ijk)=0.
!       vart(ijk)=0.
!       utt(ijk)=0.
!       vtt(ijk)=0.
!       wtt(ijk)=0.
!       end do
!       end do
!       end do

!       end do  ! nsa=1,nobst

      return
      end
