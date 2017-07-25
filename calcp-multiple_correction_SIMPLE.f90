!***********************************************************************
!
      subroutine calcp
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use geometry
      use coef
      use variables
      use buoy
      use time_mod
      use gradients
      use fieldmanipulation

      implicit none
!
!***********************************************************************
!
      integer, parameter :: isol=5
!
      integer :: i, j, k, inp, &
                 ine, inn, int, &
                 idew, idns, idtb
      real(prec) :: ppref,fmcor
      real(prec) :: reltol
      ! real(prec) :: magubarstar, ruaw, gragpplus, flowdirection ! za periodic channel

reltol = sor(ip)

! include 'const_massflux_correction.f90'

do icorr=1,ncorr

!.....tentative velocity gradients: 
      if (lstsq) then 
        call grad_lsq_qr(u,gradu,2)
        call grad_lsq_qr(v,gradv,2)
        call grad_lsq_qr(w,gradw,2)
      elseif (gauss) then
        call grad_gauss(u,gradu(1,:),gradu(2,:),gradu(3,:))
        call grad_gauss(v,gradv(1,:),gradv(2,:),gradv(3,:))
        call grad_gauss(w,gradw(1,:),gradw(2,:),gradw(3,:))
      endif

!
!.....assemble off diagonal entries of system matrix
      call fluxmass(nj,1,nij, &
                    ar1x,ar1y,ar1z, &
                    fx,ae,aw,f1)
!.....north cell - face
      call fluxmass(1,nij,nj, &
                    ar2x,ar2y,ar2z, &
                    fy,an,as,f2)
!.....top   cell - face
      call fluxmass(nij,nj,1, &
                    ar3x,ar3y,ar3z, &
                    fz,at,ab,f3)

      if(.not.const_mflux) call outbc
!
!.....assemble the main diagonal entry and rhs vector stored in su array
      idew=nj
      idns=1
      idtb=nij

      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      su(inp)=-f1(inp)+f1(inp-idew)-f2(inp)+f2(inp-idns)-f3(inp)+f3(inp-idtb)
      ap(inp)=ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp)

      enddo
      enddo
      enddo

!=====multiple pressure corrections=====================================
      do ipcorr=1,npcor
!
!.....initialize pressure correction
      pp=0.0d0


      ! if(icorr.ne.ncorr .or. ipcorr.ne.npcor) then
      ! if(icorr.ne.ncorr) then
      !   sor(ip) = 0.1
      ! else
      !   sor(ip) = reltol
      ! endif

!
!.....solving pressure correction equation
      if(isol.eq.1) call sipsol(pp,ip)
      if(isol.eq.2) call cgstab(pp,ip)
      if(isol.eq.3) call iccg(pp,ip)
      if(isol.eq.4) call pcgsip(pp,ip)
      if(isol.eq.5) call cgstab_sip(pp,ip)
         

!.....second step *** corrector stage

!
!.....reference pressure correction - p'
      ppref=pp(ijkpr)
!
!.....calculate pressure-correction gradients
      if (lstsq) then   
        call grad_lsq_qr(pp,gradp,2)
      elseif (gauss) then
        call grad_gauss(pp,gradp(1,:),gradp(2,:),gradp(3,:))
      endif

!.......the source term for the non-orthogonal corrector 
!.......also the secondary mass flux correction if ipcorr.ne.npcor.........
        if(ipcorr.ne.npcor) then                                          !
                                                                          !
        do inp=icst,icen                                                  !
        su(inp) = 0.0d0                                                   !
        enddo                                                             !
                                                                          !
        idew=nj                                                           !
        idns=1                                                            !
        idtb=nij                                                          !
                                                                          !
        do k=3,nkmm                                                       !
        do i=3,nimm                                                       !
        do j=3,njmm                                                       !
                                                                          !          
        inp=lk(k)+li(i)+j                                                 !
                                                                          !
        ine=inp+idew                                                      !
        inn=inp+idns                                                      !
        int=inp+idtb                                                      !
                                                                          !
                                                                          !
!.......correct mass fluxes at inner cv-faces for second corr.            !
!.......east cell - face                                                  !
        call fluxmc(inp,nj,1,nij,fx,fmcor)                                !
        f1(inp) = f1(inp)-fmcor                                           !
        su(inp) = su(inp)+fmcor                                           !
        su(ine) = su(ine)-fmcor                                           !
!.......north cell - face                                                 !
        call fluxmc(inp,1,nij,nj,fy,fmcor)                                !
        f2(inp) = f2(inp)-fmcor                                           !
        su(inp) = su(inp)+fmcor                                           !
        su(inn) = su(inn)-fmcor                                           !
!.......top   cell - face                                                 !
        call fluxmc(inp,nij,nj,1,fz,fmcor)                                !
        f3(inp) = f3(inp)-fmcor                                           !
        su(inp) = su(inp)+fmcor                                           !
        su(int) = su(int)-fmcor                                           !
                                                                          !
        enddo                                                             !
        enddo                                                             !
        enddo                                                             !
                                                                          !
        endif                                                             !
!.........................................................................!


      ! if(ipcorr.eq.npcor) then
!
!.....loop over inner cv's
      idew=nj
      idns=1
      idtb=nij

      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      ine=inp+idew
      inn=inp+idns
      int=inp+idtb
!
!.....correct mass fluxes at inner cv-faces only (only inner flux)
!
      f1(inp)=f1(inp)-ae(inp)*(pp(ine)-pp(inp)) 
      f2(inp)=f2(inp)-an(inp)*(pp(inn)-pp(inp)) 
      f3(inp)=f3(inp)-at(inp)*(pp(int)-pp(inp))

      enddo
      enddo
      enddo   

!
!.....correct velocities and pressures
!
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      u(inp)=u(inp)-vol(inp)*apu(inp)*gradp(1,inp) 
      v(inp)=v(inp)-vol(inp)*apv(inp)*gradp(2,inp) 
      w(inp)=w(inp)-vol(inp)*apw(inp)*gradp(3,inp) 
      p(inp)=p(inp)+urf(ip)*(pp(inp)-ppref)

      ! p(inp)=p(inp)+urf(ip)*(pp(inp)-p(inp))

      !su(inp) = 0.0d0 ! clean values of rhs vector en passant

      enddo
      enddo
      enddo   

      ! endif

!
!=====end: multiple pressure corrections loop==============================
      enddo


      include 'continuityErrors.h'

!.....explicit correction of boundary conditions 
      call correctBoundaryConditions

      ! include 'const_massflux_correction.f90'

      enddo

!
!.....update pressure at boundaries
      call bpres(p)

      pp = p

      return
      end
