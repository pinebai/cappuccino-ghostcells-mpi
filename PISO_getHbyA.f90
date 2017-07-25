!***********************************************************************
!
      SUBROUTINE PISO_getHbyA
!
!***********************************************************************
!
!     Assemble A(U), and H(U) excluding pressure gradient
!     Update U according to U=1/A(U)*H(U)
!     Peric denotes such U by Utilde
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE GEOMETRY
      USE VARIABLES
      USE COEF
      USE BUOY
      USE TIME_MOD
      USE GRADIENTS
      USE HCOEF

      IMPLICIT NONE
!
!***********************************************************************
!

!     Local variables
      integer :: i, j, k, inp, ine,inw,inn,ins,int,inb
      real(prec) :: apotime, heat, &
                    sut, svt, swt

!.....Loop over inner cells
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

!.....INITIALIZE SOURCE TERMS
      sp(inp) = 0.0d0
      su(inp) = 0.0d0
      sv(inp) = 0.0d0
      sw(inp) = 0.0d0

!=====================================
!.....BUOYANCY SOURCE TERMS
!=====================================
      IF(LCAL(IEN).AND.LBUOY) THEN
!----------------------------------------------
!........[Boussinesq-ova aproximacija: ]
!----------------------------------------------
      heat=0.0d0
      if(boussinesq) then 
        heat=beta*densit*(t(inp)-tref)*vol(inp)
      else !if(boussinesq.eq.0)
        heat=(densit-den(inp))*vol(inp)
      endif
!----------------------------------------------
      su(inp)=su(inp)-gravx*heat
      sv(inp)=sv(inp)-gravy*heat
      sw(inp)=sw(inp)-gravz*heat
!=====================================
      ENDIF


!.....UNSTEADY TERM
      if(bdf) then
!    Three Level Implicit Time Integration Method:
!    in case that BTIME=0. --> Implicit Euler
      apotime=den(inp)*vol(inp)/timestep
      sut=apotime*((1+btime)*uo(inp)-0.5*btime*uoo(inp))
      svt=apotime*((1+btime)*vo(inp)-0.5*btime*voo(inp))
      swt=apotime*((1+btime)*wo(inp)-0.5*btime*woo(inp))

      su(inp)=su(inp)+sut
      sv(inp)=sv(inp)+svt
      sw(inp)=sw(inp)+swt

      endif

      end do
      end do
      end do


!.....REYNOLS STRESSES
!      CALL CALCSTRESS
!=====================================
!.....[ADDITIONAL TERMS: ]
!=====================================
!      IF(LTURB.AND.LASM) THEN
!      END IF  !!![ASM approach] 


!
!.....Assemble A(u) and H(u) and update velocities
!
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      ine=inp+nj
      inw=inp-nj
      inn=inp+1
      ins=inp-1
      int=inp+nij
      inb=inp-nij

      if(CN) then
!.....Crank-Nicolson time stepping source terms
      apotime=den(inp)*vol(inp)/timestep
      su(inp)=su(inp)+(he(inp)*uo(ine) + hw(inp)*uo(inw)+  &
                       hn(inp)*uo(inn) + hs(inp)*uo(ins)+  &
                       ht(inp)*uo(int) + hb(inp)*uo(inb))+ &
              (apotime-he(inp)-hw(inp)                     &
                      -hn(inp)-hs(inp)                     &
                      -ht(inp)-hb(inp))*uo(inp)               
!.....End of Crank-Nicolson time stepping source terms
!.....End of Crank-Nicolson stuff
      endif

!.....Assemble H(U):
      su(inp)=su(inp)+(he(inp)*u(ine)+hw(inp)*u(inw)+ &
                       hn(inp)*u(inn)+hs(inp)*u(ins)+ &
                       ht(inp)*u(int)+hb(inp)*u(inb))
      enddo
      enddo
      enddo  
  
!
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j
   
      ine=inp+nj
      inw=inp-nj
      inn=inp+1
      ins=inp-1
      int=inp+nij
      inb=inp-nij

      if(CN) then
!.....Crank-Nicolson time stepping source terms
      apotime=den(inp)*vol(inp)/timestep
      sv(inp)=sv(inp)+(he(inp)*vo(ine) + hw(inp)*vo(inw)+   &
                       hn(inp)*vo(inn) + hs(inp)*vo(ins)+   &
                       ht(inp)*vo(int) + hb(inp)*vo(inb))+  &
              (apotime-he(inp)-hw(inp)                      &
                      -hn(inp)-hs(inp)                      &
                      -ht(inp)-hb(inp))*vo(inp)
!.....End of Crank-Nicolson time stepping source terms
      endif

!.....Assemble H(U):
      sv(inp)=sv(inp)+(he(inp)*v(ine)+hw(inp)*v(inw)+ &
                       hn(inp)*v(inn)+hs(inp)*v(ins)+ &
                       ht(inp)*v(int)+hb(inp)*v(inb))
      enddo
      enddo
      enddo

!
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      ine=inp+nj
      inw=inp-nj
      inn=inp+1
      ins=inp-1
      int=inp+nij
      inb=inp-nij

      if(CN) then
!.....Crank-Nicolson time stepping source terms
      apotime=den(inp)*vol(inp)/timestep
      sw(inp)=sw(inp)+(he(inp)*wo(ine) + hw(inp)*wo(inw)+   &
                       hn(inp)*wo(inn) + hs(inp)*wo(ins)+   &
                       ht(inp)*wo(int) + hb(inp)*wo(inb))+  &
              (apotime-he(inp)-hw(inp)                      &
                      -hn(inp)-hs(inp)                      &
                      -ht(inp)-hb(inp))*wo(inp)
!.....End of Crank-Nicolson time stepping source terms
      endif

!.....Assemble H(U):
      sw(inp)=sw(inp)+(he(inp)*w(ine)+hw(inp)*w(inw)+ &
                       hn(inp)*w(inn)+hs(inp)*w(ins)+ &
                       ht(inp)*w(int)+hb(inp)*w(inb))
      enddo
      enddo
      enddo

!=====PISO UPDATE U, V, and W===========================================
!.....APU,APV,APW remain the same after being calculated during momentum predictor!
      do k=3,nkmm
        do i=3,nimm
          do j=3,njmm

            inp=lk(k)+li(i)+j

            u(inp) = su(inp) * apu(inp)
            v(inp) = sv(inp) * apv(inp)
            w(inp )= sw(inp) * apw(inp)
            
          enddo
        enddo
      enddo

!.....Velocity gradients (will be used for interpolation in peqn_flux function: 
      if (lstsq) then
!=====lstsq=============================================================  
      call grad_lsq_qr(u,gradu,2)
      call grad_lsq_qr(v,gradv,2)
      call grad_lsq_qr(w,gradw,2)

      elseif (gauss) then
!=====gauss=============================================================
      call grad_gauss(u,gradu(1,:),gradu(2,:),gradu(3,:))
      call grad_gauss(v,gradv(1,:),gradv(2,:),gradv(3,:))
      call grad_gauss(w,gradw(1,:),gradw(2,:),gradw(3,:))
      endif

      return
      end
