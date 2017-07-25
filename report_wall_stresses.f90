!***********************************************************************
!
      subroutine report_wall_stresses
!
!***********************************************************************
!
! Calculate and print wall shear stresses and nondimensional wall 
! distance y+.
! We are interested in interfaces of inner and ghost cells.
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use bc
      use boundc
      use wall

      implicit none
!
!***********************************************************************
!
      integer :: i, j, k, inp, inbc

      do inp=1,nyz
      gentw(inp) = 0.0_dp
      gente(inp) = 0.0_dp
      suedw(inp) = 0.0_dp
      suede(inp) = 0.0_dp
      enddo
      do inp=1,nxz
      gents(inp) = 0.0_dp
      gentn(inp) = 0.0_dp
      sueds(inp) = 0.0_dp
      suedn(inp) = 0.0_dp
      enddo
      do inp=1,nxy
      gentb(inp) = 0.0_dp
      gentt(inp) = 0.0_dp
      suedb(inp) = 0.0_dp
      suedt(inp) = 0.0_dp
      enddo

      if(lcal(iu)) then
!
!.....west
      do k=3,nkmm
      do j=3,njmm
      inp=lk(k)+li(3)+j
      inbc=(j-1)*nk+k
      lbhlp=lbw(inbc+jks)
      gentw(inbc) = 0.0_dp
      suedw(inbc) = 0.0_dp
      if(lbhlp.eq.4.) then
!.....wall
      ! call wallbc(inp,nj,1,nij)
      gentw(inbc)=ypl
      suedw(inbc)=tau
      endif
      end do
      end do

!.....east
      do k=3,nkmm
      do j=3,njmm
      inp=lk(k)+li(nimm)+j
      inbc=(j-1)*nk+k
      lbhlp=lbe(inbc+jks)
      gente(inbc) = 0.0_dp
      suede(inbc) = 0.0_dp
      if(lbhlp.eq.4) then
!.....wall
      ! call wallbc(inp,-nj,1,nij)
      gente(inbc)=ypl
      suede(inbc)=tau
      endif
      end do
      end do

!.....south
      do k=3,nkmm
      do i=3,nimm
      inp=lk(k)+li(i)+3
      inbc=(i-1)*nk+k
      lbhlp=lbs(inbc+iks)
      gents(inbc) = 0.0_dp
      sueds(inbc) = 0.0_dp
      if(lbhlp.eq.4) then
!.....wall
      ! call wallbc(inp,1,nij,nj)
      gents(inbc)=ypl
      sueds(inbc)=tau
      endif
      end do
      end do

!.....north
      do k=3,nkmm
      do i=3,nimm
      inp=lk(k)+li(i)+njmm
      inbc=(i-1)*nk+k
      lbhlp=lbn(inbc+iks)
      gentn(inbc) = 0.0_dp
      suedn(inbc) = 0.0_dp
      if(lbhlp.eq.4) then
!.....wall
      ! call wallbc(inp,-1,nij,nj)
      gentn(inbc)=ypl
      suedn(inbc)=tau
      endif
      end do
      end do

!.....botom
      do i=3,nimm
      do j=3,njmm
      inp=lk(3)+li(i)+j
      inbc=(i-1)*nj+j
      lbhlp=lbb(inbc+ijs)
      gentb(inbc) = 0.0_dp
      suedb(inbc) = 0.0_dp
      if(lbhlp.eq.4) then
!.....wall
      ! call wallbc(inp,nij,nj,1)
      gentb(inbc)=ypl
      suedb(inbc)=tau
      endif
      end do
      end do

!.....top
      do i=3,nimm
      do j=3,njmm
      inp=lk(nkmm)+li(i)+j
      inbc=(i-1)*nj+j
      lbhlp=lbt(inbc+ijs)
      gentt(inbc) = 0.0_dp
      suedt(inbc) = 0.0_dp
      if(lbhlp.eq.4) then
!.....wall
      ! call wallbc(inp,-nij,nj,1)
      gentt(inbc)=ypl
      suedt(inbc)=tau
      endif
      end do
      end do

!.....print
!      call print2(nj,nk,'j','k',gentw,'ypl-west  ')
!      call print2(nj,nk,'j','k',gente,'ypl-east  ')
!      call print2(ni,nk,'i','k',gents,'ypl-south ')
!      call print2(ni,nk,'i','k',gentn,'ypl-north ')
      call print2(ni,nj,'i','j',gentb,'ypl-bottom')
      call print2(ni,nj,'i','j',gentt,'ypl-top   ')
!.....tau
!      call print2(nj,nk,'j','k',suedw,'tau-west  ')
!      call print2(nj,nk,'j','k',suede,'tau-east  ')
!      call print2(ni,nk,'i','k',sueds,'tau-south ')
!      call print2(ni,nk,'i','k',suedn,'tau-north ')
      call print2(ni,nj,'i','j',suedb,'tau-bottom')
      call print2(ni,nj,'i','j',suedt,'tau-top   ')
      endif
      return
      end
