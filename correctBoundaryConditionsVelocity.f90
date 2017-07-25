      SUBROUTINE correctBoundaryConditionsVelocity
!  
!******************************************************************************
!
!     Updates values at boundaries with fixedGradient and zeroGradient B.C. and
!     periodicity after solution for velocity is obtained.
! 
!******************************************************************************
!
      use types
      use parameters
      use indexes
      use variables
      use bc

      implicit none

!
!     Local variables
!
      integer :: i, j, k, inp
      integer :: inbc


!.....W E S T
      do k=2,nkm
      do j=2,njm
      inp=lk(k)+li(2)+j
      inbc=(j-1)*nk+k
      lbhlp=lbw(inbc)
      if(lbhlp.eq.1 .or. lbhlp.eq.3) then
      ! Inlet or Symmetry boundary
      U(inp-nj) = U(inp)
      V(inp-nj) = V(inp)
      W(inp-nj) = W(inp)
      endif
      enddo
      enddo

!.....E A S T
      do k=2,nkm
      do j=2,njm
      inp=lk(k)+li(nim)+j
      inbc=(j-1)*nk+k
      lbhlp=lbe(inbc)
      if(lbhlp.eq.2 .or. lbhlp.eq.3) then
      ! Outlet or Symmetry boundary
      U(inp+nj) = U(inp)
      V(inp+nj) = V(inp)
      W(inp+nj) = W(inp)
      endif
      enddo
      enddo

!.....S O U T H 
      do k=2,nkm
      do i=2,nim
      inp=lk(k)+li(i)+2
      inbc=(i-1)*nk+k
      lbhlp=lbs(inbc)
      if(lbhlp.eq.2 .or. lbhlp.eq.3) then
      ! Outlet or Symmetry boundary
      U(inp-1) = U(inp)
      V(inp-1) = V(inp)
      W(inp-1) = W(inp)
      endif
      enddo
      enddo

!.....N O R T H
      do k=2,nkm
      do i=2,nim
      inp=lk(k)+li(i)+njm
      inbc=(i-1)*nk+k
      lbhlp=lbn(inbc)
      if(lbhlp.eq.2 .or. lbhlp.eq.3) then
      ! Outlet or Symmetry boundary
      U(inp+1) = U(inp)
      V(inp+1) = V(inp)
      W(inp+1) = W(inp)
      endif
      enddo
      enddo

!.....B O T T O M
      do i=2,nim
      do j=2,njm
      inp=lk(2)+li(i)+j
      inbc=(i-1)*nj+j
      lbhlp=lbb(inbc)
      if(lbhlp.eq.2 .or. lbhlp.eq.3) then
      ! Outlet or Symmetry boundary
      U(inp-nij) = U(inp)
      V(inp-nij) = V(inp)
      W(inp-nij) = W(inp)
      endif
      enddo
      enddo

!.....TOP
      do i=2,nim
      do j=2,njm
      inp=lk(nkm)+li(i)+j
      inbc=(i-1)*nj+j
      lbhlp=lbt(inbc)
      if(lbhlp.eq.2 .or. lbhlp.eq.3) then
      ! Outlet or Symmetry boundary
      U(inp+nij) = U(inp)
      V(inp+nij) = V(inp)
      W(inp+nij) = W(inp)
      endif
      enddo
      enddo

      RETURN
      END
