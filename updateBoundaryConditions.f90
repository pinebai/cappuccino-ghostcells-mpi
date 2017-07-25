      SUBROUTINE updateBoundaryConditions(phi)
!  
!******************************************************************************
!
!     Updates values at boundaries with fixedGradient and zeroGradient B.C.
!     after solution for $\phi$ is obtained.
! 
!******************************************************************************
!
      use types
      use parameters
      use indexes
      use bc

      implicit none

      real(prec), dimension(nxyza) :: phi
!
!     Local variables
!
      integer :: i, j, k, inp
      integer :: inbc

!.....Modify matrix coefficients to reflect presence of Boundary Conditions in PDE problem.

!.....WEST
      do k=2,nkm
      do j=2,njm
      inp=lk(k)+li(2)+j
      inbc=(j-1)*nk+k
      lbhlp=lbw(inbc+jks)
      if(lbhlp.ne.4) then
      ! Zero gradient B.C.
      phi(inp-nj) = phi(inp)
      endif
      enddo
      enddo

!.....E A S T
      do k=2,nkm
      do j=2,njm
      inp=lk(k)+li(nim)+j
      inbc=(j-1)*nk+k
      lbhlp=lbe(inbc+jks)
      if(lbhlp.ne.4) then
      ! Zero gradient B.C.
      phi(inp+nj) = phi(inp)
      endif
      enddo
      enddo

!.....S O U T H 
      do k=2,nkm
      do i=2,nim
      inp=lk(k)+li(i)+2
      inbc=(i-1)*nk+k
      lbhlp=lbs(inbc+iks)
      if(lbhlp.ne.4) then
      ! Zero gradient B.C.
      phi(inp-1) = phi(inp)
      endif
      enddo
      enddo

!.....N O R T H
      do k=2,nkm
      do i=2,nim
      inp=lk(k)+li(i)+njm
      inbc=(i-1)*nk+k
      lbhlp=lbn(inbc+iks)
      if(lbhlp.ne.4) then
      ! Zero gradient B.C.
      phi(inp+1) = phi(inp)
      endif
      enddo
      enddo

!.....BOTTOM
      do i=2,nim
      do j=2,njm
      inp=lk(2)+li(i)+j
      inbc=(i-1)*nj+j
      lbhlp=lbb(inbc+ijs)
      if(lbhlp.eq.4) then
!.....Wall  - Dirichlet (fixed value) condition at boundary
      phi(inp-nij) = 0.0d0 ! Value at fixed value boundary FIXME
      endif
      enddo
      enddo

!.....TOP
      do i=2,nim
      do j=2,njm
      inp=lk(nkm)+li(i)+j
      inbc=(i-1)*nj+j
      lbhlp=lbt(inbc+ijs)
      if(lbhlp.ne.4) then
      ! Zero gradient B.C.
      phi(inp+nij) = phi(inp)
      endif
      if(lbhlp.eq.4) then
!.....Wall  - Dirichlet (fixed value) condition at boundary
      phi(inp+nij) = 0.0d0 ! Value at fixed value boundary FIXME
      endif
      enddo
      enddo

      RETURN
      END
