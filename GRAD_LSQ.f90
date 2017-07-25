      subroutine grad_lsq(fi,gradfi,istage)
!
!***********************************************************************
!
!      Purpose:
!      Calculates cell-centered gradients using UNWEIGHTED Least-Squares approach.
!
!      Description:
!      Approach taken from PhD thesis of Bojan Niceno, TU Delft, 2000.
!
!      FI - field variable which gradient we look for.
!      GradFi - cell centered gradient - a three component gradient vector.
!      ISTAGE - integer. If ISTAGE=1 calculates and stores only geometrical
!      parameters - a system matrix for least square problem at every cell. 
!      Usually it is called with ISTAGE=1 at the beggining of simulation.
!      If 2 it doesn't calculate system matrix, just RHS and solves system.
!
!      D - least-squares matrix, calculated once at the beggining.
!      NXYZA - Integer the length of FI array in the calling routine.
!      NX - Length of LI
!      NZ - Length of LK
!      NKMM,NIMM,NJMM - number of inner cells in Z,X and Y direction respectively.
!      LK - How many cells is there below K-th floor.
!      LI - How many cells is there left for I-th vertical cross-section.
!      IDEW,IDNS,IDTB - Integer that can bring you to E,W,N,S,T,B cells from current cell P.
!      XC,YC,ZC - coordinates of cell centers    
!
!***********************************************************************
!
      use parameters
      use indexes
      use geometry
      use gradients
      use matrix_module

      implicit none

      integer, intent(in) :: istage
      real(prec),dimension(nxyza), intent(in)   :: fi
      real(prec),dimension(3,nxyza), intent(inout) :: gradfi

!
!    Locals
!
      integer :: i,j,k,inp
      integer :: ine, inw, inn, ins, int, inb           ! integer indexes of neighbour cells
      real(prec) :: d11,d12,d13,d21,d22,d23,d31,d32,d33 ! system matrix elements
      real(prec) :: b1,b2,b3                            ! rhs vector components
      real(prec) :: tmp
   

      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

!     indexes of east, south, north, bottom, and top cell
      ine = inp + nj
      inw = inp - nj
      inn = inp + 1
      ins = inp - 1 
      int = inp + nij
      inb = inp - nij
 

      if(istage.eq.1) then

!     coefficient matrix - should be calculated only once 
      d(1,inp) = (xc(ine)-xc(inp))**2+(xc(inw)-xc(inp))**2  &
                +(xc(inn)-xc(inp))**2+(xc(ins)-xc(inp))**2  &
                +(xc(int)-xc(inp))**2+(xc(inb)-xc(inp))**2

      d(4,inp) = (yc(ine)-yc(inp))**2+(yc(inw)-yc(inp))**2  &
                +(yc(inn)-yc(inp))**2+(yc(ins)-yc(inp))**2  &
                +(yc(int)-yc(inp))**2+(yc(inb)-yc(inp))**2

      d(6,inp) = (zc(ine)-zc(inp))**2+(zc(inw)-zc(inp))**2  &
                +(zc(inn)-zc(inp))**2+(zc(ins)-zc(inp))**2  &
                +(zc(int)-zc(inp))**2+(zc(inb)-zc(inp))**2

      d(2,inp) = (xc(ine)-xc(inp))*(yc(ine)-yc(inp)) + (xc(inw)-xc(inp))*(yc(inw)-yc(inp))  &
                +(xc(inn)-xc(inp))*(yc(inn)-yc(inp)) + (xc(ins)-xc(inp))*(yc(ins)-yc(inp))  &
                +(xc(int)-xc(inp))*(yc(int)-yc(inp)) + (xc(inb)-xc(inp))*(yc(inb)-yc(inp))

      d(3,inp) = (xc(ine)-xc(inp))*(zc(ine)-zc(inp)) + (xc(inw)-xc(inp))*(zc(inw)-zc(inp))  &
                +(xc(inn)-xc(inp))*(zc(inn)-zc(inp)) + (xc(ins)-xc(inp))*(zc(ins)-zc(inp))  &
                +(xc(int)-xc(inp))*(zc(int)-zc(inp)) + (xc(inb)-xc(inp))*(zc(inb)-zc(inp))
 
      d(5,inp) = (yc(ine)-yc(inp))*(zc(ine)-zc(inp)) + (yc(inw)-yc(inp))*(zc(inw)-zc(inp))  &
                +(yc(inn)-yc(inp))*(zc(inn)-zc(inp)) + (yc(ins)-yc(inp))*(zc(ins)-zc(inp))  &
                +(yc(int)-yc(inp))*(zc(int)-zc(inp)) + (yc(inb)-yc(inp))*(zc(inb)-zc(inp))

      elseif(istage.eq.2) then

!     copy from coefficient matrix 
      d11 = d(1,inp)
      d12 = d(2,inp)
      d13 = d(3,inp)

      d22 = d(4,inp)
      d23 = d(5,inp)
      d33 = d(6,inp)

!     symmetric part
      d21 = d12
      d31 = d13
      d32 = d23

!     rhs vector
      b1 = (fi(ine)-fi(inp))*(xc(ine)-xc(inp)) + (fi(inw)-fi(inp))*(xc(inw)-xc(inp)) & 
         + (fi(inn)-fi(inp))*(xc(inn)-xc(inp)) + (fi(ins)-fi(inp))*(xc(ins)-xc(inp)) & 
         + (fi(int)-fi(inp))*(xc(int)-xc(inp)) + (fi(inb)-fi(inp))*(xc(inb)-xc(inp))  

      b2 = (fi(ine)-fi(inp))*(yc(ine)-yc(inp)) + (fi(inw)-fi(inp))*(yc(inw)-yc(inp)) & 
         + (fi(inn)-fi(inp))*(yc(inn)-yc(inp)) + (fi(ins)-fi(inp))*(yc(ins)-yc(inp)) & 
         + (fi(int)-fi(inp))*(yc(int)-yc(inp)) + (fi(inb)-fi(inp))*(yc(inb)-yc(inp))  

      b3 = (fi(ine)-fi(inp))*(zc(ine)-zc(inp)) + (fi(inw)-fi(inp))*(zc(inw)-zc(inp)) & 
         + (fi(inn)-fi(inp))*(zc(inn)-zc(inp)) + (fi(ins)-fi(inp))*(zc(ins)-zc(inp)) & 
         + (fi(int)-fi(inp))*(zc(int)-zc(inp)) + (fi(inb)-fi(inp))*(zc(inb)-zc(inp))  
 
!     denominator used troughout
      tmp = 1./(d11*d22*d33 - d11*d23*d32 - d12*d21*d33 + d12*d23*d31 + d13*d21*d32 - d13*d22*d31)

!     solve system 
      gradFi(1,inp) = ( (b1*(d22*d33 - d23*d32)) - (b2*(d21*d33 - d23*d31)) + (b3*(d21*d32 - d22*d31)) ) * tmp
      gradFi(2,inp) = ( (b2*(d11*d33 - d13*d31)) - (b1*(d12*d33 - d13*d32)) - (b3*(d11*d32 - d12*d31)) ) * tmp
      gradFi(3,inp) = ( (b1*(d12*d23 - d13*d22)) - (b2*(d11*d23 - d13*d21)) + (b3*(d11*d22 - d12*d21)) ) * tmp

      endif

      enddo
      enddo
      enddo

      return
      end
