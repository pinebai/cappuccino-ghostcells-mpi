      subroutine grad_gauss(u,dudx,dudy,dudz)
!
!***********************************************************************
!
!     Calculates cell centered gradient using Gauss theorem
!     PARAMETERS
!     U - field, the gradient of which we are looking for
!     DUDX,DUDY,DUDZ - arrays where the gradient components are stored
!
!     Gauss gradient rule:
!     ------->                                 ->
!     grad(u) = 1/vol * sum_{i=1}^{i=nf} (u)_f*Sf
!     where:
!     grad(u) - cell centered gradient vector
!     (u)_f   - face interpolated value of scalar u
!     vol     - cell volume
!     Sf      - cell face area vector
!     nf      - number of facec in a cell
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
      real(prec), dimension(nxyza), intent(in) :: u
      real(prec), dimension(nxyza), intent(inout) :: dudx,dudy,dudz

      integer :: i,j,k,inp,ine,inn,int
      real(prec) :: ue,un,ut,volr

      dudx = 0.0d0
      dudy = 0.0d0
      dudz = 0.0d0

      ! Inner faces loop
      do k=2,nkmm
      do i=2,nimm
      do j=2,njmm

      inp=lk(k)+li(i)+j

      ine=inp+nj
      inn=inp+1
      int=inp+nij


      ue = u(ine)*fx(inp)+u(inp)*(1.0d0-fx(inp))
      un = u(inn)*fy(inp)+u(inp)*(1.0d0-fy(inp))
      ut = u(int)*fz(inp)+u(inp)*(1.0d0-fz(inp))

      dudx(inp) = dudx(inp) + (ue*ar1x(inp)+un*ar2x(inp)+ut*ar3x(inp) )

      dudx(ine) = dudx(ine) - ( ue*ar1x(inp) )

      dudx(inn) = dudx(inn) - ( un*ar2x(inp) )

      dudx(int) = dudx(int) - ( ut*ar3x(inp) )

 
      dudy(inp) = dudy(inp) + ( ue*ar1y(inp)+un*ar2y(inp)+ut*ar3y(inp) )

      dudy(ine) = dudy(ine) - ( ue*ar1y(inp) )

      dudy(inn) = dudy(inn) - ( un*ar2y(inp) )

      dudy(int) = dudy(int) - ( ut*ar3y(inp) )


      dudz(inp) = dudz(inp) + ( ue*ar1z(inp)+un*ar2z(inp)+ut*ar3z(inp) )
 
      dudz(ine) = dudz(ine) - ( ue*ar1z(inp) )

      dudz(inn) = dudz(inn) - ( un*ar2z(inp) )

      dudz(int) = dudz(int) - ( ut*ar3z(inp) )

      end do
      end do 
      end do 

      ! Inner cells loop
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      volr = 1.0d0/vol(inp)

      dudx(inp) = dudx(inp)*volr
      dudy(inp) = dudy(inp)*volr
      dudz(inp) = dudz(inp)*volr

      end do
      end do 
      end do 

      return
      end
