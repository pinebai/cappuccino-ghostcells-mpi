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

      integer :: i,j,k,inp,ine,inn,int,inw,ins,inb
      real(prec) :: ue,uw,un,us,ut,ub,volr

      ! Inner cell loop
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      ine=inp+nj
      inn=inp+1
      int=inp+nij
      inw=inp-nj
      ins=inp-1
      inb=inp-nij

      volr = 1.0d0/vol(inp)

      ue = u(ine)*fx(inp)+u(inp)*(1.0d0-fx(inp))
      uw = u(inp)*fx(inw)+u(inw)*(1.0d0-fx(inw))
      un = u(inn)*fy(inp)+u(inp)*(1.0d0-fy(inp))
      us = u(inp)*fy(ins)+u(ins)*(1.0d0-fy(ins))
      ut = u(int)*fz(inp)+u(inp)*(1.0d0-fz(inp))
      ub = u(inp)*fz(inb)+u(inb)*(1.0d0-fz(inb))


      dudx(inp) = (ue*ar1x(inp)+un*ar2x(inp)+ut*ar3x(inp) &
                  -uw*ar1x(inw)-us*ar2x(ins)-ub*ar3x(inb))*volr
      dudy(inp) = (ue*ar1y(inp)+un*ar2y(inp)+ut*ar3y(inp) &
                  -uw*ar1y(inw)-us*ar2y(ins)-ub*ar3y(inb))*volr
      dudz(inp) = (ue*ar1z(inp)+un*ar2z(inp)+ut*ar3z(inp) &
                  -uw*ar1z(inw)-us*ar2z(ins)-ub*ar3z(inb))*volr

      end do
      end do 
      end do 

      return
      end
