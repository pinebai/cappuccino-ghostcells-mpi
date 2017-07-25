!***********************************************************************
!
      subroutine intfac(idew,idns,idtb,lambda)
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

      integer, intent(in) :: idew,idns,idtb
      real(dp), dimension(nxyza), intent(inout) :: lambda

      integer i,j,k,inp,inn
      real(dp) :: xf,yf,zf

!                                          __   __
!.....lambda interpolation factor lambda = pj'/ pn 
      do k=2,nkm
      do i=2,nim
      do j=2,njm
      inp=lk(k)+li(i)+j
      inn=inp+idew

!.....Treba naci j' tacku na kojoj se sece ravan face-a i linija spajanja celijskih centara => find_intersection_point
      call find_intersection_point( &
!                            plane defined by face corner, bottom and south:
                             x(inp),y(inp),z(inp),&
                             x(inp-idns),y(inp-idns),z(inp-idns), &
                             x(inp-idtb),y(inp-idtb),z(inp-idtb), &
!                            line defined by cell center and neighbour center:
                             xc(inp),yc(inp),zc(inp), &
                             xc(inp+idew),yc(inp+idew),zc(inp+idew), &
!                            intersection point:
                             xf,yf,zf &
                             )

      lambda(inp) =   sqrt( (xf-xc(inp))**2 + (yf-yc(inp))**2 + (zf-zc(inp))**2 )  &
                  / ( sqrt( (xc(inn)-xc(inp))**2 + (yc(inn)-yc(inp))**2 + (zc(inn)-zc(inp)+small)**2 ) + small )
      enddo
      enddo
      enddo

      return
      end
