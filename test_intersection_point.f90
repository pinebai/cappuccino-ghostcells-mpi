 program intersection

 implicit none

 integer, parameter :: dp = kind(1.0d0)
 real(dp) :: xjp,yjp,zjp

! Test 1. - rectilinear case
 call find_intersection_point( &
!!                           plane defined by face corner, bottom and south:
                             1._dp,1._dp,1._dp,&
                             1._dp,0._dp,1._dp, &
                             1._dp,1._dp,0._dp, &
!                            line defined by cell center and neighbour center:
                             0.5_dp,0.5_dp,0.5_dp, &
                             1.5_dp,0.5_dp,0.5_dp, &
!                            intersection point:
                             xjp,yjp,zjp &
                             )
 print*,'Test 1:'
 print*, 'Exact values: 1.0, 0.5, 0.5'
 print*, 'Calculated values: ',xjp,yjp,zjp

! Test 2. - cells with intersection point offset - same plane as before but cells are skewed but different three points determine it.
 call find_intersection_point( &
!!                           plane defined by face corner, bottom and south:
                             1._dp,1.25_dp,1._dp,&
                             1._dp,0.25_dp,1._dp, &
                             1._dp,1.25_dp,0._dp, &
!                            line defined by cell center and neighbour center:
                             0.5_dp,0.5_dp,0.5_dp, &
                             1.5_dp,0.5_dp,0.5_dp, &
!                            intersection point:
                             xjp,yjp,zjp &
                             )
 print*,'Test 2:'
 print*, 'Exact values: 1.0, 0.5, 0.5'
 print*, 'Calculated values: ',xjp,yjp,zjp


 end program

!***********************************************************************
!
      subroutine find_intersection_point( &
!!                           plane defined by face corner, bottom and south:
                             x1,y1,z1,&
                             x2,y2,z2, &
                             x3,y3,z3, &
!                            line defined by cell center and neighbour center:
                             x4,y4,z4, &
                             x5,y5,z5, &
!                            intersection point:
                             xj,yj,zj &
                             )
!
!***********************************************************************
! Find intersection point (pj={xj,yj,zj}) of 
! plane (defined by points p1={x1,y1,z1}, p2={x2,y2,z2} and p3={x3,y3,z3}),
! and line (defined by points p4={x4,y4,z4} and p5={x5,y5,z5}).
!
!
!       |1  1  1  1 |     |1  1  1  0    |
! t = - |x1 x2 x3 x4|  /  |x1 x2 x3 x5-x4|  (mind the minus sign!)
!       |y1 y2 y3 y4| /   |y1 y2 y3 y5-y4|
!       |z1 z2 z3 z4|     |z1 z2 z3 z5-z4|
!
! And intersection point is given by:
! xj = x4 +(x5-x4)*t
! yj = y4 +(y5-y4)*t
! zj = z4 +(z5-z4)*t
!
!
! Nikola Mirkov @2016
!
! example usage: 
! call find_intersection_point( &
!!                            plane defined by face corner, bottom and south:
!                             x(inp),y(inp),z(inp),&
!                             x(inp-idns),y(inp-idns),z(inp-idns), &
!                             x(inp-idtb),y(inp-idtb),z(inp-idtb), &
!!                            line defined by cell center and neighbour center:
!                             xc(inp),yc(inp),zc(inp), &
!                             xc(inp+idew),yc(inp+idew),zc(inp+idew), &
!!                            intersection point:
!                             xj,yj,zj &
!                             )
!***********************************************************************
      use types

      implicit none 
!
!***********************************************************************
!

      real(dp), intent(in) :: x1,y1,z1,&
                              x2,y2,z2, &
                              x3,y3,z3, &
                              x4,y4,z4, &
                              x5,y5,z5
      real(dp), intent(inout) :: xj,yj,zj

      real(dp) :: t

     ! izracunato u matlabu simbolika
     t =-(x2*(y3*z4-y4*z3)-x1*(y3*z4-y4*z3)-x3*(y2*z4-y4*z2)+x1*(y2*z4-y4*z2)+x3*(y1*z4-y4*z1)-x2* &
         (y1*z4-y4*z1)+x4*(y2*z3-y3*z2)-x1*(y2*z3-y3*z2)-x4*(y1*z3-y3*z1)+x2*(y1*z3-y3*z1)+x4* &
         (y1*z2-y2*z1)-x3*(y1*z2-y2*z1)) &
        /(x2*(y3*(z5-z4)-(y5-y4)*z3)-x1*(y3*(z5-z4)-(y5-y4)*z3)-x3*(y2*(z5-z4)-(y5-y4)*z2)+x1* &
         (y2*(z5-z4)-(y5-y4)*z2)+x3*(y1*(z5-z4)-(y5-y4)*z1)-x2*(y1*(z5-z4)-(y5-y4)*z1)+(x5-x4)* &
         (y2*z3-y3*z2)-(x5-x4)*(y1*z3-y3*z1)+(x5-x4)*(y1*z2-y2*z1))

      xj = x4 +(x5-x4)*t
      yj = y4 +(y5-y4)*t
      zj = z4 +(z5-z4)*t

      return
      end
