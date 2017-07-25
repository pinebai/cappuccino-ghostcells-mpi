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
                             xjp,yjp,zjp &
                             )
!
!***********************************************************************
! Find intersection point (pj={xjp,yjp,zjp}) of 
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
! xjp = x4 +(x5-x4)*t
! yjp = y4 +(y5-y4)*t
! zjp = z4 +(z5-z4)*t
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
!                             xjp,yjp,zjp &
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
      real(dp), intent(inout) :: xjp,yjp,zjp

      real(dp) :: t

     ! izracunato u matlabu simbolika
     t =-(x2*(y3*z4-y4*z3)-x1*(y3*z4-y4*z3)-x3*(y2*z4-y4*z2)+x1*(y2*z4-y4*z2)+x3*(y1*z4-y4*z1)-x2* &
         (y1*z4-y4*z1)+x4*(y2*z3-y3*z2)-x1*(y2*z3-y3*z2)-x4*(y1*z3-y3*z1)+x2*(y1*z3-y3*z1)+x4* &
         (y1*z2-y2*z1)-x3*(y1*z2-y2*z1)) &
        /(x2*(y3*(z5-z4)-(y5-y4)*z3)-x1*(y3*(z5-z4)-(y5-y4)*z3)-x3*(y2*(z5-z4)-(y5-y4)*z2)+x1* &
         (y2*(z5-z4)-(y5-y4)*z2)+x3*(y1*(z5-z4)-(y5-y4)*z1)-x2*(y1*(z5-z4)-(y5-y4)*z1)+(x5-x4)* &
         (y2*z3-y3*z2)-(x5-x4)*(y1*z3-y3*z1)+(x5-x4)*(y1*z2-y2*z1))

      xjp = x4 +(x5-x4)*t
      yjp = y4 +(y5-y4)*t
      zjp = z4 +(z5-z4)*t

      return
      end
