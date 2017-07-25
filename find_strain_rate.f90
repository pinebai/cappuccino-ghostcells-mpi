      subroutine find_strain_rate
!
!***********************************************************************
!
      use types
      use parameters
      use indexes    
      use geometry 
      use coef
      use variables
      use gradients

      implicit none
!
!***********************************************************************
!
      integer :: i,j,k,inp
      real(prec) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz ! gradients
      real(prec) :: s11,s12,s13,s21,s22,s23,s31,s32,s33,w12,w13,w23

!.....velocity gradients: 
      if (lstsq) then 
        call grad_lsq_qr(u,gradu,2)
        call grad_lsq_qr(v,gradv,2)
        call grad_lsq_qr(w,gradw,2)
      elseif (gauss) then
        call grad_gauss(u,gradu(1,:),gradu(2,:),gradu(3,:))
        call grad_gauss(v,gradv(1,:),gradv(2,:),gradv(3,:))
        call grad_gauss(w,gradw(1,:),gradw(2,:),gradw(3,:))
      endif

      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

!==================================================================
      dudx = gradu(1,inp)
      dudy = gradu(2,inp)
      dudz = gradu(3,inp)
      
      dvdx = gradv(1,inp)
      dvdy = gradv(2,inp)
      dvdz = gradv(3,inp)

      dwdx = gradw(1,inp)
      dwdy = gradw(2,inp)
      dwdz = gradw(3,inp)

!.....find strain rate tensor
!     [s_ij]: |s_ij|=sqrt[2s_ij s_ij]
      s11=dudx
      s12=0.5*(dudy+dvdx)
      s13=0.5*(dudz+dwdx)
      s22=dvdy
      s23=0.5*(dvdz+dwdy) 
      s33=dwdz
      s21=s12  ! symmetry
      s31=s13
      s32=s23

!.....find antisymmetric part of velocity gradient tensor
!     [om_ij]: |om_ij|=sqrt[2 om_ij om_ij] 
      w12=0.5*(dudy - dvdx)
      w13=0.5*(dudz - dwdx)
      w23=0.5*(dvdz - dwdy)

!==================================================================

!.....find strain rate s = sqrt (2*sij*sij)
      strain(inp)=dsqrt(2*(s11**2+s22**2+s33**2 + 2*(s12**2+s13**2+s23**2)))

!.....find scalar invariant of antisymmetric part of velocity gradient tensor
      vorticity(inp) = dsqrt(w12**2 + w23**2 + w13**2)

      enddo
      enddo
      enddo

      end subroutine
