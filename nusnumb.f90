!***********************************************************************
!
      subroutine nusnumb
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use geometry
      use variables
      use nusselt

      implicit none
!
!***********************************************************************
!

      integer :: j, k, ijk, jk, lkk
      real(prec) :: dx, dy, dz
      real(prec) :: bnusoveralle, bnusoverallw
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     [side heated cavities: west-east thermally active]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=2,njm
      bnusmeanw(j)=0.
      bnusmeane(j)=0.
      do k=2,nkm
      lkk=lk(k)
      jk=(k-1)*nj+j
      ijk=lkk+li(2)+j
      dx=0.5*(x(ijk)-x(ijk-nj))
      dz=z(ijk)-z(ijk-nij)
      bnuslocalw(jk)=(t(ijk-nj)-t(ijk))/dx
!--------------------------------------------------------
!--------[xmax/dt-xdirection]:2.2m / (85 - 15)
!--------------------------------------------------------
      bnuslocalw(jk)=bnuslocalw(jk)*1./1.
!      bnuslocalw(jk)=bnuslocalw(jk)*2.2/70.
      bnusmeanw(j)=bnusmeanw(j)+bnuslocalw(jk)*dz/1.
      ijk=lkk+li(nim)+j

      dx=0.5*(x(ijk)-x(ijk-nj))
      dz=z(ijk)-z(ijk-nij)
      bnuslocale(jk)=(t(ijk+nj)-t(ijk))/dx
!--------------------------------------------------------
!--------[xmax/dt-xdirection]:2.2m x (85 - 15)
!--------------------------------------------------------
!      bnuslocale(jk)=bnuslocale(jk)*2.2/70.
!      bnusmeane(j)=bnusmeane(j)+bnuslocale(jk)*dz/2.2
      bnuslocale(jk)=bnuslocale(jk)*1./1.
      bnusmeane(j)=bnusmeane(j)+bnuslocale(jk)*dz/1.

      bnuselt1(jk)=bnuslocalw(jk)
      bnuselt2(jk)=bnuslocale(jk)

      end do
      end do
!      write(6,*)
!      write(6,*)'======================================'
!      write(6,*)' mean nusselt number: (west)'
!      write(6,*)'======================================'
!      write(6,*) bnusmeanw
!      write(6,*)'======================================'
!      write(6,*)' mean nusselt number: (east)'
!      write(6,*)'======================================'
!      write(6,*) bnusmeane
!      write(6,*)

      bnusoverallw=0.
      bnusoveralle=0.

      do j=2,njm
      dy=y(j)-y(j-1)
      bnusoverallw=bnusoverallw+bnusmeanw(j)*dy/1.
      bnusoveralle=bnusoveralle+bnusmeane(j)*dy/1.
!      bnusoverallw=bnusoverallw+bnusmeanw(j)*dy/1.5
!      bnusoveralle=bnusoveralle+bnusmeane(j)*dy/1.5
      end do
      write(6,*)'=============================================='
      write(6,*)'bnusoverallw= ',bnusoverallw
      write(6,*)'=============================================='
      write(6,*)'=============================================='
      write(6,*)'bnusoveralle= ',bnusoveralle
      write(6,*)'=============================================='
      bns1=bnusoverallw
      bns2=bnusoveralle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     [slucaj grijanja odozdo:]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      do j=2,njm
!      bnusmeanb(j)=0.
!      bnusmeant(j)=0.
!      do i=2,nim
!      ij=(i-1)*nj+j
!      ijk=lk(2)+li(i)+j
!      dz=0.5*(z(ijk)-z(ijk-nij))
!c      if(i.eq.50.and.j.eq.50) then
!c      write(6,*)'dz1= ',dz
!c      write(6,*)'t(ijk-nij)= ',t(ijk-nij)
!c      write(6,*)'t(ijk)= ',t(ijk)
!c      end if
!      dx=x(ijk)-x(ijk-nj)
!c      dy=y(ijk)-y(ijk-1)
!      bnuslocalb(ij)=(t(ijk-nij)-t(ijk))/dz
!      bnusmeanb(j)=bnusmeanb(j)+bnuslocalb(ij)*dx/8.
!      ijk=lk(nkm)+li(i)+j
!      dz=0.5*(z(ijk)-z(ijk-nij))
!c      if(i.eq.50.and.j.eq.50) then
!c      write(6,*)'dz2= ',dz
!c      write(6,*)'t(ijk+nij)= ',t(ijk+nij)
!c      write(6,*)'t(ijk)= ',t(ijk)
!c      end if
!      dx=x(ijk)-x(ijk-nj)
!c      dy=y(ijk)-y(ijk-1)
!      bnuslocalt(ij)=(t(ijk+nij)-t(ijk))/dz
!      bnusmeant(j)=bnusmeant(j)+bnuslocalt(ij)*dx/8.
!      bnuselt1(ij)=bnuslocalb(ij)
!      bnuselt2(ij)=bnuslocalt(ij)
!      end do
!      end do
!      write(6,*)
!      write(6,*)'======================================'
!      write(6,*)' mean nusselt number: (bottom)'
!      write(6,*)'======================================'
!      write(6,*) bnusmeanb
!      write(6,*)'======================================'
!      write(6,*)' mean nusselt number: (top)'
!      write(6,*)'======================================'
!      write(6,*) bnusmeant
!      write(6,*)

!      bnusoverallb=0.
!      bnusoverallt=0.
!      do j=2,njm
!c      ijk=i*nj
!c      dx=x(ijk)-x(ijk-nj)
!      dy=y(j)-y(j-1)
!c      bnusoverallb=bnusoverallb+bnusmeanb(j)*dy/8.
!      bnusoverallt=bnusoverallt+bnusmeant(j)*dy/8.
!      end do
!      write(6,*)'=============================================='
!      write(6,*)'bnusoverallb= ',bnusoverallb
!      write(6,*)'=============================================='
!      write(6,*)'=============================================='
!      write(6,*)'bnusoverallt= ',bnusoverallt
!      write(6,*)'=============================================='

!      bns1=bnusoverallb
!     bns2=bnusoverallt

!      write(6,*)'bns1= ',bns1
!      write(6,*)'bns2= ',bns2

      return
      end
