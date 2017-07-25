!***********************************************************************
!
      subroutine printi(fi,hedfi)
!
!***********************************************************************
!
!     print field values in surfaces of constant  i
!
!***********************************************************************
!
      use types
      use parameters
      use indexes

     implicit none
!
!***********************************************************************
!

     real(prec), dimension(nxyza), intent(in) :: fi
     character(len=10):: hedfi

!     local variables
      integer :: i, j, k, ii, js, je

      write(66,20) hedfi
      do 100 i=1,ni
      ii=li(i)
      write(66,21) i
      js=-11
   10 js=js+12
      je=js+11
      je=min0(je,nj)
      write(66,22) (j,j=js,je)
      write(66,*) '  k'
      do 11 k=nk,1,-1
   11 write(66,23) k,(fi(lk(k)+ii+j),j=js,je)
      if(je.lt.nj) go to 10
  100 continue
      return
   20 format(//,40x,a10,/,40x,10('*'),/)
   21 format(/,2x,26('*-'),7x,' i =',i3,' ',7x,26('-*'))
   22 format(3x,'j = ',i3,11i10)
   23 format(1x,i3,12e10.2)
      end

!***********************************************************************
!
      subroutine printj(fi,hedfi)
!
!***********************************************************************
!
!     print field values in surfaces of constant  j
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      
     implicit none
!
!***********************************************************************
!

     real(prec), dimension(nxyza), intent(in) :: fi
     character(len=10):: hedfi


!     local variables
      integer :: i, j, k, is, ie 

      write(66,20) hedfi
!      do 100 j=1,nj,10
      do 100 j=2,2,1
      write(66,21) j
      is=-11
   10 is=is+12
      ie=is+11
      ie=min0(ie,ni)
      write(66,22) (i,i=is,ie)
      write(66,*) '  k'
      do 11 k=nk,1,-1
   11 write(66,23) k,(fi(lk(k)+li(i)+j),i=is,ie)
      if(ie.lt.ni) go to 10
  100 continue
      return
   20 format(//,40x,a10,/,40x,10('*'),/)
   21 format(/,2x,26('*-'),7x,' j =',i3,' ',7x,26('-*'))
   22 format(3x,'i = ',i3,11i10)
   23 format(1x,i3,12e10.2)
      end


!***********************************************************************
!
      subroutine printk(fi,hedfi)
!
!***********************************************************************
!
!    print field values in surfaces of constant  k
!
!***********************************************************************
!
     use types
     use parameters
     use indexes

     implicit none
!
!***********************************************************************
!

     real(prec), dimension(nxyza), intent(in) :: fi
     character(len=10):: hedfi
!
!     local variables
!
      integer :: i, j, k, lkk, is, ie

      write(66,20) hedfi
      do 100 k=1,nk
      write(66,21) k
      lkk=lk(k)
      is=-11
   10 is=is+12
      ie=is+11
      ie=min0(ie,ni)
      write(66,22) (i,i=is,ie)
      write(66,*) '  j'
      do 11 j=nj,1,-1
   11 write(66,23) j,(fi(lkk+li(i)+j),i=is,ie)
      if(ie.lt.ni) go to 10
  100 continue
      return
   20 format(//,40x,a10,/,40x,10('*'),/)
   21 format(/,2x,26('*-'),7x,' k =',i3,' ',7x,26('-*'))
   22 format(3x,'i = ',i3,11i10)
   23 format(1x,i3,12e10.2)
      end



!***********************************************************************
!
      subroutine print2(nih,njh,hedi,hedj,phi,hedphi)
!
!***********************************************************************
!
      use types
      use parameters

      implicit none
!
!***********************************************************************
!
      character(10):: hedphi
      character(1) :: hedi, hedj
      integer :: nih, njh
      real(prec), dimension(nxyza) :: phi

!     local variables
      integer :: i, j, jj, is, ie

      write(66,20) hedphi
      is=-11
  100 is=is+12
      ie=is+11
      ie=min0(nih,ie)
      write(66,21) hedi,(i,i=is,ie)
      write(66,22) hedj
      do 101 jj=1,njh
      j=njh+1-jj
  101 write(66,23) j,(phi((i-1)*njh+j),i=is,ie)
      if(ie.lt.nih) go to 100
   20 format(2x,25('*-'),7x,a10,7x,25('-*'))
   21 format(3x,a1,' = ',i3,11i10)
   22 format(2x,a1)
   23 format(1x,i3,12e10.2)
      return
      end

