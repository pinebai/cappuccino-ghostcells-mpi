!***********************************************************************
!
      subroutine setind(ngr)
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
      integer, intent(in) :: ngr
!
!     local variables
!
      integer :: i, k
      integer :: iipr, jjpr, kkpr
!.....
      ni=nigit(ngr)
      nj=njgit(ngr)
      nk=nkgit(ngr)
!.....
      nim=ni-1
      njm=nj-1
      nkm=nk-1

      nimm = nim-1
      njmm = njm-1
      nkmm = nkm-1

      nij=ni*nj
      nik=ni*nk
      njk=nj*nk
      nijk=nij*nk
      do i=1,ni
        li(i)=(i-1)*nj+ijkgit(ngr)
      end do
      do k=1,nk
        lk(k)=(k-1)*nij
      end do
!
      ijs=isbij(ngr)
      iks=isbik(ngr)
      jks=isbjk(ngr)
      ijks=ijkgit(ngr)
!
      icst=ijks+1
      icen=ijks+nijk
!-----------------------------------------
!.....monitoring point
!-----------------------------------------
      iim=2**(ngr-1)*(imon-1)+1
      jjm=2**(ngr-1)*(jmon-1)+1
      kkm=2**(ngr-1)*(kmon-1)+1
      ijkmon=lk(kkm)+li(iim)+jjm
!-----------------------------------------
!.....pressure ref. point
!-----------------------------------------
      iipr=2**(ngr-1)*(ipr-1)+1
      jjpr=2**(ngr-1)*(jpr-1)+1
      kkpr=2**(ngr-1)*(kpr-1)+1
      ijkpr=lk(kkpr)+li(iipr)+jjpr

      return
      end
