!***********************************************************************
!
subroutine setind
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

!
! Local variables
!
  integer :: i, k

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
    li(i)=(i-1)*nj
  end do

  do k=1,nk
    lk(k)=(k-1)*nij
  end do

  ! Should be removed, not useful
  ijs=0
  iks=0
  jks=0
  ijks=0

  icst=1
  icen=nijk
!-----------------------------------------
!.....monitoring point
!-----------------------------------------
  ijkmon=lk(kmon)+li(imon)+jmon
  
!-----------------------------------------
!.....pressure ref. point
!-----------------------------------------
  ijkpr=lk(kpr)+li(ipr)+jpr

  end subroutine
