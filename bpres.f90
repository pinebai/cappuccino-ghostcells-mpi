!***********************************************************************
! 
      subroutine bpres(p)
!
!***********************************************************************
! 
!     Set pressures at ghost cells by extrapolation from inner cells.
!
!***********************************************************************
! 
      use types
      use parameters
      use indexes
      use geometry
      use obstacle

      implicit none
!
!***********************************************************************
! 
      real(prec), dimension(nxyza) :: p

!     Locals:
      integer :: i, j, k, ijk
      ! integer :: ij, ik, jk, lkk
      ! integer :: nsa, ioss, ioes, joss, joes, koss, koes


!---------------------------------------------------------
!.....bottom boundary
!---------------------------------------------------------
      do i=2,nim
      do j=2,njm
      ijk=lk(3)+li(i)+j
      p(ijk-nij)=p(ijk)-(p(ijk+nij)-p(ijk))*fz(ijk)
      enddo
      enddo
!---------------------------------------------------------
!.....top boundary
!---------------------------------------------------------
      do i=2,nim
      do j=2,njm
      ijk=lk(nkmm)+li(i)+j
      p(ijk+nij)=p(ijk)+(p(ijk)-p(ijk-nij))*(1.0_dp-fz(ijk-nij))
      enddo
      enddo
!---------------------------------------------------------
!.....west boundary
!---------------------------------------------------------
      do j=2,njm
      do k=2,nkm
      ijk=lk(k)+li(3)+j
      p(ijk-nj)=p(ijk)-(p(ijk+nj)-p(ijk))*fx(ijk)
      enddo
      enddo
!---------------------------------------------------------
!.....east boundary
!---------------------------------------------------------
      do j=2,njm
      do k=2,nkm
      ijk=lk(k)+li(nimm)+j
      p(ijk+nj)=p(ijk)+(p(ijk)-p(ijk-nj))*(1.0_dp-fx(ijk-nj))
      enddo
      enddo
!---------------------------------------------------------
!.....south boundary
!---------------------------------------------------------
      do i=2,nim
      do k=2,nkm
      ijk=li(i)+lk(k)+3
      p(ijk-1)=p(ijk)-(p(ijk+1)-p(ijk))*fy(ijk)
      enddo
      enddo
!---------------------------------------------------------
!.....north boundary
!---------------------------------------------------------
      do i=2,nim
      do k=2,nkm
      ijk=lk(k)+li(i)+njmm
      p(ijk+1)=p(ijk)+(p(ijk)-p(ijk-1))*(1.0_dp-fy(ijk-1))
      enddo
      enddo
! !
! !====================================================================
! !     [OBSTACLE MODIFICATIONS: ]
! !====================================================================
!       IF(IOBST.EQ.0) RETURN
!       DO NSA=1,NOBST
!       IOSS=IOS(NSA)
!       IOES=IOE(NSA)
!       JOSS=JOS(NSA)
!       JOES=JOE(NSA)
!       KOSS=KOS(NSA)
!       KOES=KOE(NSA)

!       DO K=KOSS,KOES
!          LKK=LK(K)
!          DO J=JOSS,JOES
!          JK=(K-KOSS)*(JOES-JOSS+1)+J-JOSS+1
! !---------------------------------
! !..........[west of obstacle wall]
! !---------------------------------
!            IF(IOSS.NE.2.AND.TYPOBW(NSA,JK).EQ.1) THEN
!            IJK=LKK+LI(IOSS-1)+J
!            P(IJK+NJ)=P(IJK)+ &
!                     (P(IJK)-P(IJK-NJ))*(1.0_dp-FX(IJK-NJ))
! !           P(IJK+NJ)=P(IJK)
!            END IF
! !---------------------------------
! !..........[east of obstacle wall]
! !---------------------------------
!            IF(IOES.NE.NIM.AND.TYPOBE(NSA,JK).EQ.1) THEN
!            IJK=LKK+LI(IOES+1)+J
!            P(IJK-NJ)=P(IJK)- &
!                      (P(IJK+NJ)-P(IJK))*FX(IJK)
! !           P(IJK-NJ)=P(IJK)
!            END IF
!          END DO
! !----------------
!          DO I=IOSS,IOES
!          IK=(I-IOSS)*(KOES-KOSS+1)+K-KOSS+1
! !----------------------------------
! !..........[south of obstacle wall]
! !----------------------------------
!            IF(JOSS.NE.2.AND.TYPOBS(NSA,IK).EQ.1) THEN
!            IJK=LKK+LI(I)+JOSS-1
!            P(IJK+1)=P(IJK)+ &
!                     (P(IJK)-P(IJK-1))*(1.0_dp-FY(IJK-1))
! !           P(IJK+1)=P(IJK)
!            END IF
! !----------------------------------
! !..........[north of obstacle wall]
! !----------------------------------
!            IF(JOES.NE.NJM.AND.TYPOBN(NSA,IK).EQ.1) THEN
!            IJK=LKK+LI(I)+JOES+1
!            P(IJK-1)=P(IJK)- &
!                     (P(IJK+1)-P(IJK))*FY(IJK)
! !           P(IJK-1)=P(IJK)
!            END IF
!          END DO
! !----------------
!       END DO
!       DO I=IOSS,IOES
!       DO J=JOSS,JOES
!       IJ=(I-IOSS)*(JOES-JOSS+1)+J-JOSS+1
! !-----------------------------------
! !..........[bottom of obstacle wall]
! !-----------------------------------
!       IF(KOSS.NE.2.AND.TYPOBB(NSA,IJ).EQ.1) THEN
!       IJK=LK(KOSS-1)+LI(I)+J
!       P(IJK+NIJ)=P(IJK)+ &
!                  (P(IJK)-P(IJK-NIJ))*(1.0_dp-FZ(IJK-NIJ))
! !      P(IJK+NIJ)=P(IJK)
!       END IF
! !-----------------------------------
! !..........[top of obstacle wall]
! !-----------------------------------
!       IF(KOES.NE.NKM.AND.TYPOBT(NSA,IJ).EQ.1) THEN
!       IJK=LK(KOES+1)+LI(I)+J
!       P(IJK-NIJ)=P(IJK)- &
!                  (P(IJK+NIJ)-P(IJK))*FZ(IJK)
! !      P(IJK-NIJ)=P(IJK)
!       END IF
!       END DO
!       END DO

! !      DO K=KOSS,KOES
! !      DO I=IOSS,IOES
! !      DO J=JOSS,JOES
! !      IJK=LK(K)+LI(I)+J
! !      P(IJK)=0.
! !      END DO
! !      END DO
! !      END DO
! !-----------------------------------
!       END DO

      return
      end
