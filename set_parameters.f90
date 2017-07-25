!***********************************************************************
!
      SUBROUTINE SET_PARAMETERS
!
!***********************************************************************
!
!     NEEDED TO ALLOCATE ARRAYS
!
!***********************************************************************
!
      USE PARAMETERS
      USE INDEXES

      IMPLICIT NONE
!
!***********************************************************************
!
      NICV = NI-2 
      NJCV = NJ-2 
      NKCV = NK-2 

      NXYO=NXO*NYO
      NXZO=NXO*NZO
      NYZO=NYO*NZO

      NX=NICV*2**(NGIT-1)+2
      NY=NJCV*2**(NGIT-1)+2
      NZ=NKCV*2**(NGIT-1)+2

      NXY=NX*NY
      NXZ=NX*NZ
      NYZ=NY*NZ
      NXYZ=NXY*NZ

      NXA=NICV*(2**NGIT-1)+2*NGIT
      NYA=NJCV*(2**NGIT-1)+2*NGIT
      NZA=NKCV*(2**NGIT-1)+2*NGIT

      NXYA=NICV*NJCV*(4**NGIT-1)/3+2*(NXA+NYA-2*NGIT)
      NXZA=NICV*NKCV*(4**NGIT-1)/3+2*(NXA+NZA-2*NGIT)
      NYZA=NJCV*NKCV*(4**NGIT-1)/3+2*(NYA+NZA-2*NGIT)

      NXYZA=NICV*NJCV*NKCV*(8**NGIT-1)/7+ &
      2*(4**NGIT-1)/3*(NICV*NJCV+NICV*NKCV+NJCV*NKCV)+ &
      4*(NXA+NYA+NZA-4*NGIT)

      NXYZM=NXYZA-NXYZ+1

      RETURN
      END
