      SUBROUTINE fvm_laplacian(mu,phi)
!  
!******************************************************************************
!
!     Fills matrix with coefficients representing implicit FVM discretization 
!     of negative Laplacian operator: -div(mu*grad(phi)).
!
!     Inteded for structured 3D grids. Non-orthogonal corrections are included.
!
!     Matrix is stored in diagonal format with seven diagonals.
!     Each diagonal is stored in an array,
!     Off-main diagonnal: ae,aw,an,as,at,ab,
!     Main diagonal: ap.
!     RHS vector is SU.
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!
      use types
      use parameters
      use indexes
      use geometry
      use coef
      use variables
      use bc

      implicit none

      real(prec), dimension(nxyza) :: mu,phi
!
!     Local variables
!
      integer :: i, j, k, inp, &
                 ine, ins, inb, inbs,idew, idns, idtb
      integer :: inbc
      real(prec) :: fxe, fxp, &
                    arx, ary, arz, are, &
                    xpn,ypn,zpn,smdpn,sfdpnr
      real(prec) :: nxx,nyy,nzz
      real(prec) :: suu,sup
      real(prec) :: dxet,dzet,dyet,dxzd,dyzd,dzzd

!.....Initialize matrix arrays
      ae=0.0d0; aw = 0.0d0; an = 0.0d0; as = 0.0d0; at=0.d0; ab = 0.0d0; sp = 0.0d0; ap=0.0d0
!
!.....EAST CELL FACES
!
      idew = nj
      idns = 1
      idtb = nij
      do k=2,nkm
      do i=1,nim
      do j=2,njm

      inp=lk(k)+li(i)+j
      ine=inp+idew
      ins=inp-idns
      inb=inp-idtb
      inbs=inb-idns

      fxe=fx(inp)
      fxp=1.-fxe

!/////Correction for triple periodicity://///////////
!@      indx=0
!.....Correction for triple periodicity:
!@      if(i.eq.nie.and.idew.eq.nj)  indx=lk(k)+li(2)+j
!@      if(j.eq.nje.and.idew.eq.1)   indx=lk(k)+li(i)+2
!@      if(k.eq.nke.and.idew.eq.nij) indx=lk(2)+li(i)+j
!@      if(indx/=0) then
!@        INE=indx
!@        FXE=0.5d0
!@        FXP=0.5d0
!@      end if
!/////END: Correction for triple periodicity://///////////

!.....GEOMETRICAL QUANTITIES...............................

!.....Distance between P and E
      xpn=xc(ine)-xc(inp)
      ypn=yc(ine)-yc(inp)
      zpn=zc(ine)-zc(inp)
!/////Correction for triple periodicity://///////////:
!@      if(i.eq.nie.and.idew.eq.nj) then
!@        xpn=x(ine)-xc(ine)+x(inp)-xc(inp); xi=xc(inp)*fxp+(xc(inp)+xpn)*fxe
!@      endif
!@      if(j.eq.nje.and.idew.eq.1)  then
!@        ypn=y(ine)-yc(ine)+y(inp)-yc(inp); yi=yc(inp)*fxp+(yc(inp)+ypn)*fxe
!@      endif
!@      if(k.eq.nke.and.idew.eq.nij) then
!@        zpn=z(ine)-zc(ine)+z(inp)-zc(inp); zi=zc(inp)*fxp+(zc(inp)+zpn)*fxe
!@      endif
!/////END: Correction for triple periodicity://///////////

!.....SECOND
      DXET=.5*(X(INP)-X(INS)+X(INB)-X(INBS))
      DYET=.5*(Y(INP)-Y(INS)+Y(INB)-Y(INBS))
      DZET=.5*(Z(INP)-Z(INS)+Z(INB)-Z(INBS))
!.....THIRD
      DXZD=.5*(X(INP)-X(INB)+X(INS)-X(INBS))
      DYZD=.5*(Y(INP)-Y(INB)+Y(INS)-Y(INBS))
      DZZD=.5*(Z(INP)-Z(INB)+Z(INS)-Z(INBS))
!.....SECOND X THIRD DEFINE  ALWAYS CF AREA
      ARX=DYET*DZZD-DYZD*DZET
      ARY=DXZD*DZET-DXET*DZZD
      ARZ=DXET*DYZD-DYET*DXZD
!.....Precomputed face areas
!      ARX=AR1X(INP)
!      ARY=AR1Y(INP)
!      ARZ=AR1Z(INP)

!.....CELL FACE AREA
      ARE=sqrt(ARX**2+ARY**2+ARZ**2)

!.....Unit vectors of the normal
      nxx=arx/are
      nyy=ary/are
      nzz=arz/are

!.....VOLUME
      sfdpnr=1./(ARX*XPN*nxx+ARY*YPN*nyy+ARZ*ZPN*nzz)
      smdpn=(arx*arx+ary*ary+arz*arz)*sfdpnr
!.....END: GEOMETRICAL QUANTITIES...............................
!
!.....Coefficients of discretized Laplace equation
      ae(inp)=(fxp*mu(inp)+fxe*mu(ine))*smdpn
      aw(ine)=ae(inp) 

      END DO
      END DO
      END DO

!
!.....NORTH CELL FACES
!
      idew = 1
      idns = nij
      idtb = nj
      do k=2,nkm
      do i=2,nim
      do j=1,njm

      inp=lk(k)+li(i)+j
      ine=inp+idew
      ins=inp-idns
      inb=inp-idtb
      inbs=inb-idns

      fxe=fy(inp)
      fxp=1.-fxe

!/////Correction for triple periodicity://///////////
!@      indx=0
!.....Correction for triple periodicity:
!@      if(i.eq.nie.and.idew.eq.nj)  indx=lk(k)+li(2)+j
!@      if(j.eq.nje.and.idew.eq.1)   indx=lk(k)+li(i)+2
!@      if(k.eq.nke.and.idew.eq.nij) indx=lk(2)+li(i)+j
!@      if(indx/=0) then
!@        INE=indx
!@        FXE=0.5d0
!@        FXP=0.5d0
!@      end if
!/////END: Correction for triple periodicity://///////////

!.....GEOMETRICAL QUANTITIES...............................

!.....Distance between P and E
      xpn=xc(ine)-xc(inp)
      ypn=yc(ine)-yc(inp)
      zpn=zc(ine)-zc(inp)
!/////Correction for triple periodicity://///////////:
!@      if(i.eq.nie.and.idew.eq.nj) then
!@        xpn=x(ine)-xc(ine)+x(inp)-xc(inp); xi=xc(inp)*fxp+(xc(inp)+xpn)*fxe
!@      endif
!@      if(j.eq.nje.and.idew.eq.1)  then
!@        ypn=y(ine)-yc(ine)+y(inp)-yc(inp); yi=yc(inp)*fxp+(yc(inp)+ypn)*fxe
!@      endif
!@      if(k.eq.nke.and.idew.eq.nij) then
!@        zpn=z(ine)-zc(ine)+z(inp)-zc(inp); zi=zc(inp)*fxp+(zc(inp)+zpn)*fxe
!@      endif
!/////END: Correction for triple periodicity://///////////

!.....SECOND
      DXET=.5*(X(INP)-X(INS)+X(INB)-X(INBS))
      DYET=.5*(Y(INP)-Y(INS)+Y(INB)-Y(INBS))
      DZET=.5*(Z(INP)-Z(INS)+Z(INB)-Z(INBS))
!.....THIRD
      DXZD=.5*(X(INP)-X(INB)+X(INS)-X(INBS))
      DYZD=.5*(Y(INP)-Y(INB)+Y(INS)-Y(INBS))
      DZZD=.5*(Z(INP)-Z(INB)+Z(INS)-Z(INBS))
!.....SECOND X THIRD DEFINE  ALWAYS CF AREA
      ARX=DYET*DZZD-DYZD*DZET
      ARY=DXZD*DZET-DXET*DZZD
      ARZ=DXET*DYZD-DYET*DXZD
!.....Precomputed face areas
!      ARX=AR2X(INP)
!      ARY=AR2Y(INP)
!      ARZ=AR2Z(INP)
      !if (j==1)print*,'are: ',arx,ary,arz
!.....CELL FACE AREA
      ARE=sqrt(ARX**2+ARY**2+ARZ**2)

!.....Unit vectors of the normal
      nxx=arx/are
      nyy=ary/are
      nzz=arz/are

!.....VOLUME
      sfdpnr=1./(ARX*XPN*nxx+ARY*YPN*nyy+ARZ*ZPN*nzz)
      smdpn=(arx*arx+ary*ary+arz*arz)*sfdpnr

!.....END: GEOMETRICAL QUANTITIES...............................
!
!.....Coefficients of discretized Laplace equation
      an(inp)=(fxp*mu(inp)+fxe*mu(ine))*smdpn
      as(ine)=an(inp) 

      END DO
      END DO
      END DO

!
!.....TOP CELL FACES
!
      idew = nij
      idns = nj
      idtb = 1
      do k=1,nkm
      do i=2,nim
      do j=2,njm

      inp=lk(k)+li(i)+j
      ine=inp+idew
      ins=inp-idns
      inb=inp-idtb
      inbs=inb-idns

      fxe=fz(inp)
      fxp=1.-fxe

!/////Correction for triple periodicity://///////////
!@      indx=0
!.....Correction for triple periodicity:
!@      if(i.eq.nie.and.idew.eq.nj)  indx=lk(k)+li(2)+j
!@      if(j.eq.nje.and.idew.eq.1)   indx=lk(k)+li(i)+2
!@      if(k.eq.nke.and.idew.eq.nij) indx=lk(2)+li(i)+j
!@      if(indx/=0) then
!@        INE=indx
!@        FXE=0.5d0
!@        FXP=0.5d0
!@      end if
!/////END: Correction for triple periodicity://///////////

!.....GEOMETRICAL QUANTITIES...............................

!.....Distance between P and E
      xpn=xc(ine)-xc(inp)
      ypn=yc(ine)-yc(inp)
      zpn=zc(ine)-zc(inp)
      !print*,'dpn: ',xpn,ypn,zpn
!/////Correction for triple periodicity://///////////:
!@      if(i.eq.nie.and.idew.eq.nj) then
!@        xpn=x(ine)-xc(ine)+x(inp)-xc(inp); xi=xc(inp)*fxp+(xc(inp)+xpn)*fxe
!@      endif
!@      if(j.eq.nje.and.idew.eq.1)  then
!@        ypn=y(ine)-yc(ine)+y(inp)-yc(inp); yi=yc(inp)*fxp+(yc(inp)+ypn)*fxe
!@      endif
!@      if(k.eq.nke.and.idew.eq.nij) then
!@        zpn=z(ine)-zc(ine)+z(inp)-zc(inp); zi=zc(inp)*fxp+(zc(inp)+zpn)*fxe
!@      endif
!/////END: Correction for triple periodicity://///////////

!.....SECOND
      DXET=.5*(X(INP)-X(INS)+X(INB)-X(INBS))
      DYET=.5*(Y(INP)-Y(INS)+Y(INB)-Y(INBS))
      DZET=.5*(Z(INP)-Z(INS)+Z(INB)-Z(INBS))
!.....THIRD
      DXZD=.5*(X(INP)-X(INB)+X(INS)-X(INBS))
      DYZD=.5*(Y(INP)-Y(INB)+Y(INS)-Y(INBS))
      DZZD=.5*(Z(INP)-Z(INB)+Z(INS)-Z(INBS))
!.....SECOND X THIRD DEFINE  ALWAYS CF AREA
      ARX=DYET*DZZD-DYZD*DZET
      ARY=DXZD*DZET-DXET*DZZD
      ARZ=DXET*DYZD-DYET*DXZD
!.....Precomputed face areas
!      ARX=AR3X(INP)
!      ARY=AR3Y(INP)
!      ARZ=AR3Z(INP)
!      if(k==1)print*,'are: ',arx,ary,arz
!.....CELL FACE AREA
      ARE=sqrt(ARX**2+ARY**2+ARZ**2)

!.....Unit vectors of the normal
      nxx=arx/are
      nyy=ary/are
      nzz=arz/are
      !print*,'n: ',nxx,nyy,nzz
!.....VOLUME
      sfdpnr=1./(ARX*XPN*nxx+ARY*YPN*nyy+ARZ*ZPN*nzz)
      smdpn=(arx*arx+ary*ary+arz*arz)*sfdpnr

!.....END: GEOMETRICAL QUANTITIES...............................
!
!.....Coefficients of discretized Laplace equation
      at(inp)=(fxp*mu(inp)+fxe*mu(ine))*smdpn
      ab(ine)=at(inp) 

      END DO
      END DO
      END DO

!.....Assemble main diagonal coefficient
      do k=2,nkm
      do i=2,nim
      do j=2,njm
      inp=lk(k)+li(i)+j
      ap(inp)=ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp)
      enddo
      enddo
      enddo

!.....Modify matrix coefficients to reflect presence of Boundary Conditions in PDE problem.

!.....WEST
      do k=2,nkm
      do j=2,njm
      suu=0.0d0
      sup=0.0d0
      inp=lk(k)+li(2)+j
      inbc=(j-1)*nk+k
      lbhlp=lbw(inbc+jks)

!!      if(lbhlp.ne.4) then
!!      ! Zero gradient B.C.
!!      su(inp)=su(inp)+suu
!!      sp(inp)=sp(inp)+sup
!!      aw(inp)=0.0d0
!!      endif

!!      if(lbhlp.eq.4) then
!.....Wall  - Dirichlet (fixed value) condition at boundary
      su(inp) = su(inp)+aw(inp)*phi(inp-nj)
      aw(inp) = 0.0d0
!!      endif

      enddo
      enddo

!.....E A S T
      do k=2,nkm
      do j=2,njm
      suu=0.0d0
      sup=0.0d0
      inp=lk(k)+li(nim)+j
      inbc=(j-1)*nk+k
      lbhlp=lbe(inbc+jks)

!!      if(lbhlp.ne.4) then
!!      ! Zero gradient B.C.
!!      su(inp)=su(inp)+suu
!!      sp(inp)=sp(inp)+sup
!!      ae(inp)=0.0d0
!!      endif

!!      if(lbhlp.eq.4) then
!.....Wall  - Dirichlet (fixed value) condition at boundary
      su(inp) = su(inp)+ae(inp)*phi(inp+nj)
      ae(inp) = 0.0d0
!!      endif

      enddo
      enddo

!.....S O U T H 
      do k=2,nkm
      do i=2,nim
      suu=0.
      sup=0.
      inp=lk(k)+li(i)+2
      inbc=(i-1)*nk+k
      lbhlp=lbs(inbc+iks)

!!      if(lbhlp.ne.4) then
!!      ! Zero gradient B.C.
!!      su(inp)=su(inp)+suu
!!      sp(inp)=sp(inp)+sup
!!      as(inp)=0.0d0
!!      endif

!!      if(lbhlp.eq.4) then
!.....Wall  - Dirichlet (fixed value) condition at boundary
      su(inp) = su(inp)+as(inp)*phi(inp-1)
      as(inp) = 0.0d0
!!      endif

      enddo
      enddo

!.....N O R T H
      do k=2,nkm
      do i=2,nim
      suu=0.0d0
      sup=0.0d0
      inp=lk(k)+li(i)+njm
      inbc=(i-1)*nk+k
      lbhlp=lbn(inbc+iks)

!!      if(lbhlp.ne.4) then
!!      ! Zero gradient B.C.
!!      su(inp)=su(inp)+suu
!!      sp(inp)=sp(inp)+sup
!!      an(inp)=0.0d0
!!      endif

!!      if(lbhlp.eq.4) then
!.....Wall  - Dirichlet (fixed value) condition at boundary
      su(inp) = su(inp)+an(inp)*phi(inp+1)
      an(inp) = 0.0d0
!!      endif

      enddo
      enddo

!.....BOTTOM
      do i=2,nim
      do j=2,njm
      inp=lk(2)+li(i)+j
      inbc=(i-1)*nj+j
      lbhlp=lbb(inbc+ijs)

!!      if(lbhlp.eq.4) then
!.....Wall  - Dirichlet (fixed value) condition at boundary
      su(inp) = su(inp)+ab(inp)*phi(inp-nij)
      ab(inp) = 0.0d0
!!      endif

      enddo
      enddo

!.....TOP
      do i=2,nim
      do j=2,njm
      suu=0.0d0
      sup=0.0d0
      inp=lk(nkm)+li(i)+j
      inbc=(i-1)*nj+j
      lbhlp=lbt(inbc+ijs)

!!      if(lbhlp.ne.4) then
!!      ! Zero gradient B.C.
!!      su(inp)=su(inp)+suu
!!      sp(inp)=sp(inp)+sup
!!      at(inp)=0.0d0
!!      endif

!!      if(lbhlp.eq.4) then
!.....Wall  - Dirichlet (fixed value) condition at boundary
      su(inp) = su(inp)+at(inp)*phi(inp+nij)
      at(inp) = 0.0d0
!!      endif

      enddo
      enddo

      RETURN
      END
