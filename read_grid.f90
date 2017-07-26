subroutine read_grid
!********************************************************************************
!  Open & read grid file
!********************************************************************************
  use types
  use parameters
  use indexes
  use geometry
  use bc
  use boundc
  use title_mod
  use buoy
  use time_mod
  use printing
  use utils

  implicit none

  ! include 'mpif.h'

  integer :: ik,jk,ij,ijk
  integer :: grid_unit
  character( len = 5) :: nproc_char 



  ! NOTE: nproc_char <- this (=myid + 1) zapisan levo u vidu stringa.
  call i4_to_s_left ( this, nproc_char )
  print*, 'Opening grid file: ',adjustl(trim(grid_file))//'-'//trim(nproc_char)

  call get_unit ( grid_unit )

  open(unit=grid_unit,file=adjustl(trim(grid_file))//'-'//trim(nproc_char),form='unformatted')
  rewind grid_unit

  ! read( grid_unit ) &
  !       ni,nj,nk,nijk, &
  !       ninl,nout,nsym,npru,nwal,noc, &
  !       nwali,nwala,nwalf

! !-----------------------------------------------------------------------
! ! stop reading file for a second and do some useful work
! !-----------------------------------------------------------------------
! call set_parameters                                              
! call allocate_arrays                                             
! !-----------------------------------------------------------------------


!   ! read( grid_unit ) &
!   !       (li(i),i=1,ni),(lk(k),k=1,nk)

!   if(ninl.gt.0) then     
!   read( grid_unit ) &
!         (iji(i)  ,i=1,ninl),(ijpi(i) ,i=1,ninl), &
!         (xni(i),i=1,ninl),(yni(i),i=1,ninl),(zni(i),i=1,ninl), &
!         (xfi(i),i=1,ninl),(yfi(i),i=1,ninl),(zfi(i),i=1,ninl)
!   endif

!   if(nout.gt.0) then 
!   read( grid_unit ) &
!         (ijo(i)  ,i=1,nout),(ijpo(i) ,i=1,nout), &
!         (xno(i),i=1,nout),(yno(i),i=1,nout),(zno(i),i=1,nout), &
!         (xfo(i),i=1,nout),(yfo(i),i=1,nout),(zfo(i),i=1,nout)
!   endif

!   if(nwal.gt.0) then 
!   read( grid_unit ) &
!         (ijw(i)  ,i=1,nwal),(ijpw(i) ,i=1,nwal), &
!         (srdw(i),i=1,nwal),(dnw(i),i=1,nwal), &
!         (xnw(i),i=1,nwal),(ynw(i),i=1,nwal),(znw(i),i=1,nwal), &
!         (xfw(i),i=1,nwal),(yfw(i),i=1,nwal),(zfw(i),i=1,nwal)
!   endif

!   if(nsym.gt.0) then 
!   read( grid_unit ) &
!         (ijs(i)  ,i=1,nsym),(ijps(i) ,i=1,nsym), &
!         (srds(i),i=1,nsym),(dns(i),i=1,nsym), &
!         (xns(i),i=1,nsym),(yns(i),i=1,nsym),(zns(i),i=1,nsym), &
!         (xfs(i),i=1,nsym),(yfs(i),i=1,nsym),(zfs(i),i=1,nsym)
!   endif

!   if(npru.gt.0) then 
!   read( grid_unit ) &                    
!         (iju(i)  ,i=1,npru),(ijpu(i) ,i=1,npru), &
!         (xnpr(i),i=1,npru),(ynpr(i),i=1,npru),(znpr(i),i=1,npru), &
!         (xfpr(i),i=1,npru),(yfpr(i),i=1,npru),(zfpr(i),i=1,npru)
!   endif

!   if(noc.gt.0) then 
!   read( grid_unit ) &
!         (ijl(i)  ,i=1,noc),(ijr(i)  ,i=1,noc), &
!         (srdoc(i), i=1,noc),(foc(i), i=1,noc), &
!         (xnoc(i),i=1,noc),(ynoc(i),i=1,noc),(znoc(i),i=1,noc), &
!         (xfoc(i),i=1,noc),(yfoc(i),i=1,noc),(zfoc(i),i=1,noc)
!   endif

!   read( grid_unit ) &
! !
! !           cell data-cell centers and volumes
! !
!         (xc(i) ,i=1,nijk),(yc(i) ,i=1,nijk),(zc(i) ,i=1,nijk), &
!         (vol(i),i=1,nijk), &
! !
! !             interpolation factors in three directions
! !
!         (fx(i) ,i=1,nijk),(fy(i) ,i=1,nijk),(fz(i) ,i=1,nijk), &
! !
! !             face normal vector - its components are face area projections
! !
!         (ar1x(i),i=1,nijk),(ar1y(i),i=1,nijk),(ar1z(i),i=1,nijk), &
!         (ar2x(i),i=1,nijk),(ar2y(i),i=1,nijk),(ar2z(i),i=1,nijk), &
!         (ar3x(i),i=1,nijk),(ar3y(i),i=1,nijk),(ar3z(i),i=1,nijk), &
! !
! !             face centers for inner faces
! !
!         (xf1(i),i=1,nijk),(yf1(i),i=1,nijk),(zf1(i),i=1,nijk), &
!         (xf2(i),i=1,nijk),(yf2(i),i=1,nijk),(zf2(i),i=1,nijk), &
!         (xf3(i),i=1,nijk),(yf3(i),i=1,nijk),(zf3(i),i=1,nijk)

! !.....Close mesh file
!   close (4)



  read( grid_unit )  ni,nj,nk

!-----------------------------------------------------------------------
! stop reading file for a second and do some useful work
!-----------------------------------------------------------------------
call set_parameters                                              
call allocate_arrays                                             
!-----------------------------------------------------------------------


  read( grid_unit )  nij,nik,njk,nijk,ijs,iks,jks,ijks, &
    (lbw(jk),jk=1,njk),(lbe(jk),jk=1,njk), &
    (lbs(ik),ik=1,nik),(lbn(ik),ik=1,nik), &
    (lbb(ij),ij=1,nij),(lbt(ij),ij=1,nij), &
    low,loe,los,lon,lob,lot,lout,lthermbc,lconbc

  read( grid_unit ) (x(ijk),ijk=1,nijk)
  read( grid_unit ) (y(ijk),ijk=1,nijk)
  read( grid_unit ) (z(ijk),ijk=1,nijk)
  read( grid_unit ) (xc(ijk),ijk=1,nijk)
  read( grid_unit ) (yc(ijk),ijk=1,nijk)
  read( grid_unit ) (zc(ijk),ijk=1,nijk)
  read( grid_unit ) (vol(ijk),ijk=1,nijk)
  read( grid_unit ) (fx(ijk),ijk=1,nijk)
  read( grid_unit ) (fy(ijk),ijk=1,nijk)
  read( grid_unit ) (fz(ijk),ijk=1,nijk)
  read( grid_unit ) (ar1x(ijk),ijk=1,nijk)
  read( grid_unit ) (ar1y(ijk),ijk=1,nijk)
  read( grid_unit ) (ar1z(ijk),ijk=1,nijk)
  read( grid_unit ) (ar2x(ijk),ijk=1,nijk)
  read( grid_unit ) (ar2y(ijk),ijk=1,nijk)
  read( grid_unit ) (ar2z(ijk),ijk=1,nijk)
  read( grid_unit ) (ar3x(ijk),ijk=1,nijk)
  read( grid_unit ) (ar3y(ijk),ijk=1,nijk)
  read( grid_unit ) (ar3z(ijk),ijk=1,nijk)

  ! !----------------------------------------------------
  ! !     [thermal-bc informations: ]
  ! !----------------------------------------------------
  ! if(lthermbc.eq.1) then
  ! read( grid_unit )(lbwt(jk),jk=1,njk),(lbet(jk),jk=1,njk), &
  !        (lbst(ik),ik=1,nik),(lbnt(ik),ik=1,nik), &
  !        (lbbt(ij),ij=1,nij),(lbtt(ij),ij=1,nij)
  ! read( grid_unit ) twest,teast,tsouth,tnorth,tbottom,ttop
  ! read( grid_unit ) qflxwest,qflxeast,qflxsouth,qflxnorth,qflxbottom,qflxtop
  ! end if

  ! !----------------------------------------------------
  ! !     [concentration-bc informations: ]
  ! !----------------------------------------------------
  ! if(lconbc.eq.1) then
  ! read( grid_unit )(lbwc(jk),jk=1,njk),(lbec(jk),jk=1,njk), &
  !        (lbsc(ik),ik=1,nik),(lbnc(ik),ik=1,nik), &
  !        (lbbc(ij),ij=1,nij),(lbtc(ij),ij=1,nij)
  ! read( grid_unit ) conwest,coneast,consouth,connorth,conbottom,contop
  ! end if

  ! !----------------------------------------------------
  ! !     [obstacle informations: ]
  ! !----------------------------------------------------
  ! !----------------------------------------------------
  ! !    [velocity: ]
  ! !----------------------------------------------------
  ! if(iobst.eq.1) then
  ! read( grid_unit ) ios,ioe,jos,joe,kos,koe, &
  !         dnow,dnoe,dnos,dnon,dnob,dnot, &
  !         typobw,typobe,typobs,typobn, &
  !         typobb,typobt
  ! end if
  ! !----------------------------------------------------
  ! !    [thermal: ]
  ! !----------------------------------------------------
  ! if(iobst.eq.1.and.lthermbc.eq.1) then
  ! read( grid_unit ) typobtw,typobte,typobts,typobtn, &
  !         typobtb,typobtt, &
  !         tobwest,tobeast,tobsouth,tobnorth, &
  !         tobbottom,tobtop
  ! end if
  ! !----------------------------------------------------
  ! !    [concentration: ]
  ! !----------------------------------------------------
  ! if(iobst.eq.1.and.lconbc.eq.1) then
  ! read( grid_unit ) typobcw,typobce,typobcs,typobcn, &
  !         typobcb,typobct, &
  !         cobwest,cobeast,cobsouth,cobnorth, &
  !         cobbottom,cobtop
  ! end if
  ! !----------------------------------------------------

  close ( grid_unit )

  ! Set grid loop parameters
  call setind

print*,'Here',myid

end subroutine read_grid