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

  implicit none

  include 'mpif.h'

  integer :: ik,jk,ij,ijk

 ! NOTE: nproc_char - nproc zapisan levo u vidu stringa!!!1

  open(unit=4,file=adjustl(trim(grid_file))//nproc_char,form='unformatted')
  rewind 4

  read(4)  ni,nj,nk

  !.....stop reading file for a second...
  !.....set parameters                                                   
  call set_parameters
  !.....allocate arrays                                                  
  call allocate_arrays                                             
  !.....end:setting parameters and allocating arrays.....................


  read(4)  nij,nik,njk,nijk,ijs,iks,jks,ijks, &
          (lbw(jk),jk=jks+1,jks+njk),(lbe(jk),jk=jks+1,jks+njk), &
          (lbs(ik),ik=iks+1,iks+nik),(lbn(ik),ik=iks+1,iks+nik), &
          (lbb(ij),ij=ijs+1,ijs+nij),(lbt(ij),ij=ijs+1,ijs+nij), &
          low,loe,los,lon,lob,lot,lout,lthermbc,lconbc
  read(4) (x(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (y(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (z(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (xc(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (yc(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (zc(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (vol(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (fx(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (fy(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (fz(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (ar1x(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (ar1y(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (ar1z(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (ar2x(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (ar2y(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (ar2z(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (ar3x(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (ar3y(ijk),ijk=ijks+1,ijks+nijk)
  read(4) (ar3z(ijk),ijk=ijks+1,ijks+nijk)

  !----------------------------------------------------
  !     [thermal-bc informations: ]
  !----------------------------------------------------
  if(lthermbc.eq.1) then
  read(4)(lbwt(jk),jk=jks+1,jks+njk),(lbet(jk),jk=jks+1,jks+njk), &
         (lbst(ik),ik=iks+1,iks+nik),(lbnt(ik),ik=iks+1,iks+nik), &
         (lbbt(ij),ij=ijs+1,ijs+nij),(lbtt(ij),ij=ijs+1,ijs+nij)
  read(4) twest,teast,tsouth,tnorth,tbottom,ttop
  read(4) qflxwest,qflxeast,qflxsouth,qflxnorth,qflxbottom,qflxtop
  end if

  !----------------------------------------------------
  !     [concentration-bc informations: ]
  !----------------------------------------------------
  if(lconbc.eq.1) then
  read(4)(lbwc(jk),jk=jks+1,jks+njk),(lbec(jk),jk=jks+1,jks+njk), &
         (lbsc(ik),ik=iks+1,iks+nik),(lbnc(ik),ik=iks+1,iks+nik), &
         (lbbc(ij),ij=ijs+1,ijs+nij),(lbtc(ij),ij=ijs+1,ijs+nij)
  read(4) conwest,coneast,consouth,connorth,conbottom,contop
  end if

  ! !----------------------------------------------------
  ! !     [obstacle informations: ]
  ! !----------------------------------------------------
  ! !----------------------------------------------------
  ! !    [velocity: ]
  ! !----------------------------------------------------
  ! if(iobst.eq.1) then
  ! read(4) ios,ioe,jos,joe,kos,koe, &
  !         dnow,dnoe,dnos,dnon,dnob,dnot, &
  !         typobw,typobe,typobs,typobn, &
  !         typobb,typobt
  ! end if
  ! !----------------------------------------------------
  ! !    [thermal: ]
  ! !----------------------------------------------------
  ! if(iobst.eq.1.and.lthermbc.eq.1) then
  ! read(4) typobtw,typobte,typobts,typobtn, &
  !         typobtb,typobtt, &
  !         tobwest,tobeast,tobsouth,tobnorth, &
  !         tobbottom,tobtop
  ! end if
  ! !----------------------------------------------------
  ! !    [concentration: ]
  ! !----------------------------------------------------
  ! if(iobst.eq.1.and.lconbc.eq.1) then
  ! read(4) typobcw,typobce,typobcs,typobcn, &
  !         typobcb,typobct, &
  !         cobwest,cobeast,cobsouth,cobnorth, &
  !         cobbottom,cobtop
  ! end if
  ! !----------------------------------------------------

  close (4)

  ! set control parameters
  call setcon



end subroutine read_grid