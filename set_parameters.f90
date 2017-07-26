!***********************************************************************
!
subroutine set_parameters
!
!***********************************************************************
!
!  Needed to allocate arrays, most of this is because of the old
!  multigrid code.   
!
!***********************************************************************
!
  use parameters
  use indexes

  implicit none

  integer, parameter :: ngit = 1
!
!***********************************************************************
!

  nicv = ni-2
  njcv = nj-2 
  nkcv = nk-2 

  nxyo=nxo*nyo
  nxzo=nxo*nzo
  nyzo=nyo*nzo

  nx=nicv*2**(ngit-1)+2
  ny=njcv*2**(ngit-1)+2
  nz=nkcv*2**(ngit-1)+2

  nxy=nx*ny
  nxz=nx*nz
  nyz=ny*nz
  nxyz=nxy*nz

  nxa=nicv*(2**ngit-1)+2*ngit
  nya=njcv*(2**ngit-1)+2*ngit
  nza=nkcv*(2**ngit-1)+2*ngit

  nxya=nicv*njcv*(4**ngit-1)/3+2*(nxa+nya-2*ngit)
  nxza=nicv*nkcv*(4**ngit-1)/3+2*(nxa+nza-2*ngit)
  nyza=njcv*nkcv*(4**ngit-1)/3+2*(nya+nza-2*ngit)

  nxyza=nicv*njcv*nkcv*(8**ngit-1)/7+ &
  2*(4**ngit-1)/3*(nicv*njcv+nicv*nkcv+njcv*nkcv)+ &
  4*(nxa+nya+nza-4*ngit)

  nxyzm=nxyza-nxyz+1

  return
  end
