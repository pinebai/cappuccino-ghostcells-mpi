!***********************************************************************
!
  subroutine writefiles
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use variables
  use title_mod
  use wall
  use buoy
  use time_mod
  use gradients
  use omega_turb_models
  use statistics

  implicit none
!
!***********************************************************************
!
  integer :: i, j, k, inp, ijk, inbc
  integer :: ingc, iface
  real(dp) :: fxp,fxe,ground_zero,dnw,are,nxf,nyf,nzf,vnp,xtp,ytp,ztp,vtp,ut2,tau,utot2,utan
  real(prec) :: QVortex

!.....WRITE RESTART FILE................................................
  open(unit=3,file=restart_file,form='unformatted')
  rewind 3

  write(3) itime,time
  if(const_mflux) write(3) gradpcmf
  write(3)(f1(inp),inp=icst,icen)
  write(3)(f2(inp),inp=icst,icen)
  write(3)(f3(inp),inp=icst,icen)
  write(3)(u(inp),inp=icst,icen)
  write(3)(v(inp),inp=icst,icen)
  write(3)(w(inp),inp=icst,icen)
  write(3)(p(inp),inp=icst,icen)
  write(3)(te(inp),inp=icst,icen)
  write(3)(ed(inp),inp=icst,icen)
  write(3)(t(inp),inp=icst,icen)
  write(3)(vis(inp),inp=icst,icen)
!  write(3)(visob(inp),inp=icst,icen)
!  write(3)(vart(inp),inp=icst,icen)
!  write(3)(edd(inp),inp=icst,icen)
!  write(3)(ret(inp),inp=icst,icen)
!  write(3)(den(inp),inp=icst,icen)
!  write(3)(utt(inp),inp=icst,icen)
!  write(3)(vtt(inp),inp=icst,icen)
!  write(3)(wtt(inp),inp=icst,icen)
  write(3)(uu(inp),inp=icst,icen)
  write(3)(vv(inp),inp=icst,icen)
  write(3)(ww(inp),inp=icst,icen)
  write(3)(uv(inp),inp=icst,icen)
  write(3)(uw(inp),inp=icst,icen)
  write(3)(vw(inp),inp=icst,icen)
  write(3)(uo(inp),inp=icst,icen)
  write(3)(vo(inp),inp=icst,icen)
  write(3)(wo(inp),inp=icst,icen)
!  write(3)(to(inp),inp=icst,icen)
  write(3)(teo(inp),inp=icst,icen)
  write(3)(edo(inp),inp=icst,icen)
!  write(3)(varto(inp),inp=icst,icen)
!  write(3)(con(inp),inp=icst,icen)
!  write(3)(cono(inp),inp=icst,icen)
!  write(3)(alph(inp),inp=icst,icen)
  rewind 3
  close (3)
!.....END:WRITE RESTART FILE............................................

  write(66,*)'=*=*= Simulation restart files have been written! =*=*='



!----------------------------------------------------------
!     FILE U_END
!     Write the variable values at outlet ghost cells
!----------------------------------------------------------
  ! open(unit=80,file=trim(out_folder_path)//'/u_end')
  ! rewind 80

  ! do k=2,nkm
  ! do j=2,njm
  !     ijk=lk(k)+li(nim)+j
  !    write(80,'(2x,1p7e14.5,2x,a,i3,a,i3)') u(ijk),v(ijk),w(ijk),p(ijk),te(ijk),ed(ijk),t(ijk),'|west => j= ',j,' k= ',k
  ! enddo
  ! enddo

  ! close(80)


!----------------------------------------------------------
!     FILE PLOFILE
!     Write file with a vertical profiles of different vars
!----------------------------------------------------------

  open(unit=81,file=trim(out_folder_path)//'/PLOTFILE')
  rewind 81

  ! Define i,j index position:
  i=nim-10
  j=10
        
      ! inp = lk(3)+li(i)+j 
      ! CALL WALLBC(INP,NIJ,NJ,1)

  do k=2,nkm
  ijk=lk(k)+li(i)+j

  inbc=lk(3)+li(i)+j 
  ingc = inbc - nij

  iface  = ingc

  fxp = fz(ingc)
  fxe = 1.0_dp-fxp

  ! At what height is the wall
  ground_zero = zc(inbc)*fxp+zc(ingc)*fxe
  ! write(*,*) 'ground zero: ', ground_zero

  are = sqrt(ar3x(iface)**2+ar3y(iface)**2+ar3z(iface)**2)
  ! write(*,*) 'Face ares components: ',ar3x(iface),ar3y(iface),ar3z(iface)

  ! Face normals
  nxf = ar3x(iface)/are
  nyf = ar3y(iface)/are
  nzf = ar3z(iface)/are
  ! write(*,*) 'Normals: ', nxf,nyf,nzf

  ! Magnitude of a cell center velocity projected on boundary face normal
  Vnp = U(inbc)*nxf+V(inbc)*nyf+W(inbc)*nzf
  ! write(*,*) 'Normal velocity magnitude Vnp = ',vnp

  ! Tangential velocity components 
  xtp = U(inbc)-Vnp*nxf
  ytp = V(inbc)-Vnp*nyf
  ztp = W(inbc)-Vnp*nzf
  ! write(*,*) 'Tangential velocity components = ',xtp,ytp,ztp
! 
  ! Its magnitude
  Vtp = sqrt(xtp*xtp+ytp*ytp+ztp*ztp)
  ! write(*,*) 'Tangential velocity magnitude Vtp = ',vtp

  ! Tangent direction
  xtp = xtp/vtp
  ytp = ytp/vtp
  ztp = ztp/vtp
  ! write(*,*) 'Tangential direction components: ',xtp,ytp,ztp


  ! projektovanje razlike brzina na pravac tangente u cell centru ijp
  Ut2 = abs( U(inbc)*xtp + V(inbc)*ytp + W(inbc)*ztp) 
  write(*,*) 'ut2 = ', ut2

  ! Distance from wall
  dnw = abs(zc(ijk)-ground_zero)
  ! write(*,*) 'Normal distance of the cell center - dnw = ',dnw

!...
      Utot2=U(inbc)*U(inbc) + V(inbc)*V(inbc) + W(inbc)*W(inbc)  
      !==> Vnp    Un = ( U(inhlp)*ALF + V(inhlp)*BET + W(inhlp)*GAM )              
      !==> Vnp*Vnp   Un2 = Un*Un                                                      
      Utan = sqrt(Utot2 - Vnp*Vnp)    
      write(*,*) 'utan = ', utan

      Utau=sqrt(Utan*VISCOS/(DENSIT*dnw))  
      write(*,*) 'utau = ', utau

      yplus(ijk)=sqrt(UTAU*UTAU)*DNW*DENSIT/VISCOS
       write(*,*) 'ypl = ',yplus(ijk)
!...

  Tau = viscos*Ut2/(zc(inbc)-ground_zero)
  ! write(*,*) 'Tau = ',Tau

  ! u_tau iz mat. def. u_tau=sqrt(Tau_w/rho):                                                                   
  utau=sqrt(Tau/densit)
  write(*,*) 'utau = ',utau

  ! Ypl = rho* \Delta y * u_tau / \mu
  yplus(ijk) = densit*utau*dnw/viscos
  write(*,*) 'ypl = ',yplus(ijk)

  upl(ijk)=u(ijk)/utau  
  ! write(*,*) 'ypl = ',yplus(ijk),'Upl = ',upl(ijk)

  kplus(ijk)=te(ijk)/(utau**2)

  write(81,'(2x,1p9e14.5,2x)') zc(ijk),u(ijk),te(ijk),ed(ijk),yplus(ijk),upl(ijk),kplus(ijk),uw(ijk),utau
  enddo
  close(81)

!----------------------------------------------------------
!     [TECPLOT FILE FORMAT: ]
!----------------------------------------------------------
  open(unit=82,file=trim(out_folder_path)//'/tecplot_file.plt')
  rewind 82

  write(82,*) 'title     = " "'
  ! write(82,*) 'variables = "x"'
  write(82,*) '"y"'
  write(82,*) '"z"'
  write(82,*) '"u"'
  write(82,*) '"v"'
  write(82,*) '"w"'
  write(82,*) '"p"'
  write(82,*) '"tke"'
  write(82,*) '"ed"'
  write(82,*) '"vis"'
  write(82,*) '"uu"'
  write(82,*) '"vv"'
  write(82,*) '"ww"'
  write(82,*) '"QVortex"'
  write(82,*) 'zone t=" "'
  write(82,*) 'i=',nimm, ' ,j=',njmm, ' ,k=',nkmm,', f=point'

  ! Loop over all cells including ghost cells and write their values at
  ! cell centers but, present that to tecplot like they are values at mesh nodes
  ! and create a structured mesh with node centered data.
  do k=2,nkm
  do j=2,njm
  do i=2,nim

  ijk=lk(k)+li(i)+j
  inbc=lk(2)+li(i)+j 

  ! Q criteria for vortex identification
  QVortex = 0.5*(Vorticity(ijk)-Strain(ijk)) 

  write(82,*) xc(ijk),yc(ijk),zc(ijk),u(ijk),v(ijk),w(ijk), &
              p(ijk),te(ijk),ed(ijk),vis(ijk), &
              uu(ijk),vv(ijk),ww(ijk),qvortex
  end do
  end do
  end do

  rewind 82
  close (82)

!----------------------------------------------------------
!     [VTK FILE FORMAT: ]
!----------------------------------------------------------
  call plot_3D_field_vtk (88, trim(out_folder_path)//'/TKE_scalar_field', 'scalar', &
                         'vtk', 'TKE_scalar_field', 'turb.kin.en.',                     &
                         NI, NJ, NK, 1, 1, 1, XC, YC, ZC, TE, 0.0, 0.0)
  
  call plot_3D_field_vtk (84, trim(out_folder_path)//'/Velocity_vector_field', 'vector', &
                         'vtk', 'Velocity_vector_field', 'velocity',                     &
                          NI, NJ, NK, 1, 1, 1, XC, YC, ZC, U, V, W)
!----------------------------------------------------------

      if (ltransient) then
!--------------------------------------------------------------
!    [ writing of the statistics restart file]
!--------------------------------------------------------------
  open(unit=85,file=trim(out_folder_path)//'/statistics1')   ! <- n_sample is here, statistics restart file 1
  open(unit=86,file=trim(out_folder_path)//'/statistics2')   ! <- u_aver, v_aver,... are here, statistics restart file 2
  rewind 85
  rewind 86

  write(85,*) n_sample
  write(86,*) u_aver,v_aver,w_aver, &
              uu_aver,vv_aver,ww_aver, &
              uv_aver,uw_aver,vw_aver,te_aver, &
              te_aver
  close (85)
  close (86)

!--------------------------------------------------------------
!    [ writing of the statistics tecplot file]
!----------------------------------------------------------
  open(unit=87,file=trim(out_folder_path)//'/tecplot_stat.plt') ! <- statistics postprocessing file
  rewind 87

  write(87,*) 'title     = " "'
  write(87,*) 'variables = "x"'
  write(87,*) '"y"'
  write(87,*) '"z"'
  write(87,*) '"u_aver"'
  write(87,*) '"v_aver"'
  write(87,*) '"w_aver"'
  write(87,*) '"uu_aver"'
  write(87,*) '"vv_aver"'
  write(87,*) '"ww_aver"'
  write(87,*) '"uv_aver"'
  write(87,*) '"uw_aver"'
  write(87,*) '"vw_aver"'
  write(87,*) '"te_aver"'
  write(87,*) 'zone t=" "'
  write(87,*) 'i=',nimm, ' ,j=',njmm, ' ,k=',nkmm,', f=point'

  do k=2,nkm
  do j=2,njm
  do i=2,nim

  ijk=lk(k)+li(i)+j

  write(87,*) xc(ijk),yc(ijk),zc(ijk), &
              u_aver(ijk),v_aver(ijk),w_aver(ijk), &   
              uu_aver(ijk),vv_aver(ijk),ww_aver(ijk), &
              uv_aver(ijk),uw_aver(ijk),vw_aver(ijk), &
              te_aver(ijk)
  end do
  end do
  end do

  close(87)
!--------------------------------------------------------------
  endif

  end subroutine
