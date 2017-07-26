!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
program cappuccino_ghostcells_mpi
!
!***********************************************************************
!
!  Three dimensional finite volume solver for Navier-Stokes equations
!  for structured grids, with collocated variable arrangement.
!  Boundary conditions are treated using ghost cells.
!  
! ./channel <input_file> <inlet_file> <grid_file> <monitor_file> <restart_file> <out_folder_path>      
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  use types
  use parameters
  use indexes
  use geometry
  use coef
  use variables
  use title_mod
  use buoy
  use time_mod
  use gradients
  use fieldManipulation
  use utils

  implicit none

  include 'mpif.h'

  integer :: ijk, i, j, k, inp, itimes, itimee
  integer :: narg
  real(prec) :: source
  real(prec) :: magUbarStar, rUAw, gragPplus, flowDirection
  real(prec) :: suma,dt
  real :: start, finish
  character(len=2) :: trpn
  character(len=6) :: timechar
!                                                                       
!***********************************************************************
!

  !----------------------------------------------------------------------
  ! MPI start up

    call MPI_INIT(ierr)                   ! <- ierr=0 if successful

    call MPI_COMM_SIZE(MPI_COMM_WORLD, &  ! <- ???
                                 nproc,&  ! <- Get the total number of PEs
                                 ierr  )  ! <- ierr=0 if successful

    call MPI_COMM_RANK(MPI_COMM_WORLD, &  ! <- ???
                                 myid, &  ! <- Get the ID of the PE
                                 ierr  )  ! <- ierr=0 if successful
    
    this = myid + 1

    if(nproc == 1) then
      nproc = 0
      this = 0
    endif

    write(*,'(2(a,i2))') ' np = ', nproc, ' myid = ', myid
  !----------------------------------------------------------------------


  ! Check command line arguments
  narg=command_argument_count()
  if (narg==0.or.narg<5) then
    write(*,*) 'Not enough arguments - exiting!'
    stop
  endif
  call get_command_argument(1,input_file)
  call get_command_argument(2,inlet_file)
  call get_command_argument(3,grid_file)
  call get_command_argument(4,monitor_file)
  call get_command_argument(5,restart_file)
  call get_command_argument(6,out_folder_path)


  if (myid .eq. 0) then

    ! Open files
    call openfiles

    ! Print cappuccino logo to log file.
    call show_logo

  endif


  ! Read input file
  call read_input

  ! Read grid file
  call read_grid


      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      ierr = 1
      CALL MPI_FINALIZE(ierr)
      STOP

  ! Initialisation
  call init

  ! Open files for data at monitoring points 
  if(ltransient) then
  open(unit=89,file=trim(out_folder_path)//'/transient_monitoring_points')
  open(unit=91,file=trim(out_folder_path)//'/transient_monitoring_points_names')
  rewind 89
  rewind 91
    do imon=1,mpoints
      read(91, *) trpn
      open(91+imon,file=trim(out_folder_path)//"/transient_monitor_point_"//trpn,access='append')
      ! rewind(91+imon)
    end do
  end if

  ! Initial output
  call print_header

!
!--------------------------------------------------------------------------------------------------
!     [t i m e   l o o p : ]
!--------------------------------------------------------------------------------------------------
!
  if(.not.lread) time=0.0d0
  itimes=itime+1 
  itimee=itime+numstep

  time_loop: do itime=itimes,itimee

    ! Shift variables in time:
    do ijk=1,nijk
       if(bdf) then
       uoo(ijk)=uo(ijk)
       voo(ijk)=vo(ijk)
       woo(ijk)=wo(ijk)
       too(ijk)=to(ijk)
       teoo(ijk)=teo(ijk)
       edoo(ijk)=edo(ijk)
       vartoo(ijk)=varto(ijk)
       conoo(ijk)=cono(ijk)
       endif
       uo(ijk)=u(ijk)
       vo(ijk)=v(ijk)
       wo(ijk)=w(ijk)
       to(ijk)=t(ijk)
       teo(ijk)=te(ijk)
       edo(ijk)=ed(ijk)            
       varto(ijk)=vart(ijk)
       cono(ijk)=con(ijk)
    end do

    ! Set inlet boundary conditions
    if(itime.eq.itimes) call bcin

    ! Courant number report:
    include 'CourantNo.h'

! 
!--------------------------------------------------------------------------------------------------
!.....ITERATION LOOP
!--------------------------------------------------------------------------------------------------
!
    iteration_loop: do iter=1,maxit

      call cpu_time(start)

      ! Update values at ghost cells
      call update_values_at_ghostcells

      ! Calculate velocities. Momentum predictor for PISO.
      if(lcal(iu))    call calcuvw
      ! CALL CORVEL
      ! Update OUTLET BC.
      if(lcal(iu).and..not.const_mflux)    call outbc  
      ! Pressure-velocity coupling. Two options: SIMPLE and PISO
      if(lcal(ip).and.simple)   CALL CALCP
      if(lcal(ip).and.piso)     CALL PISO_multiple_correction
      if(lcal(ip).and.pimple)   CALL PIMPLE_multiple_correction
      ! Turbulence equations
      if(lcal(ite))   call calcscm(te,gradte  ,ite) 
      if(lcal(ied))   call calcscm(ed,graded  ,ied)
      ! Temperature , temperature variance, and concentration eqs.
      if(lcal(ien))   call calcsct(t,   gradt     ,ien)
      if(lcal(ivart)) call calcsct(vart,gradvart  ,ivart)
      if(lcal(icon))  call calcsct(con, gradcon   ,icon)
      ! Update eddy viscosity
      if( lcal(ivis) .and. .not. lles)  call modvis
      ! The SGS viscosity of les model
      if(lcal(ivis) .and. lles)         call calc_vis_les

      call cpu_time(finish)
      write(66,'(a,g0.3,a)') 'ExecutionTime = ',finish-start,' s'
      write(66,*)

      !---------------------------------------------------------------------------------------------
      ! Residual normalization, convergence check  
      !---------------------------------------------------------------------------------------------
      do i=1,nphi
      resor(i)=resor(i)*snorin(i)
      end do 
      ! resor(ied)=resor(ied)*small
      ! resor(icon)=resor(icon)*small

      ! Write to monitor file
      ! include 'simpleMonitorResiduals.h'

      source=max(resor(iu),resor(iv),resor(iw),resor(ip)) 
      if(source.gt.slarge) then
          write(66,"(//,10x,a)") "*** Program terminated -  iterations diverge ***" 
          stop ! zavrsi program
      endif

      !---------------------------------------------------------------------------------------------
      ! False time steping: jump to next line after the end of time loop
      !---------------------------------------------------------------------------------------------
      if(.not.ltransient) then
          if(source.lt.sormax) exit time_loop
      end if
!      
!--------------------------------------------------------------------------------------------------
!     nestacionar: 
!     imas dve mogucnosti ovde:
! 1) nije konvergirao a potrosio je predvidje ni broj iteracija za ovaj vremenski korak
! 2) konvergirao je pre nego sto je zavrsio sa svim predvidjenim iteracijama 
!--------------------------------------------------------------------------------------------------

      if(ltransient) then 

          ! konverigao u okviru timestep-a ili potrosio sve iteracije
          if(source.lt.sormax.or.iter.eq.maxit) then 

            if(const_mflux) then
             include 'const_massflux_correction.f90'
            endif

            if(mod(itime,nzapis).eq.0) call writefiles
            call writehistory !<- write monitoring points
            call calc_statistics 

            if(mod(itime,50).eq.0) then 

               Open(Unit=87,File=Trim(Out_Folder_Path)//'/tecplot-vel.plt') 
               Rewind 87
               Write(87,*) 'Title     = " Velocity field - snapshot"'
               Write(87,*) 'Variables = "X"'
               Write(87,*) '"Y"'
               Write(87,*) '"Z"'
               Write(87,*) '"U"'
               Write(87,*) '"V"'
               Write(87,*) '"W"'
               Write(87,*) 'Zone T=" "'
               Write(87,*) 'I=',Nimm, ' ,J=',Njmm, ' ,K=',Nkmm,', F=Point'
               do k=2,nkm
               do j=2,njm
               do i=2,nim
               Inp=Lk(K)+Li(I)+J
               Write(87,*) Xc(Inp),Yc(Inp),Zc(Inp),U(Inp),V(Inp),W(Inp)
               Enddo
               Enddo
               Enddo 
               Close(87)

               write(timechar,'(i6)') itime
               call execute_command_line("tec360 -b -p makro-film.mcr")
               call execute_command_line("mv untitled100.png "//timechar//".png")
            endif
           
            cycle time_loop 
          endif

      end if 

    end do iteration_loop

    ! Zapis poslije svakih [nzapis] iteracija:
    if(.not.ltransient) then
      if(mod(itime,nzapis).eq.0.and.itime.ne.numstep) call writefiles
    endif

    if(ltransient) call flush(66)
 
  end do time_loop

  ! False time stepping comes here after time loop with exit command

  call corvel
  if(loute) call output               !<- 
  if(loute) call report_wall_stresses !<- vrati kad ti bude trebalo y+, Tau na zidovima


!--------------------------------------------------------------------------------------------------
!.....Calculate Nusselt number distributions:
!--------------------------------------------------------------------------------------------------
  ! if(lthermbc.eq.1) then
  ! call nusnumb
  ! write(90,*) time,bnuselt1,bnuselt2
  ! rewind 90
  ! end if

  ! Write field values for the next run 
  if(lwrite) call writefiles

 ! call deallocate_arrays

  !----------------------------------------------------------------------
  ! MPI final call
    call MPI_Finalize(ierr) ! ierr=0 if successful

end program

