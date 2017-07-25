 program rk4projectionCappuccino
!
!*******************************************************************************
! Application:
! Transient solver for Navier-Stokes equations based on Runge-Kutta explicit
! timestepping and projection algorithm for velocity
!
! Usage
!    ./channel <input_file_path> <inlet_file_path> <grid_file_path> -->>
!-->> <monitor_file_path> <restart_file_path> <out_folder_path>      
!
!*******************************************************************************

      use types
      use parameters
      use indexes
      use geometry
      use coef
      use variables
      use title_mod
      use buoy
      use time_mod
      use obstacle
      use bc
      use gradients
      use fieldManipulation ! Volume Weighted Average function

      implicit none

      integer :: istart, iend, istage
      integer :: i, j, k, inp, inpr
      integer :: idew, idns, idtb
      integer :: narg 
      integer, parameter :: ngr = 1  ! For multiblock in future
      integer, parameter :: isol = 5 ! Pressure eq. solver
      integer, parameter :: nstages = 1 ! number of RK4 stages

      real :: start, finish
      real(prec) :: suma
      real(prec) :: magUbarStar, rUAw, gragPplus, flowDirection 
      real(prec) :: ppref
      real(prec) :: dt

      real(prec), dimension(:), allocatable :: Uc,Vc,Wc

      character(len=2) :: trpn
      character(len=5) :: timechar
!                                                                       
!***********************************************************************
!

      ! Check if any command line arguments are found
      narg=command_argument_count()
      if (narg==0.or.narg<5) write(*,'(a)') 'Not enough arguments!'
      call get_command_argument(1,input_file)
      call get_command_argument(2,inlet_file)
      call get_command_argument(3,grid_file)
      call get_command_argument(4,monitor_file)
      call get_command_argument(5,restart_file)
      call get_command_argument(6,out_folder_path)


      ! Open Files
      open(unit=66,file=monitor_file)
      open(unit=7,file=inlet_file)
      open(unit=80,file=trim(out_folder_path)//'/u_end')
      open(unit=81,file=trim(out_folder_path)//'/plotfile')
      rewind 66
      rewind 7
      rewind 80
      rewind 81


      ! Initialize run
      call modinp


      ! Open Files For Data At Monitoring Points 
      open(unit=89,file=trim(out_folder_path)//'/transient_monitoring_points')
      rewind 89
      open(unit=91,file=trim(out_folder_path)//'/transient_monitoring_points_names')
      rewind 91
      do imon=1,mpoints
          read(91, *) trpn
          open(91+imon,file=trim(out_folder_path)//"/transient_monitor_point_"//trpn,access='append')
          !rewind(91+imon)
      end do

      ! Print summary to monitor file
      call print_header

      ! Allocate tentative velocities
      allocate( Uc(nxyza) )  
      allocate( Vc(nxyza) )  
      allocate( Wc(nxyza) )

!===============================================================================
!     [T I M E   L O O P : ]
!===============================================================================
      if(.not.lread) time=0.0d0
      istart = itime+1 
      iend   = itime+numstep

      write(66,*)'Starting time loop'


      time_loop: do itime=1,numstep

      ! Store old variables
      Uo = U
      Vo = V
      Wo = W

      Uc = U
      Vc = V
      Wc = W
      
      To = T
      Teo = Te
      Edo = Ed            

      ! Update time: runTime++;
      time=time+timestep

      dt = timestep

      ! Set inlet boundary conditions at every timestep? ...or only at start
      If(itime.Eq.istart) Call Bcin

      write(66,*)
      write(66,'(a,i0,a,f10.7)') ' Time step no. : ',itime,' Time = ',time
      write(66,*)

      ! Courant number report:
      include 'CourantNo.h'

      call cpu_time(start)


      do istage = 1,nstages

        ! Velocity gradients: 

        if (lstsq) then
          call grad_lsq_qr(u,gradu,2)
          call grad_lsq_qr(v,gradv,2)
          call grad_lsq_qr(w,gradw,2)
        elseif (gauss) then
          call grad_gauss(u,gradu(1,:),gradu(2,:),gradu(3,:))
          call grad_gauss(v,gradv(1,:),gradv(2,:),gradv(3,:))
          call grad_gauss(w,gradw(1,:),gradw(2,:),gradw(3,:))
        endif


        ! Evaluate fluxmass

        ! Cell-face loop:
        ! East Cell-Face
        call fluxmass_plain(nj,1,nij, ar1x,ar1y,ar1z, fx,f1)
        ! North Cell-Face
        call fluxmass_plain(1,nij,nj, ar2x,ar2y,ar2z, fy,f2)
        ! Top   Cell-Face
        call fluxmass_plain(nij,nj,1, ar3x,ar3y,ar3z, fz,f3)

        !// Adjusts the inlet and outlet fluxes to obey continuity, which is necessary for creating a well-posed
        !// problem where a solution for pressure exists.
        if(.not.const_mflux) call outbc


        ! Initialize ths vectors

        su = 0.0_dp
        sv = 0.0_dp 
        sw = 0.0_dp        


        ! Evaluate explicit velocity fluxes

        ! Cell-face loop:
        ! East Cell-Face
        call fluxuvw_explicit(nj,1,nij, ar1x,ar1y,ar1z, fx,f1)
        ! North Cell-Face
        call fluxuvw_explicit(1,nij,nj, ar2x,ar2y,ar2z, fy,f2)
        ! Top   Cell-Face
        call fluxuvw_explicit(nij,nj,1, ar3x,ar3y,ar3z, fz,f3)


        ! Runge-Kutta update for scpecific stage
        call update_velocity_rk4(uc,vc,wc,istage)

        include 'pressureCorrection.f90'

        ! Update values at ghost cells
        call update_values_at_ghostcells

      enddo


      ! Update fluxes using new velocities
      include 'evaluate_fluxmass.f90'

      ! Continuity errors for present time step
      include 'continuityErrors.h'
      

      ! Turbulence equations - Update RANS or SGS model
      if(lles) then
         ! The SGS viscosity of LES model
         call calc_vis_les
      else
         if(lcal(ite))   call calcscm(kgrid,te,gradte  ,ite) 
         if(lcal(ied))   call calcscm(kgrid,ed,graded  ,ied)
         ! Update turbulent eddy viscosity
         if(lcal(ivis))  call modvis
      endif

      ! Temperature eq.
      if(lcal(ien))   call calcsct(kgrid,t,   gradt     ,ien)

      ! Temperature variance eq.
      if(lcal(ivart)) call calcsct(kgrid,vart,gradvart  ,ivart)

      ! Concentration eq.
      if(lcal(icon))  call calcsct(kgrid,con, gradcon   ,icon)





      call cpu_time(finish)
      write(66,'(a,g0.3,a)') 'ExecutionTime = ',finish-start,' s'
      write(66,*)


      if(const_mflux) then
         include 'const_massflux_correction.f90'
      endif

      if(mod(itime,nzapis).eq.0) call writefiles
      call writehistory !<- write monitoring points
      call calc_statistics 

      if(mod(itime,50).eq.0) then 

            Open(Unit=87,File=Trim(Out_Folder_Path)//'/tecplot-vel.plt') 
            Rewind 87
            Write(87,*) 'Title     = " "'
            Write(87,*) 'Variables = "X"'
            Write(87,*) '"Y"'
            Write(87,*) '"Z"'
            Write(87,*) '"U"'
            Write(87,*) '"V"'
            Write(87,*) '"W"'
            Write(87,*) 'Zone T=" "'
            Write(87,*) 'I=',Nimm, ' ,J=',Njmm, ' ,K=',Nkmm,', F=Point'
            Do k=2,nkm
            do j=2,njm
            do i=2,nim
            Inp=Lk(K)+Li(I)+J
            Write(87,*) Xc(Inp),Yc(Inp),Zc(Inp),U(Inp),V(Inp),W(Inp)
            Enddo
            Enddo
            Enddo 
            Close(87)

           write(timechar,'(i5)') itime
           call execute_command_line("tec360 -b -p makro-film.mcr")
           call execute_command_line("mv untitled100.png "//timechar//".png")
      endif 

      call flush(66)

!===============================================================================
      end do time_loop

      !if(loute) call report_wall_stresses
      if(lwrite) call writefiles


      call deallocate_arrays

      if (allocated(Uc)) deallocate(Uc)
      if (allocated(Vc)) deallocate(Vc)
      if (allocated(Wc)) deallocate(Wc)

      end program rk4projectionCappuccino