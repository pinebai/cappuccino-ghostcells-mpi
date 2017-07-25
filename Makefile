#
# Makefile for CHANNEL program
#

# Compiler flags:
LFLAGS = -llapack 
F90FLAGS = -Wall -O2 

# Compiler:
F90 = mpif90


LINEAR_SOLVER_FILES=\
      sipsol.f90 \
      cgstab.f90 \
      cgstab_sip.f90 \
      pcg-jacobi.f90 \
      iccg.f90 \
      pcg-sip.f90 \

F90FILES=\
      modules_allocatable.f90 \
      fieldManipulation.f90 \
      set_parameters.f90 \
      allocate.f90 \
      asm_terms.f90 \
      bcin.f90 \
      bpres.f90 \
      calcheatflux.f90 \
      calcp-multiple_correction_SIMPLE.f90 \
      calcscm.f90 \
      calcsct.f90 \
      calcstress.f90 \
      calcuvw.f90 \
      calculate_stuff_for_sst.f90 \
      calculate_stuff_for_earsm.f90 \
      calc_statistics.f90 \
      calc_vis_les.f90 \
      correctBoundaryConditions.f90 \
      fluxmass2.f90 \
      fluxscm-overrelaxed-correction.f90 \
      fluxsct.f90 \
      fluxuvw-overrelaxed-correction.f90 \
      fvm_laplacian.f90 \
      corvel.f90 \
      find_strain_rate.f90 \
      find_intersection_point.f90 \
      fluxmc.f90 \
      intfac.f90 \
      matrix.f90 \
      modinp.f90 \
      modvis.f90 \
      print_header.f90 \
      report_wall_stresses.f90 \
      nusnumb.f90 \
      openfiles.f90 \
      outbc.f90 \
      output.f90 \
      PISO_multiple_correction.f90 \
      PISO_assemble_pressure_eq.f90 \
      PISO_getHbyA.f90 \
      PIMPLE_multiple_correction.f90 \
      print.f90 \
      plot_vtk.f90 \
      readfiles.f90 \
      read_input.f90 \
      read_grid.f90 \
      random_seed.f90 \
      setind.f90 \
      setcon.f90 \
      show_logo.f90 \
      update_values_at_ghostcells_nonperiodic.f90 \
      writefiles.f90 \
      writehistory.f90 \
      GRAD_LSQ.f90 \
      GRAD_LSQ_QR.f90 \
      GRAD_GAUSS.f90 \
      wallbc-blended.f90 \
      main.f90 


#
# How to create object files:
#
LINEAR_SOLVERS = ${LINEAR_SOLVER_FILES:.f90=.o}
F90OBJS = ${F90FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

all: cappuccino_ghostcells_mpi

cappuccino_ghostcells_mpi: ${F90OBJS} ${LINEAR_SOLVERS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} ${LINEAR_SOLVERS} ${LFLAGS} ${INCS} -o cappuccino_ghostcells_mpi 

clean:
	@rm  *.o *.mod

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 

.f90.o:
	${F90} ${F90FLAGS} -c ${INCS}  ${@:.o=.f90}
