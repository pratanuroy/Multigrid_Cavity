# Set the compiler
#FC = ifort
FC = mpif90
#FC = mpiifort
#LIBDIR = /usr/local/intel/mkl91023/lib/em64t
#LIBS = -L${LIBDIR} -lmkl_em64t -lmkl_lapack -lguide -lpthread
#LIBS = -mkl -lmkl_lapack95_lp64

# Flag for optimization
#OptFLAG = -g -xHost
OptFLAG = -O3

# Array bound check and other debugging flags
#ArrayBoundFlag = -check bounds  -check all -traceback -override-limits
#ArrayBoundFlag = -O2 -check all -traceback -warn all
#ArrayBoundFlag = -O2  -stand f90 -assume realloc_lhs  -check all  -traceback  -warn all  -fstack-protector  -assume protect_parens  -implicitnone

# List of programs
NS_PROG = main.f90 global_variables.f90 read_data.f90 prelim_data.f90 multigrid_levels.f90 VTR_mod.f90 kind_parameters.f90\
setup_num_levels.f90 initial_conditions.f90 properties.f90 boundary_conditions.f90 setup_mg_variables.f90 setup_mg_geometry.f90 \
coordinates.f90 print_results.f90 reset_pressure.f90 apfun.f90 line_gs_solver.f90 tdma.f90 velocity_coefficients.f90\
momentum_interpolation.f90 pcoefficients.f90 correction.f90 convergence.f90 convergencep.f90 monitor_residuals.f90 parameters.f90 \
unsteady_functions.f90 source_terms.f90 steady_solver_sg.f90 broadcast_data.f90 check_convergence.f90 Heaviside.f90 write_data.f90
#
MG_PROG = multigrid_vcycle.f90 mg_residual.f90 fine_to_coarse.f90 coarse_to_fine.f90 fine_to_coarse_res.f90 mg_source.f90 mg_error.f90 mg_correction.f90 face_value.f90 \
fine_to_coarse_uface.f90 fine_to_coarse_vface.f90 momentum_interpolation_mg.f90 coarse_to_fine_uface.f90 coarse_to_fine_vface.f90 fine_to_coarse_uface_sum.f90 \
fine_to_coarse_vface_sum.f90 mg_correction_face_vel.f90 correction_coarse.f90 calculate_face_velocity.f90 \
 calculate_face_velocity_interp.f90 steady_solver_mg.f90 deallocate_mg_variables.f90

PARALLEL_PROG = parallel_variables.f90 map_task.f90 exchange_data.f90 
SOURCES = $(NS_PROG) $(MG_PROG) $(PARALLEL_PROG)
OBJECTS = $(SOURCES:.f90=.o)
EXECUTABLE = multigrid.exe


#all:$(EXECUBTALE)
$(EXECUTABLE): $(OBJECTS)
	$(FC) $(OptFLAG) $(OBJECTS) -o $(EXECUTABLE)
%.mod:%.o %.f90
	$(FC) $(OptFLAG) $(ArrayBoundFlag) -c -o $@ $<

%.o:%.f90
	$(FC) -c  $(OptFLAG) $(ArrayBoundFlag) $<
main.o: kind_parameters.mod global_variables.mod multigrid_levels.mod VTR_mod.mod  parallel_variables.mod parameters.mod main.f90

	$(FC) -c $(OptFLAG) $(ArrayBoundFlag) main.f90

%.o:%.mod %.f90
	$(FC) -c $(OptFLAG) $(ArrayBoundFlag) $<

clean:
	rm *.o
cleanall:
	rm *.mod *.o *.x
run:
	echo input_data.txt | ./$(EXECUTABLE)
runOutput:
	(echo input_data.txt | ./$(EXECUTABLE)> output.txt)&
mpirun:
	srun -n 4 -p pdebug ./$(EXECUTABLE)
# End of the makefile
