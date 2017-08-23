subroutine broadcast_data
        USE parameters
        USE multigrid_levels
        USE parallel_variables
        use kind_parameters
        implicit none
        include "mpif.h"

        call MPI_Bcast(Lx, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(Ly, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(icv, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(jcv, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(reynolds_number, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(prop_flag, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

        call MPI_Bcast(u_west, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(u_east, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(u_south, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(u_north, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

        call MPI_Bcast(v_west, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(v_east, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(v_south, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(v_north, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

        call MPI_Bcast(t_west, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(t_east, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(t_south, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(t_north, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

        call MPI_Bcast(alphau, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(alphav, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(alphap, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(alphat, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

        call MPI_Bcast(epsilonu, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(epsilonv, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(epsilonp, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(epsilont, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

        call MPI_Bcast(itermax, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(time_step, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(max_time, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(num_level, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(mg_handle, 2, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)

        call MPI_Bcast(debug, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(parallel_debug, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(write_VTK, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(output_step, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(monitor_step, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(prog_name, 50, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)


end subroutine
