program main

        ! Written by: Pratanu Roy
        ! Start Date: Dec 22, 2013
        ! Last Modified : Aug 22, 2017
        ! This program solves the 2D Navier-Stokes equations
        ! Domain = [0, Lx] x [0, Ly]
        ! Use multigrid technique with parallelization

        use global_variables
        use parameters
        use multigrid_levels
        use parallel_variables
        use kind_parameters
        implicit none
        include "mpif.h"

        integer:: ilevel
        real(kind=dp) :: total_cpu_time

        alphapc = 1.0d0
        alpha_mg = 0.8d0

        call cpu_time(tstart)

        call MPI_INIT(ierr)

        if (ierr .ne. MPI_SUCCESS) print *, 'Error Initializing ...'

        call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

        if (ierr .ne. MPI_SUCCESS) print *, 'Error in number of processors ...'

        call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)

        if (ierr .ne. MPI_SUCCESS) print *, 'Error in taskid ...'

        if (taskid .eq. root) write(*,*) 'Number of processor = ', numtasks

        root = 0

        ! Read data by root processor (processor no. 0)

        call read_data

        ! Set Lx and Ly from the input file

        xmin = prob_lo(1)
        xmax = prob_hi(1)

        ymin = prob_lo(2)
        ymax = prob_hi(2)

        Lx = xmax - xmin
        Ly = ymax - ymin

        ! Broadcast input data to all the processors

        call broadcast_data

        ! Set the max indices for u, v, p, t

        ipmax = icv + 2
        jpmax = jcv + 2

        call setup_mg_variables(ipmax, jpmax, num_level)
        call setup_mg_geometry

        ! For non-periodic or regular boundary condition set the boundary processor index

        inextid = taskid + 1
        iprevid = taskid - 1

        start_x = mg_level(num_level)%xi
        end_x = mg_level(num_level)%xf

        !! For periodic boundary condition
        !! The leftmost processor sends the left boundary values to the right boundary value of rightmost procesor and vice versa
        !
        !if ( taskid .eq. 0) then
        !
        ! iprevid = numtasks-1
        !
        !else if (taskid .eq. numtasks - 1) then
        !
        ! inextid = 0
        !
        !end if
        !
        ! For usual boundary condition
        ! The leftmost processor has no left neighbour and the rightmost processor has no right neighbour

        if ( taskid .eq. 0) then

                iprevid = MPI_PROC_NULL

        else if (taskid .eq. numtasks-1) then

                inextid = MPI_PROC_NULL

        end if

        if (parallel_debug .and. taskid .eq. root) write(*,*) 'start_x ', start_x
        if (parallel_debug .and. taskid .eq. root) write(*,*) 'end_x ', end_x
        if (parallel_debug .and. taskid .eq. root) write(*,*) 'taskid = ', taskid, 'iprevid =', iprevid, 'inextid = ', inextid

        call initial_conditions(num_level)

        if (taskid .eq. root) then

                call prelim_data

        end if

        !tcounter = 0
        total_cpu_time  = 0.0d0

        !do time = 0.0d0, max_time, time_step

        tstart = MPI_WTime()

        residualu = 1.0d0
        residualv = 1.0d0
        residualp = 1.0d0
        residualt = 0.0d0

        is_converged = 0

        do ilevel = 1,num_level

        mg_level(ilevel)%u0(:,:) = mg_level(ilevel)%u(:,:)
        mg_level(ilevel)%v0(:,:) = mg_level(ilevel)%v(:,:)
        !mg_level(ilevel)%t0(:,:) = mg_level(ilevel)%t(:,:)
        !    mg_level(ilevel)%Tprev(:,:) = mg_level(ilevel)%T(:,:)
        mg_level(ilevel)%uf0(:,:) = mg_level(ilevel)%uf(:,:)
        mg_level(ilevel)%vf0(:,:) = mg_level(ilevel)%vf(:,:)


        !!!! This part is for unsteady simulations. Currently turned off
        !!!! ------------------------------------------------------------
        if (1 .eq. 0) then 
        call unsteady_functions(mg_level(ilevel)%ap_zero, &
                mg_level(ilevel)%u0,      &
                mg_level(ilevel)%v0,      &
                mg_level(ilevel)%f_uvel,      &
                mg_level(ilevel)%f_vvel,      &
                mg_level(ilevel)%aw,      &
                mg_level(ilevel)%ae,      &
                mg_level(ilevel)%as,      &
                mg_level(ilevel)%an,      &
                mg_level(ilevel)%p_uf,      &
                mg_level(ilevel)%p_vf,      &
                mg_level(ilevel)%su,      &
                mg_level(ilevel)%sv,      &
                mg_level(ilevel)%fu,      &
                mg_level(ilevel)%fv,      &
                mg_level(ilevel)%rho,      &
                mg_level(ilevel)%dx,      &
                mg_level(ilevel)%dy,      &
                mg_level(ilevel)%nx,      &
                mg_level(ilevel)%ny,      &
                mg_level(ilevel)%xi,      &
                mg_level(ilevel)%xf       )
        end if
        !!!!------------------------------------------------------------
        call properties( mg_level(ilevel)%T,    &
                mg_level(ilevel)%rho,    &
                mg_level(ilevel)%mu,    &
                mg_level(ilevel)%rho_uf,    &
                mg_level(ilevel)%rho_vf,    &
                mg_level(ilevel)%mu_uf,    &
                mg_level(ilevel)%mu_vf,    &
                mg_level(ilevel)%kcond,    &
                mg_level(ilevel)%cp,   &
                mg_level(ilevel)%gammat,   &
                mg_level(ilevel)%nx,   &
                mg_level(ilevel)%ny,   &
                mg_level(ilevel)%xi,     &
                mg_level(ilevel)%xf      )

        end do

        ! Start from the coarses grid for Full Multigrid (FMG) method. 

        ilevel = 1

        do iteration = 1,itermax

        call steady_solver_mg(ilevel)
        call monitor_residuals

        call check_convergence

        if (is_converged) exit

        !    call steady_solver_sg

        end do

        tcounter = tcounter + 1

        if (taskid .eq. 0) then

                write(*,*) 'time step = ', time

        end if

        tend = MPI_WTime()
        if (parallel_debug) write(*,*) 'Total time = ', tend-tstart, 'secs'
        total_cpu_time = total_cpu_time + (tend-tstart)

        !end do

        call write_data(num_level)

        if (taskid .eq. root .and. debug) then

                call print_results(ipmax, jpmax, start_x, end_x, u, v, p, uf, vf, xp, yp)

        end if

        call deallocate_mg_variables()

        if (taskid .eq. root) then 
                write(*,*) 'taskid = ', taskid, 'End of program'
                write(*,*) 'total CPU time = ', total_cpu_time

                !write(*,*) 'Communication time = ', communication_time, 'secs'
                !write(*,*) 'Computation time = ', (tend-tstart) - communication_time, 'secs'

        end if

        call MPI_FINALIZE(ierr)

end program main
