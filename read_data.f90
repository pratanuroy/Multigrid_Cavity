subroutine read_data
        use global_variables
        use parameters
        use multigrid_levels
        use parallel_variables
        use kind_parameters
        implicit none

        logical:: file_exists
        integer:: fstatus = 0

        !read(5,*) input_file

        input_file = 'input_data.txt'

        Inquire(file=input_file,exist = file_exists)
        if (taskid .eq. root) write(*,*) 'Does the input file exist?', file_exists

        if (file_exists .eq. .FALSE.) then

                if (taskid .eq. root) write(*,*)'The input file does not exist. Aborting program!!'
                call exit(fstatus)

        end if

        if (file_exists .and. taskid .eq. root) then

                namelist /geometry_data/ dim_prob, prob_lo, prob_hi, icv, jcv, grid_type
                namelist /property_data/ prop_flag, variable_property

                namelist /boundary_data/ u_west, u_east, u_south, u_north, &
                        v_west, v_east, v_south, v_north, &
                        t_west, t_east, t_south, t_north 

                namelist /parameter_data/ alphau, alphav, alphap, alphat, &
                        epsilonu, epsilonv, epsilonp, epsilont, &
                        itermax, prog_name, write_VTK, time_step, max_time, reynolds_number, gravity, &
                        output_step, monitor_step

                namelist /debug_data / debug, parallel_debug
                namelist /multigrid_data / mg_handle, num_level

                open(unit=8, file=input_file, status = 'old')

                read(8, geometry_data)
                read(8, property_data)
                read(8, boundary_data)
                read(8, parameter_data)
                read(8, debug_data)
                read(8, multigrid_data)

                close(8)

                if (debug) write(*,geometry_data)
                if (debug) write(*,property_data)
                if (debug) write(*,boundary_data)
                if (debug) write(*,parameter_data)

        end if

end subroutine
