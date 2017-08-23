subroutine prelim_data
        use global_variables
        use parameters
use kind_parameters
        implicit none

        print *, 'Cavity width = ',Lx
        print *, 'Cavity height = ',Ly
        print *, ''
        print *, 'No. of x control volumes = ',icv
        print *, 'No. of y control volumes = ',jcv
        print *, ''

        print *, 'Property flag (1 for Water, 2 for Air) = ' , prop_flag
        print *, 'Reynolds Number = ',reynolds_number
        print *, ''

        print *, 'Boundary Conditions: '
        print *, ''
        print *, 'u-velocity at left boundary: ',u_west
        print *, 'u-velocity at right boundary: ',u_east
        print *, 'u-velocity at bottom boundary: ',u_south
        if (lid_velocity .gt. 0.0d0) then 

                print *, 'u-velocity at top boundary: ',lid_velocity

        else 
                print *, 'u-velocity at top boundary: ',u_north
        end if

        print *, ''
        print *, 'v-velocity at left boundary: ',v_west
        print *, 'v-velocity at right boundary: ',v_east
        print *, 'v-velocity at bottom boundary: ',v_south
        print *, 'v-vleocity at top boundary: ',v_north
        print *, ''
        print *, 'temperature at left boundary: ',t_west
        print *, 'temperature at right boundary: ',t_east
        print *, 'temperature at bottom boundary: ',t_south
        print *, 'temperature at top boundary: ',t_north

        print *, ''
        print *, 'relaxation factor for u = ',alphau
        print *, 'relaxation factor for v = ',alphav
        print *, 'relaxation factor for p = ',alphap
        print *, 'relaxation factor for t = ',alphat
        print *, ''

        print *, 'convergence epsilon for u = ',epsilonu
        print *, 'convergence epsilon for v = ',epsilonv
        print *, 'convergence epsilon for p = ',epsilonp
        print *, 'convergence epsilon for t = ',epsilont
        print *, ''

        print *, 'Maximum Iteration = ',itermax
        print *, 'Program name = ', prog_name
        print *, ''

        print *, 'Is the vtk file writing on or off ?',write_vtk
        print *, 'Is variable property on or off ? ',variable_property
        print *, ''
        !print *, ipmax, jpmax, iumax, jumax, ivmax, jvmax

        print *, 'Reynolds Number = ',reynolds_number
        print *, 'Lid velocity = ',lid_velocity, 'm/s'
        print *, 'Hydraulic Diameter = ', hydraulic_dia, 'm'
        write(*,*) "Rayeligh Number = ", Rayleigh
        print *, ''

        print *, 'Time step = ', time_step
        print *, 'Maximum time = ', max_time
        print *, 'A21 = ', A21
        print *, 'A22 = ', A22
        print *, 'Multigrid Handle = ', mg_handle

        !print *, 'dx = '
        !print *, dx
        !print *, ''
        !print *, 'dy = '
        !print *, dy
        !print *, ''
        !print *, 'delx = '
        !print *, delx
        !print *, ''
        !print *, 'dely = '
        !print *, dely
        !print *, ''
        !
        !print *, 'xp = '
        !print *, xp
        !print *, ''
        !print *, 'yp = '
        !print *, yp
        !print *, ''

end subroutine
