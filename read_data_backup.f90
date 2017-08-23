subroutine read_data
use global_variables
use parameters
use multigrid_levels
use kind_parameters
implicit none

logical:: file_exists
integer:: fstatus = 0

!read(5,*) input_file

input_file = 'input_data.txt'

Inquire(file=input_file,exist = file_exists)
write(*,*) 'Does the input file exist?', file_exists

if (file_exists .eq. .FALSE.) then

    write(*,*)'The input file does not exist. Aborting program!!'
    call exit(fstatus)

end if

open(unit=8, file=input_file, form = 'formatted')

! Read Grid Domain Parameters
read(8,'(4/, E10.6)') Lx
read(8,'(E10.6)') Ly
read(8,'(3/, I4)') icv
read(8,'(I4)') jcv

! Read Flow Properties
read(8,'(3/, I2)') prop_flag
read(8,'(3/, F)') reynolds_number

! Read Boundary Conditions
! Boundary conditions for u-velocity

read(8,'(4/, F12.10)') u_west
read(8,'(F12.10)') u_east
read(8,'(F12.10)') u_south
read(8,'(F12.10)') u_north

! Boundary conditions for v-velocity

read(8,'(2/, F12.10)') v_west
read(8,'(F12.10)') v_east
read(8,'(F12.10)') v_south
read(8,'(F12.10)') v_north

! Boundary conditions for temperature

read(8,'(2/, F12.10)') t_west
read(8,'(F12.10)') t_east
read(8,'(F12.10)') t_south
read(8,'(F12.10)') t_north

! Relaxation Factors

read(8,'(2/, E4.2)') alphau
read(8,'(E4.2)') alphav
read(8,'(E4.2)') alphap
read(8,'(E4.2)') alphat

! Convergence Criteria

read(8,'(2/, E6.1)') epsilonu
read(8,'(E6.1)') epsilonv
read(8,'(F6.1)') epsilonp
read(8,'(F6.1)') epsilont

! Maximum number of iterations for Navier Stokes Solver (steady state / inner loop)
read(8,'(2/, I6)') itermax

! program name (for output file)
read(8,'(2/, A)') prog_name

! does input data exist?
read(8,'(2/, A)') is_input_data
write(*,*) 'does input data (u, v, t, p) exist?', is_input_data

! Writing VTK file (on/off)
read(8,'(2/, A)') write_vtk

! Variable property with temperature (on/off)
read(8,'(2/, A)') variable_property

! Grid Type (Uniform/ Nonuniform)
read(8,'(2/, A)') grid_type

! Value of q; Grid size = 2^q * 2^q
!read(8,'(2/, I)') qval

! Unsteady parameters
! Time Step
read(8,'(3/, E6.1)') time_step

! Maximum Time
read(8,'(2/, E6.1)') max_time

! Maximum Time
read(8,'(2/, A)') mg_handle

! Maximum Time
read(8,'(2/, I)') num_level

close(8)

    if (file_exists .and. taskid .eq. 0) then

    namelist /geometry_data/ dim_prob, prob_lo, prob_hi, icv, jcv, grid_type
    namelist /property_data/ prop_flag, variable_property, &
        rho_l, mu_l, kcond_l, cp_l, &
        rho_g, mu_g, kcond_g, cp_g
        v_west, v_east, v_south, v_north &
        t_west, t_east, t_south, t_north &

    namelist /parameter_data/ alphau, alphav, alphap, alphat &
        epsilonu, epsilonv, epsilonp, epsilont &
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

    close(8)

    if (debug) write(*,geometry_data)
    if (debug) write(*,property_data)
    if (debug) write(*,boundary_data)
    if (debug) write(*,parameter_data)

! Overwrite the original Lx and Ly from input file

xmin = prob_lo(1)
xmax = prob_hi(1)

ymin = prob_lo(2)
ymax = prob_hi(2)

Lx = xmax - xmin
Ly = ymax - ymin

end if

end subroutine
