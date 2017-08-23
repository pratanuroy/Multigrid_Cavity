module parameters
use kind_parameters
 !! Constant Parameters

 real(kind=dp), parameter:: pi = 3.1415926535897932384626433832795028841971693993751058209d0

 !! Iteration variables
 integer:: itermax, counter, maxcase, qval
 integer*4:: icv, jcv, ipmax, jpmax, N, isol, iteration
 integer, allocatable:: iter_num(:)

 !! Other parameters
 real(kind=dp):: Lx, Ly, reynolds_number, inlet_velocity, lid_velocity, hydraulic_dia, totalmb, mass_ratio, v_ref, rho_ref, time_ref, mu_ref, bond_number
 real(kind=dp):: beta, tambient, Rayleigh, alpha, gravity

 !! For specifying boundary conditions
 real(kind=dp):: u_west, u_east, u_south, u_north, v_west, v_east, v_south, v_north
 real(kind=dp):: t_west, t_east, t_south, t_north

!! Property Flag
 integer*4:: prop_flag

!! Under-relaxation parameters, residual variables

 real(kind=dp):: alphau, alphav, alphap, alphat, alphapc, alpha_mg
 real(kind=dp):: epsilonu, epsilonv, epsilonp, epsilont

! Residual Variables
 real(kind=dp):: residualu, residualv, residualp, residualt
 integer*4 :: is_converged

! Monitor convergence and output steps variables

integer*4 :: output_step, monitor_step
logical :: write_VTK

! Debug variables

integer*4 :: debug, parallel_debug

!! CPU Time computing variables

 real(kind=dp):: tstart, tend

  !! Print files variables
 character*50:: filename, string, input_file, prog_name,  grid_type, mg_handle
 logical:: is_input_data, variable_property

!! Unsteady variables

 real(kind=dp):: time_step, max_time, time, tlocal1, tlocal2
 integer:: tcounter

! Constants for Runge-Kutta 2nd order implicit method
real(kind=dp):: A21 = 0.0d0
real(kind=dp):: A22 = 1.0d0

! Parameters for level set method

real(kind=dp):: CFL, epsilon, delta_phi, surface_tension
real(kind=dp) :: rho_g, rho_l, mu_g, mu_l, kcond_g, kcond_l, cp_g, cp_l
real(kind=dp):: radius, center_x, center_y, L_perimeter, global_error, volume_error, error_sum, H_phi_exact, H_phi_numerical, vol_exact, vol_sum, vol_error
real(kind=dp):: alphax1, alphax2, alphax3, alphay1, alphay2, alphay3, phi_xx, phi_xy, phi_yy
real(kind=dp):: tau_redistance, phi_heaviside, phi_delta, gamma_ls, beta_ls
integer:: iter_redistance

end module parameters
