module global_variables
use kind_parameters

 !! Geometry variables

 integer*4 :: dim_prob
 real(kind=dp), dimension(1:3):: prob_lo, prob_hi
 real(kind=dp), dimension(:), pointer:: dx, dy, delxu, delyu, delxv, delyv, delx, dely
 real(kind=dp), dimension(:), pointer:: xu, yu, xv, yv, xp, yp
 real(kind=dp) :: xmax, xmin, ymax, ymin

 !! Velocity variables

 real(kind=dp), dimension(:,:), pointer:: u, v, aw, ae, as, an, ap, bu, bv, d_u, d_v, fu, fv, resu, resv, su, sv, s_uf, s_vf
 real(kind=dp), dimension(:,:), pointer:: uf, vf, auf_bar, avf_bar

 real(kind=dp), dimension(:,:), pointer:: p, pc, awp, aep, asp, anp, app, bp, fp
 real(kind=dp), dimension(:,:), pointer:: p_uf, pc_uf, p_vf, pc_vf, res_uf, res_vf, f_uf, f_vf

!real(kind=dp), dimension(:,:), pointer:: t, awt, aet, ast, ant, apt, bt, ft
real(kind=dp), dimension(:,:), pointer:: rho, mu, kcond, cp, gammat
real(kind=dp), dimension(:,:), pointer:: rho_uf, mu_uf, rho_vf, mu_vf, gammat_uf, gammat_vf

real(kind=dp), allocatable, dimension(:,:):: uprint, vprint, pprint, tprint, rho_print, uexact, vexact, u_error, v_error

 ! Variables for unsteady calculations

 real(kind=dp), dimension(:,:), pointer:: ap_zero
 real(kind=dp), dimension(:,:), pointer:: atp_zero, f_uvel, f_vvel, f_temperature, u0, v0, t0, rho0, uf0, vf0

 !! TDMA variables

 real(kind=dp), allocatable:: a(:,:), b(:)

! Constants for Runge-Kutta 2nd order implicit method
!real(kind=dp):: A21 = 0.0d0
!real(kind=dp):: A22 = 1.0d0

!Variables for level set method
!real(kind=dp), dimension(:,:), pointer:: phi_x, phi_y, nv_x, nv_y, lambda, H_phi, grad_phi
real(kind=dp), dimension(:,:), pointer:: T, Tprev
real(kind=dp), allocatable, dimension(:,:):: y_rt

end module global_variables
