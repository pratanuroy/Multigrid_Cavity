module multigrid_levels
use kind_parameters
 implicit none

 type :: level_tag

 integer :: nx
 integer :: ny

 ! Property variables

 real(kind=dp), dimension(:,:), pointer :: rho, mu, kcond, cp, gammat, rho_uf, rho_vf, mu_uf, mu_vf, gammat_uf, gammat_vf

 ! Velocity variables

 real(kind=dp), dimension(:,:), pointer :: u, v, uf, vf
 real(kind=dp), dimension(:,:), pointer :: ap, aw, ae, as, an, bu, bv, b, d_u, d_v, fu, fv, auf_bar, avf_bar, f_uf, f_vf, su, sv, s_uf, s_vf
 real(kind=dp), dimension(:,:), pointer :: resu, resv, res_uf, res_vf

 ! Pressure correction variables

 real(kind=dp), dimension(:,:), pointer :: p, pc, p_uf, p_vf, pc_uf, pc_vf
 real(kind=dp), dimension(:,:), pointer :: app, awp, aep, asp, anp, bp, fp
 real(kind=dp), dimension(:,:), pointer :: respc

  ! Temperature variables

! real(kind=dp), dimension(:,:), pointer :: t
! real(kind=dp), dimension(:,:), pointer :: apt, awt, aet, ast, ant, bt, ft
! real(kind=dp), dimension(:,:), pointer :: rest

 ! For storing previous mg iteration
 real(kind=dp), dimension(:,:), pointer :: uprev, vprev, pcprev, ufprev, vfprev

 ! Geometry variables

 real(kind=dp), dimension(:), pointer :: dx, dy
 real(kind=dp), dimension(:), pointer :: delx, dely
 real(kind=dp), dimension(:), pointer:: xp, yp

! Variables for unsteady calculations

 real(kind=dp), dimension(:,:), pointer:: ap_zero, apt_zero, u0, v0, t0, auf_zero, avf_zero, uf0, vf0, f_uvel, f_vvel

 ! Other variables

 real(kind=dp), dimension(:,:), pointer :: error_u, error_v, error_uf, error_vf

! For parallel variables

integer:: xi, xf

!Variables for level set method
real(kind=dp), dimension(:,:), pointer:: phi_x, phi_y, nv_x, nv_y, grad_phi, lambda, H_phi
real(kind=dp), dimension(:,:), pointer:: T, Tprev, T_exact, T_RK1, T_RK2, alpha_x, alpha_y, phi_sign, kappa, error_check

 end type level_tag

 integer:: num_level
 type(level_tag), dimension(:), allocatable:: mg_level

end module multigrid_levels
