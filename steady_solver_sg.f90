 subroutine steady_solver_sg
USE global_Variables
 USE parameters
 USE multigrid_levels
 USE parallel_variables
use kind_parameters
 implicit none
  include "mpif.h"

        call reset_pressure(pc, ipmax, jpmax, start_x, end_x)

        call velocity_coefficients(aw, ae, as, an, ap, ap_zero, d_u, d_v, uf, vf, rho_uf, mu_uf, rho_vf, mu_vf, gammat_uf, gammat_vf, dx, dy, delx, dely, ipmax, jpmax, start_x, end_x)
        call source_terms(u0, v0, uf0, vf0, ap_zero, f_uvel, f_vvel, su, sv, bu, bv, s_uf, s_vf, rho_vf, p_uf, p_vf, dx, dy, ipmax, jpmax, start_x, end_x)
        call line_gs_solver(u, ap, aw, ae, as, an, bu, fu, ipmax, jpmax, alphau, start_x, end_x)
        call line_gs_solver(v, ap, aw, ae, as, an, bv, fv, ipmax, jpmax, alphav, start_x, end_x)

        call  momentum_interpolation(u, v, p, aw, ae, as, an, ap, su, sv, f_uf, f_vf, s_uf, s_vf, res_uf, res_vf, auf_bar, avf_bar, uf, vf, dx, dy, ipmax, jpmax, start_x, end_x)

        call pcoefficients(uf, vf, auf_bar, avf_bar, awp, aep, asp, anp, app, bp, dx, dy, ipmax, jpmax, start_x, end_x)

        call line_gs_solver(pc, app, awp, aep, asp, anp, bp, fp, ipmax, jpmax, alphapc, start_x, end_x)

        call correction(u, v, uf, vf, pc, p, auf_bar, avf_bar, p_uf, p_vf, d_u, d_v, dx, dy, ipmax, jpmax, start_x, end_x, alphap)

        call convergence(u, aw, ae, as, an, ap, bu, fu, dx, dy, ipmax, jpmax, start_x, end_x, residualu)
        call convergence(v, aw, ae, as, an, ap, bv, fv, dx, dy, ipmax, jpmax, start_x, end_x, residualv)
        call convergencep(bp, fp, dx, dy, ipmax, jpmax, start_x, end_x, residualp)
        !call convergencet

end subroutine
