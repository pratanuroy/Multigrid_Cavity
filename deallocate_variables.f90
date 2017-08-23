subroutine deallocate_variables
    use variables
use kind_parameters
    implicit none

    ! Deallocation of geometry variables

    write(*,*) 'Deallocating geometry variables'

    deallocate(dx, dy, dz, delx, dely, delz, xp, yp, zp, zp_par)

    ! Deallocate velocity variables at node points

    deallocate(u, bu, d_u, s_u, Hu)
    deallocate(v, bv, d_v, s_v, Hv)
    deallocate(w, bw, d_w, s_w, Hw)

    ! Deallocate velocity coefficients at node points

    deallocate(aw, ae, as, an, ab, at, ap)

    ! Deallocate velocity variables at CV faces

    deallocate(uf, vf, wf, auf_bar, avf_bar, awf_bar, Huf_bar, Hvf_bar, Hwf_bar) 

    ! Deallocate p-variables at nodal points

    deallocate(p, pc, apw, ape, aps, apn, apb, apt, app, bp)

    ! Deallocate p-variables at CV Faces

    deallocate(p_uf, p_vf, p_wf, pc_uf, pc_vf, pc_wf)
    
    ! Deallocate temperature-variables

    deallocate(temperature, atw, ate, ats, atn, atb, att, atp, bt)

    ! Deallocate concentration-variables

    deallocate(concentration, concentration_mea, concentration_co2)
    deallocate(acw, ace, acs, acn, acb, act, acp)
    deallocate(bc, bc_unsteady, bc_unsteady_co2, bc_unsteady_mea, bc_phys_mea, bc_phys_co2, bc_react_co2, bc_react_mea)
    deallocate(interfacial_diffusion,  mass_transfer_coeff, Sherwood_number)

    ! Deallocate property variables

    deallocate(rho, mu, kcond, cp, gammat) 
    deallocate(rho_uf, rho_vf, rho_wf, mu_uf, mu_vf, mu_wf, gammat_uf,  gammat_vf, gammat_wf)
    deallocate(gammac, gammac_uf, gammac_vf, gammac_wf)
   
    ! Deallocate variables for unsteady calculations

    deallocate(u0, v0, w0, ap_zero, f_uvel, f_vvel, f_wvel)
    deallocate(temperature0, f_temperature, atp_zero)
    deallocate(concentration0, f_concentration)
    deallocate(rho0)
    deallocate(uf0, vf0, wf0)

    ! Deallocate level set variables
 
    deallocate(T, Tprev, T_RK1, T_RK2, Tsurf, T_exact)
    deallocate(alpha_x, alpha_y, alpha_z)
    deallocate(nv_x, nv_y, nv_z, kappa, nv_uf, nv_vf, nv_wf)
    deallocate(grad_phi, phi_sign, H_phi, delta_phival)

    ! Deallocate contact angel variables as a function of z direction
    
    deallocate(theta_x, theta_y, theta_z, theta)

    ! Deallocate velocities for solid particle moving in viscous fluid

    !deallocate(usolid, vsolid, wsolid, usolid0, vsolid0, wsolid0)

    ! Allocate space for boundary condition type
    deallocate(bctype, bctype_conc)

end subroutine
