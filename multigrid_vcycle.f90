recursive subroutine multigrid_vcycle(current_level)
USE multigrid_levels
USE parallel_variables
use parameters
use kind_parameters
implicit none

integer:: i, current_level, l
!real(kind=dp), dimension(mg_level(current_level)%nx, mg_level(current_level)%ny) :: error_u, error_v, error_pc
integer::  imax, jmax
!real(kind=dp) :: temp_sum

if (debug) write(*,*) 'Current Level = ', current_level

imax = mg_level(current_level)%nx
jmax = mg_level(current_level)%ny


SELECT CASE(current_level)

CASE(1)

 do i = 1,iter_num(current_level)

    call velocity_coefficients(mg_level(current_level)%aw,      &
                               mg_level(current_level)%ae,      &
                               mg_level(current_level)%as,      &
                               mg_level(current_level)%an,      &
                               mg_level(current_level)%ap,      &
                               mg_level(current_level)%ap_zero, &
                               mg_level(current_level)%d_u,     &
                               mg_level(current_level)%d_v,     &
                               mg_level(current_level)%uf,      &
                               mg_level(current_level)%vf,      &
                               mg_level(current_level)%rho_uf,  &
                               mg_level(current_level)%mu_uf,   &
                               mg_level(current_level)%rho_vf,  &
                               mg_level(current_level)%mu_vf,   &
                               mg_level(current_level)%gammat_uf,  &
                               mg_level(current_level)%gammat_vf,   &
                               mg_level(current_level)%dx,      &
                               mg_level(current_level)%dy,      &
                               mg_level(current_level)%delx,    &
                               mg_level(current_level)%dely,    &
                               mg_level(current_level)%nx,      &
                               mg_level(current_level)%ny,      &
                               mg_level(current_level)%xi,      &
                               mg_level(current_level)%xf       )

    call source_terms(         mg_level(current_level)%u0,      &
                               mg_level(current_level)%v0,      &
                                mg_level(current_level)%uf0,      &
                               mg_level(current_level)%vf0,      &
                               mg_level(current_level)%ap_zero, &
                               mg_level(current_level)%f_uvel,     &
                               mg_level(current_level)%f_vvel,     &
                               mg_level(current_level)%su,      &
                               mg_level(current_level)%sv,      &
                               mg_level(current_level)%bu,      &
                               mg_level(current_level)%bv,  &
                               mg_level(current_level)%s_uf,  &
                               mg_level(current_level)%s_vf,   &
                               mg_level(current_level)%rho_vf,   &
                               mg_level(current_level)%p_uf,  &
                               mg_level(current_level)%p_vf,   &
                               mg_level(current_level)%dx,      &
                               mg_level(current_level)%dy,      &
                               mg_level(current_level)%nx,      &
                               mg_level(current_level)%ny,      &
                               mg_level(current_level)%xi,      &
                               mg_level(current_level)%xf       )

    call line_gs_solver(mg_level(current_level)%u,      &
                        mg_level(current_level)%ap,     &
                        mg_level(current_level)%aw,     &
                        mg_level(current_level)%ae,     &
                        mg_level(current_level)%as,     &
                        mg_level(current_level)%an,     &
                        mg_level(current_level)%bu,     &
                        mg_level(current_level)%fu,     &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        alphau,                         &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )

    call line_gs_solver(mg_level(current_level)%v,      &
                        mg_level(current_level)%ap,     &
                        mg_level(current_level)%aw,     &
                        mg_level(current_level)%ae,     &
                        mg_level(current_level)%as,     &
                        mg_level(current_level)%an,     &
                        mg_level(current_level)%bv,     &
                        mg_level(current_level)%fv,     &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        alphav,                         &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )

    call momentum_interpolation(mg_level(current_level)%u,      &
                                mg_level(current_level)%v,     &
                                mg_level(current_level)%p,     &
                                mg_level(current_level)%aw,     &
                                mg_level(current_level)%ae,     &
                                mg_level(current_level)%as,     &
                                mg_level(current_level)%an,     &
                                mg_level(current_level)%ap,     &
                                mg_level(current_level)%su,     &
                                mg_level(current_level)%sv,     &
                                mg_level(current_level)%f_uf,     &
                                mg_level(current_level)%f_vf,     &
                                mg_level(current_level)%s_uf,     &
                                mg_level(current_level)%s_vf,     &
                                mg_level(current_level)%res_uf,     &
                                mg_level(current_level)%res_vf,     &
                                mg_level(current_level)%auf_bar,     &
                                mg_level(current_level)%avf_bar,     &
                                mg_level(current_level)%uf,     &
                                mg_level(current_level)%vf,     &
                                mg_level(current_level)%dx,     &
                                mg_level(current_level)%dy,     &
                                mg_level(current_level)%nx,     &
                                mg_level(current_level)%ny,     &
                                mg_level(current_level)%xi,     &
                                mg_level(current_level)%xf      )

!call print_results(imax, jmax, mg_level(current_level)%fu, mg_level(current_level)%fv, mg_level(current_level)%fp, mg_level(current_level)%uf, mg_level(current_level)%vf)

    call reset_pressure(mg_level(current_level)%pc, &
                        mg_level(current_level)%nx, &
                        mg_level(current_level)%ny, &
                        mg_level(current_level)%xi, &
                        mg_level(current_level)%xf  )

    call pcoefficients( mg_level(current_level)%uf,      &
                        mg_level(current_level)%vf,      &
                        mg_level(current_level)%auf_bar,      &
                        mg_level(current_level)%avf_bar,      &
                        mg_level(current_level)%awp,      &
                        mg_level(current_level)%aep,      &
                        mg_level(current_level)%asp,      &
                        mg_level(current_level)%anp,      &
                        mg_level(current_level)%app,      &
                        mg_level(current_level)%bp,     &
                        mg_level(current_level)%dx,      &
                        mg_level(current_level)%dy,      &
                        mg_level(current_level)%nx,      &
                        mg_level(current_level)%ny,      &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )
    do l = 1,2
    call line_gs_solver(mg_level(current_level)%pc,      &
                        mg_level(current_level)%app,     &
                        mg_level(current_level)%awp,     &
                        mg_level(current_level)%aep,     &
                        mg_level(current_level)%asp,     &
                        mg_level(current_level)%anp,     &
                        mg_level(current_level)%bp,     &
                        mg_level(current_level)%fp,     &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        alphapc,                         &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )
    end do

!if (1 .eq. 0) then
!    call tcoefficients( mg_level(current_level)%t,      &
!                        mg_level(current_level)%awt,      &
!                        mg_level(current_level)%aet,      &
!                        mg_level(current_level)%ast,      &
!                        mg_level(current_level)%ant,      &
!                        mg_level(current_level)%apt,      &
!                        mg_level(current_level)%bt,     &
!                        mg_level(current_level)%uf,      &
!                        mg_level(current_level)%vf,      &
!                        mg_level(current_level)%rho_uf,      &
!                        mg_level(current_level)%rho_vf,      &
!                        mg_level(current_level)%gammat_uf,      &
!                        mg_level(current_level)%gammat_vf,      &
!                        mg_level(current_level)%dx,      &
!                        mg_level(current_level)%dy,      &
!                        mg_level(current_level)%delx,      &
!                        mg_level(current_level)%dely,      &
!                        mg_level(current_level)%nx,      &
!                        mg_level(current_level)%ny       )
!
!    call line_gs_solver(mg_level(current_level)%t,      &
!                        mg_level(current_level)%apt,     &
!                        mg_level(current_level)%awt,     &
!                        mg_level(current_level)%aet,     &
!                        mg_level(current_level)%ast,     &
!                        mg_level(current_level)%ant,     &
!                        mg_level(current_level)%bt,     &
!                        mg_level(current_level)%ft,     &
!                        mg_level(current_level)%nx,     &
!                        mg_level(current_level)%ny,     &
!                        alphat,                         &
!                        mg_level(current_level)%dx,     &
!                        mg_level(current_level)%dy,     &
!                        mg_level(current_level)%delx,   &
!                        mg_level(current_level)%dely    )
!
!       call convergence(mg_level(current_level)%t,  &
!                    mg_level(current_level)%awt,   &
!                    mg_level(current_level)%aet,   &
!                    mg_level(current_level)%ast,   &
!                    mg_level(current_level)%ant,   &
!                    mg_level(current_level)%apt,   &
!                    mg_level(current_level)%bt,    &
!                    mg_level(current_level)%ft,    &
!                    mg_level(current_level)%dx,   &
!                    mg_level(current_level)%dy,   &
!                    mg_level(current_level)%nx,   &
!                    mg_level(current_level)%ny,   &
!                    residualt                      )
!end if
   ! mg_level(current_level)%pc = mg_level(current_level)%pc + mg_level(current_level)%pc

    call correction(mg_level(current_level)%u,      &
                    mg_level(current_level)%v,      &
                    mg_level(current_level)%uf,     &
                    mg_level(current_level)%vf,     &
                    mg_level(current_level)%pc,     &
                    mg_level(current_level)%p,      &
                    mg_level(current_level)%auf_bar,&
                    mg_level(current_level)%avf_bar,&
                    mg_level(current_level)%p_uf,   &
                    mg_level(current_level)%p_vf,   &
                    mg_level(current_level)%d_u,    &
                    mg_level(current_level)%d_v,    &
                    mg_level(current_level)%dx,     &
                    mg_level(current_level)%dy,     &
                    mg_level(current_level)%nx,     &
                    mg_level(current_level)%ny,     &
                    mg_level(current_level)%xi,     &
                    mg_level(current_level)%xf,     &
                    alphap                           )

   call convergence(mg_level(current_level)%u,  &
                    mg_level(current_level)%aw,   &
                    mg_level(current_level)%ae,   &
                    mg_level(current_level)%as,   &
                    mg_level(current_level)%an,   &
                    mg_level(current_level)%ap,   &
                    mg_level(current_level)%bu,    &
                    mg_level(current_level)%fu,    &
                    mg_level(current_level)%dx,   &
                    mg_level(current_level)%dy,   &
                    mg_level(current_level)%nx,   &
                    mg_level(current_level)%ny,   &
                    mg_level(current_level)%xi,     &
                    mg_level(current_level)%xf,     &
                    residualu                      )

   call convergence(mg_level(current_level)%v,  &
                    mg_level(current_level)%aw,   &
                    mg_level(current_level)%ae,   &
                    mg_level(current_level)%as,   &
                    mg_level(current_level)%an,   &
                    mg_level(current_level)%ap,   &
                    mg_level(current_level)%bv,    &
                    mg_level(current_level)%fv,    &
                    mg_level(current_level)%dx,   &
                    mg_level(current_level)%dy,   &
                    mg_level(current_level)%nx,   &
                    mg_level(current_level)%ny,   &
                    mg_level(current_level)%xi,     &
                    mg_level(current_level)%xf,     &
                    residualv                      )

   call convergencep(mg_level(current_level)%bp,    &
                     mg_level(current_level)%fp,    &
                     mg_level(current_level)%dx,   &
                     mg_level(current_level)%dy,   &
                     mg_level(current_level)%nx,   &
                     mg_level(current_level)%ny,   &
                     mg_level(current_level)%xi,     &
                     mg_level(current_level)%xf,     &
                     residualp                    )

!    if (taskid .eq. root)  then
!
!        write(*,*) current_level, i
!        call monitor_residuals
!
!    end if
 end do

CASE(2:)

  do i = 1,iter_num(current_level)

!	call print_results(imax, jmax, mg_level(current_level)%u, mg_level(current_level)%v, mg_level(current_level)%p, mg_level(current_level)%uf, mg_level(current_level)%vf)


    call velocity_coefficients(mg_level(current_level)%aw,      &
                               mg_level(current_level)%ae,      &
                               mg_level(current_level)%as,      &
                               mg_level(current_level)%an,      &
                               mg_level(current_level)%ap,      &
                               mg_level(current_level)%ap_zero, &
                               mg_level(current_level)%d_u,     &
                               mg_level(current_level)%d_v,     &
                               mg_level(current_level)%uf,      &
                               mg_level(current_level)%vf,      &
                               mg_level(current_level)%rho_uf,  &
                               mg_level(current_level)%mu_uf,   &
                               mg_level(current_level)%rho_vf,  &
                               mg_level(current_level)%mu_vf,   &
                               mg_level(current_level)%gammat_uf,  &
                               mg_level(current_level)%gammat_vf,   &
                               mg_level(current_level)%dx,      &
                               mg_level(current_level)%dy,      &
                               mg_level(current_level)%delx,    &
                               mg_level(current_level)%dely,    &
                               mg_level(current_level)%nx,      &
                               mg_level(current_level)%ny,      &
                               mg_level(current_level)%xi,      &
                               mg_level(current_level)%xf       )

    call source_terms(         mg_level(current_level)%u0,      &
                               mg_level(current_level)%v0,      &
                                mg_level(current_level)%uf0,      &
                               mg_level(current_level)%vf0,      &
                               mg_level(current_level)%ap_zero, &
                               mg_level(current_level)%f_uvel,     &
                               mg_level(current_level)%f_vvel,     &
                               mg_level(current_level)%su,      &
                               mg_level(current_level)%sv,      &
                               mg_level(current_level)%bu,      &
                               mg_level(current_level)%bv,  &
                               mg_level(current_level)%s_uf,  &
                               mg_level(current_level)%s_vf,   &
                               mg_level(current_level)%rho_vf,   &
                               mg_level(current_level)%p_uf,  &
                               mg_level(current_level)%p_vf,   &
                               mg_level(current_level)%dx,      &
                               mg_level(current_level)%dy,      &
                               mg_level(current_level)%nx,      &
                               mg_level(current_level)%ny,      &
                               mg_level(current_level)%xi,      &
                               mg_level(current_level)%xf       )

    call line_gs_solver(mg_level(current_level)%u,      &
                        mg_level(current_level)%ap,     &
                        mg_level(current_level)%aw,     &
                        mg_level(current_level)%ae,     &
                        mg_level(current_level)%as,     &
                        mg_level(current_level)%an,     &
                        mg_level(current_level)%bu,     &
                        mg_level(current_level)%fu,     &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        alphau,                         &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )

    call line_gs_solver(mg_level(current_level)%v,      &
                        mg_level(current_level)%ap,     &
                        mg_level(current_level)%aw,     &
                        mg_level(current_level)%ae,     &
                        mg_level(current_level)%as,     &
                        mg_level(current_level)%an,     &
                        mg_level(current_level)%bv,     &
                        mg_level(current_level)%fv,     &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        alphav,                         &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )

    call momentum_interpolation(mg_level(current_level)%u,      &
                                mg_level(current_level)%v,     &
                                mg_level(current_level)%p,     &
                                mg_level(current_level)%aw,     &
                                mg_level(current_level)%ae,     &
                                mg_level(current_level)%as,     &
                                mg_level(current_level)%an,     &
                                mg_level(current_level)%ap,     &
                                mg_level(current_level)%su,     &
                                mg_level(current_level)%sv,     &
                                mg_level(current_level)%f_uf,     &
                                mg_level(current_level)%f_vf,     &
                                mg_level(current_level)%s_uf,     &
                                mg_level(current_level)%s_vf,     &
                                mg_level(current_level)%res_uf,     &
                                mg_level(current_level)%res_vf,     &
                                mg_level(current_level)%auf_bar,     &
                                mg_level(current_level)%avf_bar,     &
                                mg_level(current_level)%uf,     &
                                mg_level(current_level)%vf,     &
                                mg_level(current_level)%dx,     &
                                mg_level(current_level)%dy,     &
                                mg_level(current_level)%nx,     &
                                mg_level(current_level)%ny,     &
                                mg_level(current_level)%xi,     &
                                mg_level(current_level)%xf      )

!call print_results(imax, jmax, mg_level(current_level)%fu, mg_level(current_level)%fv, mg_level(current_level)%fp, mg_level(current_level)%uf, mg_level(current_level)%vf)

    call reset_pressure(mg_level(current_level)%pc, &
                        mg_level(current_level)%nx, &
                        mg_level(current_level)%ny, &
                        mg_level(current_level)%xi, &
                        mg_level(current_level)%xf  )

    call pcoefficients( mg_level(current_level)%uf,      &
                        mg_level(current_level)%vf,      &
                        mg_level(current_level)%auf_bar,      &
                        mg_level(current_level)%avf_bar,      &
                        mg_level(current_level)%awp,      &
                        mg_level(current_level)%aep,      &
                        mg_level(current_level)%asp,      &
                        mg_level(current_level)%anp,      &
                        mg_level(current_level)%app,      &
                        mg_level(current_level)%bp,     &
                        mg_level(current_level)%dx,      &
                        mg_level(current_level)%dy,      &
                        mg_level(current_level)%nx,      &
                        mg_level(current_level)%ny,      &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )

    do l = 1,2
    call line_gs_solver(mg_level(current_level)%pc,      &
                        mg_level(current_level)%app,     &
                        mg_level(current_level)%awp,     &
                        mg_level(current_level)%aep,     &
                        mg_level(current_level)%asp,     &
                        mg_level(current_level)%anp,     &
                        mg_level(current_level)%bp,     &
                        mg_level(current_level)%fp,     &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        alphapc,                         &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )

     end do
     
!if (1 .eq. 0) then
!    call tcoefficients( mg_level(current_level)%t,      &
!                        mg_level(current_level)%awt,      &
!                        mg_level(current_level)%aet,      &
!                        mg_level(current_level)%ast,      &
!                        mg_level(current_level)%ant,      &
!                        mg_level(current_level)%apt,      &
!                        mg_level(current_level)%bt,     &
!                        mg_level(current_level)%uf,      &
!                        mg_level(current_level)%vf,      &
!                        mg_level(current_level)%rho_uf,      &
!                        mg_level(current_level)%rho_vf,      &
!                        mg_level(current_level)%gammat_uf,      &
!                        mg_level(current_level)%gammat_vf,      &
!                        mg_level(current_level)%dx,      &
!                        mg_level(current_level)%dy,      &
!                        mg_level(current_level)%delx,      &
!                        mg_level(current_level)%dely,      &
!                        mg_level(current_level)%nx,      &
!                        mg_level(current_level)%ny       )
!
!    call line_gs_solver(mg_level(current_level)%t,      &
!                        mg_level(current_level)%apt,     &
!                        mg_level(current_level)%awt,     &
!                        mg_level(current_level)%aet,     &
!                        mg_level(current_level)%ast,     &
!                        mg_level(current_level)%ant,     &
!                        mg_level(current_level)%bt,     &
!                        mg_level(current_level)%ft,     &
!                        mg_level(current_level)%nx,     &
!                        mg_level(current_level)%ny,     &
!                        alphat                         )
!
!       call convergence(mg_level(current_level)%t,  &
!                    mg_level(current_level)%awt,   &
!                    mg_level(current_level)%aet,   &
!                    mg_level(current_level)%ast,   &
!                    mg_level(current_level)%ant,   &
!                    mg_level(current_level)%apt,   &
!                    mg_level(current_level)%bt,    &
!                    mg_level(current_level)%ft,    &
!                    mg_level(current_level)%dx,   &
!                    mg_level(current_level)%dy,   &
!                    mg_level(current_level)%nx,   &
!                    mg_level(current_level)%ny,   &
!                    residualt                      )
!end if
   ! mg_level(current_level)%pc = mg_level(current_level)%pc + mg_level(current_level)%pc

    call correction(mg_level(current_level)%u,      &
                    mg_level(current_level)%v,      &
                    mg_level(current_level)%uf,     &
                    mg_level(current_level)%vf,     &
                    mg_level(current_level)%pc,     &
                    mg_level(current_level)%p,      &
                    mg_level(current_level)%auf_bar,&
                    mg_level(current_level)%avf_bar,&
                    mg_level(current_level)%p_uf,   &
                    mg_level(current_level)%p_vf,   &
                    mg_level(current_level)%d_u,    &
                    mg_level(current_level)%d_v,    &
                    mg_level(current_level)%dx,     &
                    mg_level(current_level)%dy,     &
                    mg_level(current_level)%nx,     &
                    mg_level(current_level)%ny,     &
                    mg_level(current_level)%xi,     &
                    mg_level(current_level)%xf,     &
                    alphap                           )

   call convergence(mg_level(current_level)%u,  &
                    mg_level(current_level)%aw,   &
                    mg_level(current_level)%ae,   &
                    mg_level(current_level)%as,   &
                    mg_level(current_level)%an,   &
                    mg_level(current_level)%ap,   &
                    mg_level(current_level)%bu,    &
                    mg_level(current_level)%fu,    &
                    mg_level(current_level)%dx,   &
                    mg_level(current_level)%dy,   &
                    mg_level(current_level)%nx,   &
                    mg_level(current_level)%ny,   &
                    mg_level(current_level)%xi,     &
                    mg_level(current_level)%xf,     &
                    residualu                      )

   call convergence(mg_level(current_level)%v,  &
                    mg_level(current_level)%aw,   &
                    mg_level(current_level)%ae,   &
                    mg_level(current_level)%as,   &
                    mg_level(current_level)%an,   &
                    mg_level(current_level)%ap,   &
                    mg_level(current_level)%bv,    &
                    mg_level(current_level)%fv,    &
                    mg_level(current_level)%dx,   &
                    mg_level(current_level)%dy,   &
                    mg_level(current_level)%nx,   &
                    mg_level(current_level)%ny,   &
                    mg_level(current_level)%xi,     &
                    mg_level(current_level)%xf,     &
                    residualv                      )

   call convergencep(mg_level(current_level)%bp,    &
                     mg_level(current_level)%fp,    &
                     mg_level(current_level)%dx,   &
                     mg_level(current_level)%dy,   &
                     mg_level(current_level)%nx,   &
                     mg_level(current_level)%ny,   &
                     mg_level(current_level)%xi,     &
                     mg_level(current_level)%xf,     &
                     residualp                    )

!    if (taskid .eq. root)  then
!
!        write(*,*) current_level, i
!        call monitor_residuals
!
!    end if

 end do

!call print_results(imax, jmax, mg_level(current_level)%u, mg_level(current_level)%v, mg_level(current_level)%p, mg_level(current_level)%uf, mg_level(current_level)%vf)

 call mg_residual(  mg_level(current_level)%u, &
                    mg_level(current_level)%ap,  &
                    mg_level(current_level)%aw,  &
                    mg_level(current_level)%ae,  &
                    mg_level(current_level)%as,  &
                    mg_level(current_level)%an,  &
                    mg_level(current_level)%bu,   &
                    mg_level(current_level)%fu,   &
                    mg_level(current_level)%nx,  &
                    mg_level(current_level)%ny,  &
                    mg_level(current_level)%xi,  &
                    mg_level(current_level)%xf,  &
                    mg_level(current_level)%resu )

 call mg_residual(  mg_level(current_level)%v, &
                    mg_level(current_level)%ap,  &
                    mg_level(current_level)%aw,  &
                    mg_level(current_level)%ae,  &
                    mg_level(current_level)%as,  &
                    mg_level(current_level)%an,  &
                    mg_level(current_level)%bv,   &
                    mg_level(current_level)%fv,   &
                    mg_level(current_level)%nx,  &
                    mg_level(current_level)%ny,  &
                    mg_level(current_level)%xi,  &
                    mg_level(current_level)%xf,  &
                    mg_level(current_level)%resv )

 !call mg_residual(  mg_level(current_level)%t, &
 !                   mg_level(current_level)%apt,  &
 !                   mg_level(current_level)%awt,  &
 !                   mg_level(current_level)%aet,  &
 !!                   mg_level(current_level)%ast,  &
 !                   mg_level(current_level)%ant,  &
 !                   mg_level(current_level)%bt,   &
 !                   mg_level(current_level)%ft,   &
 !                   mg_level(current_level)%nx,  &
 !!                   mg_level(current_level)%ny,  &
 !                  mg_level(current_level)%dx,  &
 !                   mg_level(current_level)%dy,  &
 !                   mg_level(current_level)%rest )

  call fine_to_coarse_res(  mg_level(current_level)%resu ,     &
                        mg_level(current_level)%nx,         &
                        mg_level(current_level)%ny,         &
                        mg_level(current_level)%xi,         &
                        mg_level(current_level)%xf,         &
                        mg_level(current_level-1)%xi,         &
                        mg_level(current_level-1)%xf,         &
                        mg_level(current_level-1)%resu      )

  call fine_to_coarse_res(  mg_level(current_level)%resv,        &
                        mg_level(current_level)%nx,         &
                        mg_level(current_level)%ny,         &
                        mg_level(current_level)%xi,         &
                        mg_level(current_level)%xf,         &
                        mg_level(current_level-1)%xi,         &
                        mg_level(current_level-1)%xf,         &
                        mg_level(current_level-1)%resv      )

  !call fine_to_coarse_res(  mg_level(current_level)%rest,        &
  !                      mg_level(current_level)%nx,         &
  !                      mg_level(current_level)%ny,         &
  !                      mg_level(current_level-1)%rest      )
!mg_level(current_level)%u = mg_level(current_level)%u + 4.0d0*mg_level(current_level)%resu
!mg_level(current_level)%v = mg_level(current_level)%v + 4.0d0*mg_level(current_level)%resv

  call fine_to_coarse(  mg_level(current_level)%u ,     &
                        mg_level(current_level)%nx,         &
                        mg_level(current_level)%ny,         &
                        mg_level(current_level)%xi,         &
                        mg_level(current_level)%xf,         &
                        mg_level(current_level-1)%xi,         &
                        mg_level(current_level-1)%xf,         &
                        mg_level(current_level-1)%u     )

  call fine_to_coarse(  mg_level(current_level)%v,        &
                        mg_level(current_level)%nx,         &
                        mg_level(current_level)%ny,         &
                        mg_level(current_level)%xi,         &
                        mg_level(current_level)%xf,         &
                        mg_level(current_level-1)%xi,         &
                        mg_level(current_level-1)%xf,         &
                        mg_level(current_level-1)%v      )
!mg_level(current_level-1)%u = mg_level(current_level-1)%u + mg_level(current_level-1)%resu
!mg_level(current_level-1)%v = mg_level(current_level-1)%v + mg_level(current_level-1)%resv

  !call fine_to_coarse(  mg_level(current_level)%t,        &
  !                      mg_level(current_level)%nx,         &
  !                      mg_level(current_level)%ny,         &
  !                      mg_level(current_level-1)%t      )

  call fine_to_coarse_uface(  mg_level(current_level)%uf,        &
                        mg_level(current_level)%nx-1,         &
                        mg_level(current_level)%ny,         &
                        mg_level(current_level)%xi,         &
                        mg_level(current_level)%xf,         &
                        mg_level(current_level-1)%xi,         &
                        mg_level(current_level-1)%xf,         &
                        mg_level(current_level-1)%uf      )

  call fine_to_coarse_vface(  mg_level(current_level)%vf,        &
                        mg_level(current_level)%nx,         &
                        mg_level(current_level)%ny-1,         &
                        mg_level(current_level)%xi,         &
                        mg_level(current_level)%xf,         &
                        mg_level(current_level-1)%xi,         &
                        mg_level(current_level-1)%xf,         &
                        mg_level(current_level-1)%vf      )

!write(*,*) 'mg vcycle ', current_level, mg_level(current_level)%xi, mg_level(current_level)%xf

  call fine_to_coarse_uface_sum(  mg_level(current_level)%res_uf,        &
                        mg_level(current_level)%nx-1,         &
                        mg_level(current_level)%ny,         &
                        mg_level(current_level)%xi,         &
                        mg_level(current_level)%xf,         &
                        mg_level(current_level-1)%xi,         &
                        mg_level(current_level-1)%xf,         &
                        mg_level(current_level-1)%res_uf      )

  call fine_to_coarse_vface_sum(  mg_level(current_level)%res_vf,        &
                        mg_level(current_level)%nx,         &
                        mg_level(current_level)%ny-1,         &
                        mg_level(current_level)%xi,         &
                        mg_level(current_level)%xf,         &
                        mg_level(current_level-1)%xi,         &
                        mg_level(current_level-1)%xf,         &
                        mg_level(current_level-1)%res_vf      )
!mg_level(current_level-1)%uf = mg_level(current_level-1)%uf + mg_level(current_level-1)%res_uf
!mg_level(current_level-1)%vf = mg_level(current_level-1)%vf + mg_level(current_level-1)%res_vf

!  call fine_to_coarse_uface(  mg_level(current_level)%res_uf,        &
!                        mg_level(current_level)%nx-1,         &
!                        mg_level(current_level)%ny,         &
!                        mg_level(current_level-1)%res_uf      )
!
!  call fine_to_coarse_vface(  mg_level(current_level)%res_vf,        &
!                        mg_level(current_level)%nx,         &
!                        mg_level(current_level)%ny-1,         &
!                        mg_level(current_level-1)%res_vf      )

    call pcoefficients( mg_level(current_level-1)%uf,      &
                        mg_level(current_level-1)%vf,      &
                        mg_level(current_level-1)%auf_bar,      &
                        mg_level(current_level-1)%avf_bar,      &
                        mg_level(current_level-1)%awp,      &
                        mg_level(current_level-1)%aep,      &
                        mg_level(current_level-1)%asp,      &
                        mg_level(current_level-1)%anp,      &
                        mg_level(current_level-1)%app,      &
                        mg_level(current_level-1)%fp,     &
                        mg_level(current_level-1)%dx,      &
                        mg_level(current_level-1)%dy,      &
                        mg_level(current_level-1)%nx,      &
                        mg_level(current_level-1)%ny,      &
                        mg_level(current_level-1)%xi,     &
                        mg_level(current_level-1)%xf      )

call exchange_data(mg_level(current_level-1)%fp,mg_level(current_level-1)%nx, mg_level(current_level-1)%xi, mg_level(current_level-1)%xf)

    mg_level(current_level-1)%p(:,:) = 0.0d0
    mg_level(current_level-1)%p_uf(:,:) = 0.0d0
    mg_level(current_level-1)%p_vf(:,:) = 0.0d0

!       call properties( mg_level(current_level-1)%T,    &
!                        mg_level(current_level-1)%rho,    &
!                        mg_level(current_level-1)%mu,    &
!                        mg_level(current_level-1)%rho_uf,    &
!                        mg_level(current_level-1)%rho_vf,    &
!                        mg_level(current_level-1)%mu_uf,    &
!                        mg_level(current_level-1)%mu_vf,    &
!                        mg_level(current_level-1)%kcond,    &
!                        mg_level(current_level-1)%cp,   &
!                        mg_level(current_level-1)%gammat,   &
!                        mg_level(current_level-1)%nx,   &
!                        mg_level(current_level-1)%ny,   &
!                        mg_level(current_level-1)%xi,     &
!                        mg_level(current_level-1)%xf      )

    call velocity_coefficients(mg_level(current_level-1)%aw,      &
                               mg_level(current_level-1)%ae,      &
                               mg_level(current_level-1)%as,      &
                               mg_level(current_level-1)%an,      &
                               mg_level(current_level-1)%ap,      &
                               mg_level(current_level-1)%ap_zero, &
                               mg_level(current_level-1)%d_u,     &
                               mg_level(current_level-1)%d_v,     &
                               mg_level(current_level-1)%uf,      &
                               mg_level(current_level-1)%vf,      &
                               mg_level(current_level-1)%rho_uf,  &
                               mg_level(current_level-1)%mu_uf,   &
                               mg_level(current_level-1)%rho_vf,  &
                               mg_level(current_level-1)%mu_vf,   &
                               mg_level(current_level-1)%gammat_uf,  &
                               mg_level(current_level-1)%gammat_vf,   &
                               mg_level(current_level-1)%dx,      &
                               mg_level(current_level-1)%dy,      &
                               mg_level(current_level-1)%delx,    &
                               mg_level(current_level-1)%dely,    &
                               mg_level(current_level-1)%nx,      &
                               mg_level(current_level-1)%ny,      &
                               mg_level(current_level-1)%xi,      &
                               mg_level(current_level-1)%xf     )

    call source_terms(         mg_level(current_level-1)%u0,      &
                               mg_level(current_level-1)%v0,      &
                                mg_level(current_level-1)%uf0,      &
                               mg_level(current_level-1)%vf0,      &
                               mg_level(current_level-1)%ap_zero, &
                               mg_level(current_level-1)%f_uvel,     &
                               mg_level(current_level-1)%f_vvel,     &
                               mg_level(current_level-1)%su,      &
                               mg_level(current_level-1)%sv,      &
                               mg_level(current_level-1)%bu,      &
                               mg_level(current_level-1)%bv,  &
                               mg_level(current_level-1)%s_uf,  &
                               mg_level(current_level-1)%s_vf,   &
                               mg_level(current_level-1)%rho_vf,   &
                               mg_level(current_level-1)%p_uf,  &
                               mg_level(current_level-1)%p_vf,   &
                               mg_level(current_level-1)%dx,      &
                               mg_level(current_level-1)%dy,      &
                               mg_level(current_level-1)%nx,      &
                               mg_level(current_level-1)%ny,      &
                               mg_level(current_level-1)%xi,      &
                               mg_level(current_level-1)%xf       )

!  mg_level(current_level-1)%su = 0.0d0
!  mg_level(current_level-1)%sv = 0.0d0

!    mg_level(current_level-1)%res_uf = 0.0d0
!   mg_level(current_level-1)%res_vf = 0.0d0

! Calculate the source terms for next coarser level

    call mg_residual(  mg_level(current_level-1)%u, &
                    mg_level(current_level-1)%ap,  &
                    mg_level(current_level-1)%aw,  &
                    mg_level(current_level-1)%ae,  &
                    mg_level(current_level-1)%as,  &
                    mg_level(current_level-1)%an,  &
                    mg_level(current_level-1)%su,   &
                    mg_level(current_level-1)%resu,   &
                    mg_level(current_level-1)%nx,  &
                    mg_level(current_level-1)%ny,  &
                    mg_level(current_level-1)%xi,  &
                    mg_level(current_level-1)%xf,  &
                    mg_level(current_level-1)%fu )

    call mg_residual( mg_level(current_level-1)%v, &
                    mg_level(current_level-1)%ap,  &
                    mg_level(current_level-1)%aw,  &
                    mg_level(current_level-1)%ae,  &
                    mg_level(current_level-1)%as,  &
                    mg_level(current_level-1)%an,  &
                    mg_level(current_level-1)%sv,   &
                    mg_level(current_level-1)%resv,   &
                    mg_level(current_level-1)%nx,  &
                    mg_level(current_level-1)%ny,  &
                    mg_level(current_level-1)%xi,  &
                    mg_level(current_level-1)%xf,  &
                    mg_level(current_level-1)%fv )

!    call mg_source( mg_level(current_level-1)%t, &
!                    mg_level(current_level-1)%apt,  &
!                    mg_level(current_level-1)%awt,  &
!                    mg_level(current_level-1)%aet,  &
!                    mg_level(current_level-1)%ast,  &
!                    mg_level(current_level-1)%ant,  &
!                    mg_level(current_level-1)%bt,   &
!                    mg_level(current_level-1)%rest,   &
!                    mg_level(current_level-1)%nx,  &
!                    mg_level(current_level-1)%ny,  &
!                    mg_level(current_level-1)%dx,  &
!                    mg_level(current_level-1)%dy,  &
!                    mg_level(current_level-1)%ft )

!call print_results(mg_level(current_level-1)%nx, mg_level(current_level-1)%ny, mg_level(current_level-1)%resu, mg_level(current_level-1)%resv, mg_level(current_level-1)%fu, mg_level(current_level-1)%uf, mg_level(current_level-1)%vf)

    call momentum_interpolation_mg( mg_level(current_level-1)%u,      &
                                    mg_level(current_level-1)%v,     &
                                    mg_level(current_level-1)%p,     &
                                    mg_level(current_level-1)%su,     &
                                    mg_level(current_level-1)%sv,     &
                                    mg_level(current_level-1)%s_uf,     &
                                    mg_level(current_level-1)%s_vf,     &
                                    mg_level(current_level-1)%f_uf,     &
                                    mg_level(current_level-1)%f_vf,     &
                                    mg_level(current_level-1)%res_uf,     &
                                    mg_level(current_level-1)%res_vf,     &
                                    mg_level(current_level-1)%aw,     &
                                    mg_level(current_level-1)%ae,     &
                                    mg_level(current_level-1)%as,     &
                                    mg_level(current_level-1)%an,     &
                                    mg_level(current_level-1)%ap,     &
                                    mg_level(current_level-1)%uf,     &
                                    mg_level(current_level-1)%vf,     &
                                    mg_level(current_level-1)%dx,     &
                                    mg_level(current_level-1)%dy,     &
                                    mg_level(current_level-1)%nx,     &
                                    mg_level(current_level-1)%ny,     &
                                    mg_level(current_level-1)%xi,     &
                                    mg_level(current_level-1)%xf	  )

!    mg_level(current_level-1)%p = 0.0d0
    mg_level(current_level-1)%fp = -1.0d0*mg_level(current_level-1)%fp
!    mg_level(current_level-1)%fp = 0.0d0


 mg_level(current_level-1)%uprev = mg_level(current_level-1)%u
 mg_level(current_level-1)%vprev = mg_level(current_level-1)%v

!  call print_results(imax/2+1, jmax/2+1, mg_level(current_level-1)%fu, mg_level(current_level-1)%fv, mg_level(current_level-1)%fp, mg_level(current_level-1)%ufprev, mg_level(current_level-1)%vfprev)

!do jp = 2,mg_level(current_level-1)%ny-1
!
! do ip = 2,mg_level(current_level-1)%nx-1
!
!	mg_level(current_level-1)%fp = -1.0d0*((mg_level(current_level-1)%uf(ip-1,jp)-mg_level(current_level-1)%uf(ip,jp))/mg_level(current_level-1)%dx(ip) + (mg_level(current_level-1)%vf(ip,jp-1)-mg_level(current_level-1)%vf(ip,jp))/mg_level(current_level-1)%dy(jp))
!
! end do
!
!end do

  call multigrid_vcycle(current_level-1)

  call mg_error(mg_level(current_level-1)%uprev, &
                mg_level(current_level-1)%u,     &
                mg_level(current_level-1)%nx,    &
                mg_level(current_level-1)%ny,    &
                mg_level(current_level-1)%xi,    &
                mg_level(current_level-1)%xf,    &
                mg_level(current_level-1)%error_u)

  call mg_error(mg_level(current_level-1)%vprev, &
                mg_level(current_level-1)%v,     &
                mg_level(current_level-1)%nx,    &
                mg_level(current_level-1)%ny,    &
                mg_level(current_level-1)%xi,    &
                mg_level(current_level-1)%xf,    &
                mg_level(current_level-1)%error_v)

  !call mg_error(mg_level(current_level-1)%tprev, &
  !              mg_level(current_level-1)%t,     &
  !              mg_level(current_level-1)%nx,    &
  !              mg_level(current_level-1)%ny,    &
  !              mg_level(current_level-1)%error_t)

!  call print_results(imax/2+1, jmax/2+1, mg_level(current_level-1)%fu, mg_level(current_level-1)%fv, mg_level(current_level-1)%fp, mg_level(current_level-1)%uf, mg_level(current_level-1)%vf)
!call print_results(imax, jmax, mg_level(current_level)%bu, mg_level(current_level)%bv, mg_level(current_level)%bp, mg_level(current_level)%uf, mg_level(current_level)%vf)

  call coarse_to_fine(  mg_level(current_level-1)%error_u,  &
                        mg_level(current_level-1)%nx,       &
                        mg_level(current_level-1)%ny,       &
                        mg_level(current_level-1)%xi,         &
                        mg_level(current_level-1)%xf,         &
                        mg_level(current_level)%xi,         &
                        mg_level(current_level)%xf,         &
                        mg_level(current_level)%error_u     )

  call coarse_to_fine(  mg_level(current_level-1)%error_v,  &
                        mg_level(current_level-1)%nx,       &
                        mg_level(current_level-1)%ny,       &
                        mg_level(current_level-1)%xi,         &
                        mg_level(current_level-1)%xf,         &
                        mg_level(current_level)%xi,         &
                        mg_level(current_level)%xf,         &
                        mg_level(current_level)%error_v     )

  call coarse_to_fine(  mg_level(current_level-1)%p,  &
                        mg_level(current_level-1)%nx,       &
                        mg_level(current_level-1)%ny,       &
                        mg_level(current_level-1)%xi,         &
                        mg_level(current_level-1)%xf,         &
                        mg_level(current_level)%xi,         &
                        mg_level(current_level)%xf,         &
                        mg_level(current_level)%pc     )
!
!      call coarse_to_fine_ls(mg_level(current_level-1)%T, &
!                    mg_level(current_level-1)%nx,    &
!                    mg_level(current_level-1)%ny,    &
!                    mg_level(current_level-1)%xi,    &
!                    mg_level(current_level-1)%xf,    &
!                    mg_level(current_level)%xi,    &
!                    mg_level(current_level)%xf,    &
!                    mg_level(current_level)%T,    &
!                    mg_level(current_level)%delx,    &
!                    mg_level(current_level)%dely  )
!

 ! call coarse_to_fine(  mg_level(current_level-1)%error_t,  &
 !                       mg_level(current_level-1)%nx,       &
 !                       mg_level(current_level-1)%ny,       &
 !                       mg_level(current_level)%error_t     )

if (1 .eq. 1) then

  call mg_correction(   mg_level(current_level)%u,      &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        mg_level(current_level)%xi,    &
                        mg_level(current_level)%xf,    &
                        mg_level(current_level)%error_u,&
                        alpha_mg )

  call mg_correction(   mg_level(current_level)%v,      &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        mg_level(current_level)%xi,    &
                        mg_level(current_level)%xf,    &
                        mg_level(current_level)%error_v, &
                        alpha_mg  )

  call mg_correction(   mg_level(current_level)%p,      &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        mg_level(current_level)%xi,    &
                        mg_level(current_level)%xf,    &
                        mg_level(current_level)%pc,     &
                        alpha_mg  )

 ! call mg_correction(   mg_level(current_level)%t,      &
 !                       mg_level(current_level)%nx,     &
 !                       mg_level(current_level)%ny,     &
 !                       mg_level(current_level)%error_t,     &
 !                       alphat  )

end if
!  call mg_correction_face_vel(  mg_level(current_level)%uf,      &
!  								mg_level(current_level)%vf,      &
!  								mg_level(current_level)%pc,      &
!  								mg_level(current_level)%auf_bar,      &
!  								mg_level(current_level)%avf_bar,      &
!`								mg_level(current_level)%dx,      &
!								mg_level(current_level)%dy,      &
!                        		mg_level(current_level)%nx,     &
 !                       		mg_level(current_level)%ny)


!call print_results(imax, jmax, mg_level(current_level)%fu, mg_level(current_level)%fv, mg_level(current_level)%fp, mg_level(current_level)%error_uf, mg_level(current_level)%error_vf)

	call face_value(   mg_level(current_level)%p,      &
                        mg_level(current_level)%p_uf,     &
                        mg_level(current_level)%p_vf,     &
                        mg_level(current_level)%dx,     &
                        mg_level(current_level)%dy,     &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        mg_level(current_level)%xi,    &
                        mg_level(current_level)%xf      )

 do i = 1,iter_num(current_level)

!
!        do jp = mg_level(current_level)%xi,mg_level(current_level)%xf
!
!            do ip = 2,mg_level(current_level)%nx-1
!
!                mg_level(current_level)%alpha_x(ip,jp) = 0.50d0*(mg_level(current_level)%u(ip,jp) + mg_level(current_level)%u0(ip,jp))*time_step
!                mg_level(current_level)%alpha_y(ip,jp) = 0.50d0*(mg_level(current_level)%v(ip,jp) + mg_level(current_level)%v0(ip,jp))*time_step
!
!            end do
!
!        end do
!

    call velocity_coefficients(mg_level(current_level)%aw,      &
                               mg_level(current_level)%ae,      &
                               mg_level(current_level)%as,      &
                               mg_level(current_level)%an,      &
                               mg_level(current_level)%ap,      &
                               mg_level(current_level)%ap_zero, &
                               mg_level(current_level)%d_u,     &
                               mg_level(current_level)%d_v,     &
                               mg_level(current_level)%uf,      &
                               mg_level(current_level)%vf,      &
                               mg_level(current_level)%rho_uf,  &
                               mg_level(current_level)%mu_uf,   &
                               mg_level(current_level)%rho_vf,  &
                               mg_level(current_level)%mu_vf,   &
                               mg_level(current_level)%gammat_uf,  &
                               mg_level(current_level)%gammat_vf,   &
                               mg_level(current_level)%dx,      &
                               mg_level(current_level)%dy,      &
                               mg_level(current_level)%delx,    &
                               mg_level(current_level)%dely,    &
                               mg_level(current_level)%nx,      &
                               mg_level(current_level)%ny,      &
                               mg_level(current_level)%xi,      &
                               mg_level(current_level)%xf       )

    call source_terms(         mg_level(current_level)%u0,      &
                               mg_level(current_level)%v0,      &
                               mg_level(current_level)%uf0,      &
                               mg_level(current_level)%vf0,      &
                               mg_level(current_level)%ap_zero, &
                               mg_level(current_level)%f_uvel,     &
                               mg_level(current_level)%f_vvel,     &
                               mg_level(current_level)%su,      &
                               mg_level(current_level)%sv,      &
                               mg_level(current_level)%bu,      &
                               mg_level(current_level)%bv,  &
                               mg_level(current_level)%s_uf,  &
                               mg_level(current_level)%s_vf,   &
                               mg_level(current_level)%rho_vf,   &
                               mg_level(current_level)%p_uf,  &
                               mg_level(current_level)%p_vf,   &
                               mg_level(current_level)%dx,      &
                               mg_level(current_level)%dy,      &
                               mg_level(current_level)%nx,      &
                               mg_level(current_level)%ny,      &
                               mg_level(current_level)%xi,      &
                               mg_level(current_level)%xf       )

    call line_gs_solver(mg_level(current_level)%u,      &
                        mg_level(current_level)%ap,     &
                        mg_level(current_level)%aw,     &
                        mg_level(current_level)%ae,     &
                        mg_level(current_level)%as,     &
                        mg_level(current_level)%an,     &
                        mg_level(current_level)%bu,     &
                        mg_level(current_level)%fu,     &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        alphau,                         &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )

    call line_gs_solver(mg_level(current_level)%v,      &
                        mg_level(current_level)%ap,     &
                        mg_level(current_level)%aw,     &
                        mg_level(current_level)%ae,     &
                        mg_level(current_level)%as,     &
                        mg_level(current_level)%an,     &
                        mg_level(current_level)%bv,     &
                        mg_level(current_level)%fv,     &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        alphav,                         &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )

    call momentum_interpolation(mg_level(current_level)%u,      &
                                mg_level(current_level)%v,     &
                                mg_level(current_level)%p,     &
                                mg_level(current_level)%aw,     &
                                mg_level(current_level)%ae,     &
                                mg_level(current_level)%as,     &
                                mg_level(current_level)%an,     &
                                mg_level(current_level)%ap,     &
                                mg_level(current_level)%su,     &
                                mg_level(current_level)%sv,     &
                                mg_level(current_level)%f_uf,     &
                                mg_level(current_level)%f_vf,     &
                                mg_level(current_level)%s_uf,     &
                                mg_level(current_level)%s_vf,     &
                                mg_level(current_level)%res_uf,     &
                                mg_level(current_level)%res_vf,     &
                                mg_level(current_level)%auf_bar,     &
                                mg_level(current_level)%avf_bar,     &
                                mg_level(current_level)%uf,     &
                                mg_level(current_level)%vf,     &
                                mg_level(current_level)%dx,     &
                                mg_level(current_level)%dy,     &
                                mg_level(current_level)%nx,     &
                                mg_level(current_level)%ny,     &
                                mg_level(current_level)%xi,     &
                                mg_level(current_level)%xf      )

!call print_results(imax, jmax, mg_level(current_level)%fu, mg_level(current_level)%fv, mg_level(current_level)%fp, mg_level(current_level)%uf, mg_level(current_level)%vf)

    call reset_pressure(mg_level(current_level)%pc, &
                        mg_level(current_level)%nx, &
                        mg_level(current_level)%ny, &
                        mg_level(current_level)%xi, &
                        mg_level(current_level)%xf  )

    call pcoefficients( mg_level(current_level)%uf,      &
                        mg_level(current_level)%vf,      &
                        mg_level(current_level)%auf_bar,      &
                        mg_level(current_level)%avf_bar,      &
                        mg_level(current_level)%awp,      &
                        mg_level(current_level)%aep,      &
                        mg_level(current_level)%asp,      &
                        mg_level(current_level)%anp,      &
                        mg_level(current_level)%app,      &
                        mg_level(current_level)%bp,     &
                        mg_level(current_level)%dx,      &
                        mg_level(current_level)%dy,      &
                        mg_level(current_level)%nx,      &
                        mg_level(current_level)%ny,      &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )
    do l = 1,2
    call line_gs_solver(mg_level(current_level)%pc,      &
                        mg_level(current_level)%app,     &
                        mg_level(current_level)%awp,     &
                        mg_level(current_level)%aep,     &
                        mg_level(current_level)%asp,     &
                        mg_level(current_level)%anp,     &
                        mg_level(current_level)%bp,     &
                        mg_level(current_level)%fp,     &
                        mg_level(current_level)%nx,     &
                        mg_level(current_level)%ny,     &
                        alphapc,                         &
                        mg_level(current_level)%xi,     &
                        mg_level(current_level)%xf      )
     end do


!if (1 .eq. 0) then
!    call tcoefficients( mg_level(current_level)%t,      &
!                        mg_level(current_level)%awt,      &
!                        mg_level(current_level)%aet,      &
!                        mg_level(current_level)%ast,      &
!                        mg_level(current_level)%ant,      &
!                        mg_level(current_level)%apt,      &
!                        mg_level(current_level)%bt,     &
!                        mg_level(current_level)%uf,      &
!                        mg_level(current_level)%vf,      &
!                        mg_level(current_level)%rho_uf,      &
!                        mg_level(current_level)%rho_vf,      &
!                        mg_level(current_level)%gammat_uf,      &
!                        mg_level(current_level)%gammat_vf,      &
!                        mg_level(current_level)%dx,      &
!                        mg_level(current_level)%dy,      &
!                        mg_level(current_level)%delx,      &
!                        mg_level(current_level)%dely,      &
!                        mg_level(current_level)%nx,      &
!                        mg_level(current_level)%ny       )
!
!    call line_gs_solver(mg_level(current_level)%t,      &
!                        mg_level(current_level)%apt,     &
!                        mg_level(current_level)%awt,     &
!                        mg_level(current_level)%aet,     &
!                        mg_level(current_level)%ast,     &
!                        mg_level(current_level)%ant,     &
!                        mg_level(current_level)%bt,     &
!                        mg_level(current_level)%ft,     &
!                        mg_level(current_level)%nx,     &
!                        mg_level(current_level)%ny,     &
!                        alphat,                         &
!                        mg_level(current_level)%dx,     &
!                        mg_level(current_level)%dy,     &
!                        mg_level(current_level)%delx,   &
!                        mg_level(current_level)%dely    )
!
!       call convergence(mg_level(current_level)%t,  &
!                    mg_level(current_level)%awt,   &
!                    mg_level(current_level)%aet,   &
!                    mg_level(current_level)%ast,   &
!                    mg_level(current_level)%ant,   &
!                    mg_level(current_level)%apt,   &
!                    mg_level(current_level)%bt,    &
!                    mg_level(current_level)%ft,    &
!                    mg_level(current_level)%dx,   &
!                    mg_level(current_level)%dy,   &
!                    mg_level(current_level)%nx,   &
!                    mg_level(current_level)%ny,   &
!                    residualt                      )
!end if
   ! mg_level(current_level)%pc = mg_level(current_level)%pc + mg_level(current_level)%pc

    call correction(mg_level(current_level)%u,      &
                    mg_level(current_level)%v,      &
                    mg_level(current_level)%uf,     &
                    mg_level(current_level)%vf,     &
                    mg_level(current_level)%pc,     &
                    mg_level(current_level)%p,      &
                    mg_level(current_level)%auf_bar,&
                    mg_level(current_level)%avf_bar,&
                    mg_level(current_level)%p_uf,   &
                    mg_level(current_level)%p_vf,   &
                    mg_level(current_level)%d_u,    &
                    mg_level(current_level)%d_v,    &
                    mg_level(current_level)%dx,     &
                    mg_level(current_level)%dy,     &
                    mg_level(current_level)%nx,     &
                    mg_level(current_level)%ny,     &
                    mg_level(current_level)%xi,     &
                    mg_level(current_level)%xf,     &
                    alphap                           )

   call convergence(mg_level(current_level)%u,  &
                    mg_level(current_level)%aw,   &
                    mg_level(current_level)%ae,   &
                    mg_level(current_level)%as,   &
                    mg_level(current_level)%an,   &
                    mg_level(current_level)%ap,   &
                    mg_level(current_level)%bu,    &
                    mg_level(current_level)%fu,    &
                    mg_level(current_level)%dx,   &
                    mg_level(current_level)%dy,   &
                    mg_level(current_level)%nx,   &
                    mg_level(current_level)%ny,   &
                    mg_level(current_level)%xi,     &
                    mg_level(current_level)%xf,     &
                    residualu                      )

   call convergence(mg_level(current_level)%v,  &
                    mg_level(current_level)%aw,   &
                    mg_level(current_level)%ae,   &
                    mg_level(current_level)%as,   &
                    mg_level(current_level)%an,   &
                    mg_level(current_level)%ap,   &
                    mg_level(current_level)%bv,    &
                    mg_level(current_level)%fv,    &
                    mg_level(current_level)%dx,   &
                    mg_level(current_level)%dy,   &
                    mg_level(current_level)%nx,   &
                    mg_level(current_level)%ny,   &
                    mg_level(current_level)%xi,     &
                    mg_level(current_level)%xf,     &
                    residualv                      )

   call convergencep(mg_level(current_level)%bp,    &
                     mg_level(current_level)%fp,    &
                     mg_level(current_level)%dx,   &
                     mg_level(current_level)%dy,   &
                     mg_level(current_level)%nx,   &
                     mg_level(current_level)%ny,   &
                     mg_level(current_level)%xi,     &
                     mg_level(current_level)%xf,     &
                     residualp                    )

!    if (taskid .eq. root)  then
!
!        write(*,*) current_level, i
!        call monitor_residuals
!
!    end if

 end do

CASE DEFAULT

   stop 'Something is wrong e.g. negative ''current_level'' '


END SELECT

end subroutine multigrid_vcycle
