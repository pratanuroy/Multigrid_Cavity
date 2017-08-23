 subroutine steady_solver_mg(ilevel)
USE global_Variables
 USE parameters
 USE multigrid_levels
 USE parallel_variables
use kind_parameters
 implicit none
  include "mpif.h"

 !integer::  ip, jp
 integer, intent(inout) :: ilevel
!if (mg_handle .eq. 'on') then

 call multigrid_vcycle(ilevel)
 
 if (residualp .lt. epsilonp .and. residualu .lt. epsilonu .and. residualv .lt. epsilonv .and. ilevel .lt. num_level) then

    residualu = 1.0d0
    residualv = 1.0d0
    residualp = 1.0d0
    residualt = 0.0d0

    ilevel = ilevel + 1

    if (taskid .eq. root+1) write(*,*) 'Going to next level = ', ilevel

    call coarse_to_fine(mg_level(ilevel-1)%u,  &
                        mg_level(ilevel-1)%nx,       &
                        mg_level(ilevel-1)%ny,       &
                        mg_level(ilevel-1)%xi,         &
                        mg_level(ilevel-1)%xf,         &
                        mg_level(ilevel)%xi,         &
                        mg_level(ilevel)%xf,         &
                        mg_level(ilevel)%u     )

    call coarse_to_fine(  mg_level(ilevel-1)%v,  &
                         mg_level(ilevel-1)%nx,       &
                         mg_level(ilevel-1)%ny,       &
                         mg_level(ilevel-1)%xi,         &
                        mg_level(ilevel-1)%xf,         &
                        mg_level(ilevel)%xi,         &
                        mg_level(ilevel)%xf,         &
                        mg_level(ilevel)%v     )

    call coarse_to_fine(  mg_level(ilevel-1)%p,  &
                        mg_level(ilevel-1)%nx,       &
                        mg_level(ilevel-1)%ny,       &
                        mg_level(ilevel-1)%xi,         &
                        mg_level(ilevel-1)%xf,         &
                        mg_level(ilevel)%xi,         &
                        mg_level(ilevel)%xf,         &
                        mg_level(ilevel)%p     )

    call face_value(   mg_level(ilevel)%p,      &
                        mg_level(ilevel)%p_uf,     &
                        mg_level(ilevel)%p_vf,     &
                        mg_level(ilevel)%dx,     &
                        mg_level(ilevel)%dy,     &
                        mg_level(ilevel)%nx,     &
                        mg_level(ilevel)%ny,     &
                        mg_level(ilevel)%xi,    &
                       mg_level(ilevel)%xf      )

!     call calculate_face_velocity(mg_level(ilevel)%u,      &
!                            mg_level(ilevel)%v,     &
!                            mg_level(ilevel)%uf,     &
!                            mg_level(ilevel)%vf,     &
!                            mg_level(ilevel)%nx,     &
!                            mg_level(ilevel)%ny,     &
!                            mg_level(ilevel)%xi,    &
!                            mg_level(ilevel)%xf     )
!
    call coarse_to_fine(mg_level(ilevel-1)%ap,  &
                        mg_level(ilevel-1)%nx,       &
                        mg_level(ilevel-1)%ny,       &
                        mg_level(ilevel-1)%xi,         &
                        mg_level(ilevel-1)%xf,         &
                        mg_level(ilevel)%xi,         &
                        mg_level(ilevel)%xf,         &
                        mg_level(ilevel)%ap     )
if (1 .eq. 0) then
    call coarse_to_fine(mg_level(ilevel-1)%bp,  &
                        mg_level(ilevel-1)%nx,       &
                        mg_level(ilevel-1)%ny,       &
                        mg_level(ilevel-1)%xi,         &
                        mg_level(ilevel-1)%xf,         &
                        mg_level(ilevel)%xi,         &
                        mg_level(ilevel)%xf,         &
                        mg_level(ilevel)%ap     )
end if
    call source_terms(         mg_level(ilevel)%u0,      &
                               mg_level(ilevel)%v0,      &
                               mg_level(ilevel)%uf0,      &
                               mg_level(ilevel)%vf0,      &
                               mg_level(ilevel)%ap_zero, &
                               mg_level(ilevel)%f_uvel,     &
                               mg_level(ilevel)%f_vvel,     &
                               mg_level(ilevel)%su,      &
                               mg_level(ilevel)%sv,      &
                               mg_level(ilevel)%bu,      &
                               mg_level(ilevel)%bv,  &
                               mg_level(ilevel)%s_uf,  &
                               mg_level(ilevel)%s_vf,   &
                               mg_level(ilevel)%rho_vf,   &
                               mg_level(ilevel)%p_uf,  &
                               mg_level(ilevel)%p_vf,   &
                               mg_level(ilevel)%dx,      &
                               mg_level(ilevel)%dy,      &
                               mg_level(ilevel)%nx,      &
                               mg_level(ilevel)%ny,      &
                               mg_level(ilevel)%xi,      &
                               mg_level(ilevel)%xf       )

     call calculate_face_velocity_interp(mg_level(ilevel)%u,      &
                            mg_level(ilevel)%v,     &
                            mg_level(ilevel)%uf,     &
                            mg_level(ilevel)%vf,     &
                            mg_level(ilevel)%p,     &
                            mg_level(ilevel)%p_uf,     &
                            mg_level(ilevel)%p_vf,     &
                            mg_level(ilevel)%s_uf,     &
                            mg_level(ilevel)%s_vf,     &
                            mg_level(ilevel)%ap,     &
                            mg_level(ilevel)%dx,     &
                            mg_level(ilevel)%dy,     &
                            mg_level(ilevel)%nx,     &
                            mg_level(ilevel)%ny,     &
                            mg_level(ilevel)%xi,    &
                            mg_level(ilevel)%xf     )
if (1 .eq. 0) then
    call correction(mg_level(ilevel)%u,      &
                    mg_level(ilevel)%v,      &
                    mg_level(ilevel)%uf,     &
                    mg_level(ilevel)%vf,     &
                    mg_level(ilevel)%pc,     &
                    mg_level(ilevel)%p,      &
                    mg_level(ilevel)%auf_bar,&
                    mg_level(ilevel)%avf_bar,&
                    mg_level(ilevel)%p_uf,   &
                    mg_level(ilevel)%p_vf,   &
                    mg_level(ilevel)%d_u,    &
                    mg_level(ilevel)%d_v,    &
                    mg_level(ilevel)%dx,     &
                    mg_level(ilevel)%dy,     &
                    mg_level(ilevel)%nx,     &
                    mg_level(ilevel)%ny,     &
                    mg_level(ilevel)%xi,     &
                    mg_level(ilevel)%xf,     &
                    alphap                           )
end if
!      call coarse_to_fine_ls(mg_level(ilevel-1)%T, &
!                    mg_level(ilevel-1)%nx,    &
!                    mg_level(ilevel-1)%ny,    &
!                    mg_level(ilevel-1)%xi,    &
!                    mg_level(ilevel-1)%xf,    &
!                    mg_level(ilevel)%xi,    &
!                    mg_level(ilevel)%xf,    &
!                    mg_level(ilevel)%T,    &
!                    mg_level(ilevel)%delx,    &
!                    mg_level(ilevel)%dely  )

 end if

end subroutine
