recursive subroutine multigrid_vcycle(current_level)
USE Multigrid_Levels
!USE Global_Variables
USE Parallel_Variables
use kind_parameters
implicit none

integer:: i, current_level
real(kind=dp), dimension(mg_level(current_level)%nx, mg_level(current_level)%xi-1:mg_level(current_level)%xf+1) :: res, error_h
real(kind=dp):: residualt
real(kind=dp):: alpha

alpha = 1.0d0

!write(*,*) 'Current Level = ', current_level

SELECT CASE(current_level)

CASE(1)

 do i = 1,5

if (1 .eq. 0) then
  call tcoefficients(mg_level(current_level)%phi,   &
                     mg_level(current_level)%ap,   &
                     mg_level(current_level)%aw,   &
                     mg_level(current_level)%ae,   &
                     mg_level(current_level)%as,   &
                     mg_level(current_level)%an,   &
                     mg_level(current_level)%b,    &
                     mg_level(current_level)%nx,   &
                     mg_level(current_level)%ny,   &
                     mg_level(current_level)%dx,   &
                     mg_level(current_level)%dy,   &
                     mg_level(current_level)%delx, &
                     mg_level(current_level)%dely, &
                     mg_level(current_level)%xi,   &
                     mg_level(current_level)%xf    )
end do
  call line_gs_solver(mg_level(current_level)%pc, &
                      mg_level(current_level)%ap,  &
                      mg_level(current_level)%aw,  &
                      mg_level(current_level)%ae,  &
                      mg_level(current_level)%as,  &
                      mg_level(current_level)%an,  &
                      mg_level(current_level)%bp,   &
                      mg_level(current_level)%nx,  &
                      mg_level(current_level)%ny,  &
                      mg_level(current_level)%dx,  &
                      mg_level(current_level)%dy,  &
                      mg_level(current_level)%xi,  &
                      mg_level(current_level)%xf,  &
                      alpha                 )

   call calc_residualt(mg_level(current_level)%pc,  &
                       mg_level(current_level)%ap,   &
                       mg_level(current_level)%aw,   &
                       mg_level(current_level)%ae,   &
                       mg_level(current_level)%as,   &
                       mg_level(current_level)%an,   &
                       mg_level(current_level)%bp,    &
                       mg_level(current_level)%nx,   &
                       mg_level(current_level)%ny,   &
                       mg_level(current_level)%dx,   &
                       mg_level(current_level)%dy,   &
                       mg_level(current_level)%xi,   &
                       mg_level(current_level)%xf,   &
                       residualt                      )

    !if (taskid .eq. root) write(*,*) current_level, iteration, residualt
    !if (taskid .eq. root+1) write(*,*) current_level, i, residualt

 end do

CASE(2:)

 do i = 1,3
if (1 .eq. 0) then
  call tcoefficients(mg_level(current_level)%phi,   &
                     mg_level(current_level)%ap,   &
                      mg_level(current_level)%aw,   &
                      mg_level(current_level)%ae,   &
                      mg_level(current_level)%as,   &
                      mg_level(current_level)%an,   &
                      mg_level(current_level)%b,    &
                      mg_level(current_level)%nx,   &
                      mg_level(current_level)%ny,   &
                      mg_level(current_level)%dx,   &
                      mg_level(current_level)%dy,   &
                      mg_level(current_level)%delx, &
                      mg_level(current_level)%dely, &
                      mg_level(current_level)%xi,   &
                      mg_level(current_level)%xf    )
end if
  call line_gs_solver(mg_level(current_level)%phi, &
                      mg_level(current_level)%ap,  &
                      mg_level(current_level)%aw,  &
                      mg_level(current_level)%ae,  &
                      mg_level(current_level)%as,  &
                      mg_level(current_level)%an,  &
                      mg_level(current_level)%b,   &
                      mg_level(current_level)%nx,  &
                      mg_level(current_level)%ny,  &
                      mg_level(current_level)%dx,  &
                      mg_level(current_level)%dy,  &
                      mg_level(current_level)%xi,  &
                      mg_level(current_level)%xf,  &
                      alpha                 )

   call calc_residualt(mg_level(current_level)%phi,  &
                       mg_level(current_level)%ap,   &
                       mg_level(current_level)%aw,   &
                       mg_level(current_level)%ae,   &
                       mg_level(current_level)%as,   &
                       mg_level(current_level)%an,   &
                       mg_level(current_level)%b,    &
                       mg_level(current_level)%nx,   &
                       mg_level(current_level)%ny,   &
                       mg_level(current_level)%dx,   &
                       mg_level(current_level)%dy,   &
                       mg_level(current_level)%xi,   &
                       mg_level(current_level)%xf,   &
                       residualt                      )

    !if (taskid .eq. root) write(*,*) current_level, iteration, residualt
    !if (taskid .eq. root+1) write(*,*) current_level, i, residualt

 end do

 call calculate_residual(mg_level(current_level)%phi, &
                         mg_level(current_level)%ap,  &
                         mg_level(current_level)%aw,  &
                         mg_level(current_level)%ae,  &
                         mg_level(current_level)%as,  &
                         mg_level(current_level)%an,  &
                         mg_level(current_level)%b,   &
                         mg_level(current_level)%nx,  &
                         mg_level(current_level)%ny,  &
                         mg_level(current_level)%dx,  &
                         mg_level(current_level)%dy,  &
                         mg_level(current_level)%xi,  &
                         mg_level(current_level)%xf,  &
                         res                          )

!  write(*,*) 'Norm of res = ', sum(abs(res))

  call fine_to_coarse(res, mg_level(current_level)%nx,  &
                      mg_level(current_level)%ny,       &
                      mg_level(current_level)%xi,      &
                      mg_level(current_level)%xf,      &
                      mg_level(current_level-1)%xi,    &
                      mg_level(current_level-1)%xf,    &
                      mg_level(current_level-1)%b      )

  mg_level(current_level-1)%phi(:,mg_level(current_level-1)%xi-1:mg_level(current_level-1)%xf+1) = 0.0d0

  call multigrid_vcycle(current_level-1)

  call coarse_to_fine(mg_level(current_level-1)%phi, &
                    mg_level(current_level-1)%nx,    &
                    mg_level(current_level-1)%ny,    &
                    mg_level(current_level-1)%xi,    &
                    mg_level(current_level-1)%xf,    &
                    mg_level(current_level)%xi,    &
                    mg_level(current_level)%xf,    &
                    error_h                          )
!write(*,*) error
  call calculate_correction(mg_level(current_level)%phi, &
                            mg_level(current_level)%nx,  &
                            mg_level(current_level)%ny,  &
                            mg_level(current_level)%xi,  &
                            mg_level(current_level)%xf,  &
                            error_h                        )

  do i = 1,3
if (1 .eq. 0) then
  call tcoefficients(mg_level(current_level)%phi,   &
                     mg_level(current_level)%ap,   &
                      mg_level(current_level)%aw,   &
                      mg_level(current_level)%ae,   &
                      mg_level(current_level)%as,   &
                      mg_level(current_level)%an,   &
                      mg_level(current_level)%b,    &
                      mg_level(current_level)%nx,   &
                      mg_level(current_level)%ny,   &
                      mg_level(current_level)%dx,   &
                      mg_level(current_level)%dy,   &
                      mg_level(current_level)%delx, &
                      mg_level(current_level)%dely, &
                      mg_level(current_level)%xi,   &
                      mg_level(current_level)%xf    )
end if
  call line_gs_solver(mg_level(current_level)%pc, &
                      mg_level(current_level)%ap,  &
                      mg_level(current_level)%aw,  &
                      mg_level(current_level)%ae,  &
                      mg_level(current_level)%as,  &
                      mg_level(current_level)%an,  &
                      mg_level(current_level)%bp,   &
                      mg_level(current_level)%nx,  &
                      mg_level(current_level)%ny,  &
                      mg_level(current_level)%dx,  &
                      mg_level(current_level)%dy,  &
                      mg_level(current_level)%xi,  &
                      mg_level(current_level)%xf,  &
                      alpha                 )

   call calc_residualt(mg_level(current_level)%pc,  &
                       mg_level(current_level)%ap,   &
                       mg_level(current_level)%aw,   &
                       mg_level(current_level)%ae,   &
                       mg_level(current_level)%as,   &
                       mg_level(current_level)%an,   &
                       mg_level(current_level)%bp,    &
                       mg_level(current_level)%nx,   &
                       mg_level(current_level)%ny,   &
                       mg_level(current_level)%dx,   &
                       mg_level(current_level)%dy,   &
                       mg_level(current_level)%xi,   &
                       mg_level(current_level)%xf,   &
                       residualt                      )

    !if (taskid .eq. root) write(*,*) current_level, iteration, residualt
    !if (taskid .eq. root+1) write(*,*) current_level, i, residualt

 end do

END SELECT

end subroutine multigrid_vcycle
