subroutine monitor_residuals
use global_variables
use parameters
use parallel_variables
use kind_parameters
implicit none

if (iteration .eq. 1 .and. taskid .eq. root) then

   write(*,*) 'iteration', ' ', 'u-residual', ' ', 'v-residual', ' ', 'continuity-residual'

end if

if (mod(iteration,monitor_step) .eq. 0 .and. taskid .eq. root) then

   !write(*,'(I7, 4(E13.5))') iteration, residualu, residualv, residualp, residualt
   write(*,'(I7, 3(E13.5))') iteration, residualu, residualv, residualp

end if

end subroutine
