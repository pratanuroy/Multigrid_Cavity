subroutine check_convergence
!use global_variables
use parameters
use parallel_variables
use kind_parameters
implicit none
include "mpif.h"

        if (residualu .lt. epsilonu .and. residualv .lt. epsilonv .and. residualp .lt. epsilonp .and. residualt .lt. epsilont) then

            if (taskid .eq. root) then

                write(*,*) 'Solution converged in ', iteration, 'iterations!'
                write(*,'(I7, 6(E13.5))') iteration, residualu, residualv, residualp, residualt
                write(*,*) 'time step number = ', tcounter

            end if

            is_converged = 1

        end if

end subroutine
