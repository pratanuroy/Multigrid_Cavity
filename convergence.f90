subroutine  convergence(phi, aw, ae, as, an, ap, b_phi, f_phi, dx, dy, imax, jmax, xi, xf, residual)
use parallel_variables
use parameters
use kind_parameters
implicit none
include "mpif.h"

integer:: imax, jmax, iu, ju, xi, xf
real(kind=dp), dimension(imax, xi-1:xf+1):: phi, aw, ae, as, an, ap, b_phi, f_phi
real(kind=dp):: dx(imax), dy(jmax)
real(kind=dp):: residual
real(kind=dp):: local_sum_num, local_sum_den, temp1, temp2, global_sum_num, global_sum_den

local_sum_num = 0.0d0
local_sum_den = 0.0d0

do ju = xi,xf

      do iu = 2,imax-1

         temp1 = aw(iu,ju)*phi(iu-1,ju) + ae(iu,ju)*phi(iu+1,ju) + as(iu,ju)*phi(iu,ju-1) + an(iu,ju)*phi(iu,ju+1)
         temp2 = ap(iu,ju)*phi(iu,ju)
         !sum_num = sum_num + dabs(temp2-temp1-b_phi(iu,ju)*dx(iu)*dy(ju)-f_phi(iu,ju)*dx(iu)*dy(ju))
         local_sum_num = local_sum_num + dabs(temp2-temp1-b_phi(iu,ju)-f_phi(iu,ju))
         local_sum_den = local_sum_den + dabs(temp2)

      end do

end do

tlocal1 = MPI_WTime()

call MPI_Allreduce(local_sum_num, global_sum_num, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
call MPI_Allreduce(local_sum_den, global_sum_den, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

tlocal2 = MPI_WTime()

communication_time = communication_time + (tlocal2 - tlocal1)

residual = global_sum_num/global_sum_den

end subroutine
