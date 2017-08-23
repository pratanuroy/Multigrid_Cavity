subroutine  convergencep(bp, fp, dx, dy, imax, jmax, xi, xf, residual)
use parameters
use kind_parameters
implicit none
include "mpif.h"

integer, intent(in):: imax, jmax, xi, xf
real(kind=dp), intent(in), dimension(imax, xi-1:xf+1):: bp, fp
real(kind=dp), intent(in):: dx(imax), dy(jmax)
real(kind=dp), intent(out):: residual

real(kind=dp):: sum_res, global_sum
integer:: ip, jp
integer:: ierr

sum_res = 0.0d0
global_sum = 0.0d0
do jp = xi,xf

    do ip = 2,imax-1

        sum_res = sum_res + dabs(bp(ip,jp)+fp(ip,jp))

    end do

end do


call MPI_Allreduce(sum_res, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

residual = global_sum/(v_ref*hydraulic_dia)

end subroutine
