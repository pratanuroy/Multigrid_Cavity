subroutine  reset_pressure(pc, imax, jmax, xi, xf)
!use variables
use kind_parameters
implicit none

integer, intent(in):: imax, jmax, xi, xf
real(kind=dp), dimension(imax,xi-1:xf+1), intent(inout):: pc

!! Local variables
integer:: ip, jp

do jp = xi-1,xf+1

    do ip = 1,imax

        pc(ip,jp) = 0.0d0

    end do

end do

end subroutine
