subroutine mg_correction(u_h, nx, ny, xi, xf, e_h, alpha)
use kind_parameters
implicit none

integer, intent(in):: nx, ny, xi, xf
real(kind=dp), intent(in), dimension(nx,xi-1:xf+1):: e_h
real(kind=dp), intent(inout), dimension(nx,xi-1:xf+1):: u_h
real(kind=dp), intent(in):: alpha
integer:: i,j
! Correct the values of u_h using e_h
! Input
! u_h   Variable in grid h
! e_h   Estimated error in grid h
! h     grid spacing
!
! Output
! u_h   Updated value of u_h

! Author: Pratanu Roy
! History:
! First Written: July 20, 2012

do j = xi,xf

    do i = 2,nx-1

        u_h(i,j) = u_h(i,j) + alpha*e_h(i,j)

    end do

end do

call exchange_data(u_h, nx, xi, xf)

end subroutine
