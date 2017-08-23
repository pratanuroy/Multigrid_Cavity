subroutine mg_error(phi_prev, phi, nx, ny, xi, xf, error_phi)
use kind_parameters
implicit none

integer, intent(in):: nx, ny, xi, xf
real(kind=dp), intent(in), dimension(nx,xi-1:xf+1):: phi_prev, phi
real(kind=dp), intent(out), dimension(nx,xi-1:xf+1):: error_phi
integer:: i, j

! Calculate the residual r = f - Au, where u is the guessed value
! Input
! u      Matrix for variables
! f      Right Hand Side
! hx     grid spacing in x direction
! hy     grid spacing in y direction
!
! Output
! r      Residual

! Author: Pratanu Roy
! History:
! First Written: July 20, 2012
! Modified: Feb 11, 2013

do j = xi,xf

    do i = 2,nx-1

        error_phi(i,j) = phi(i,j) - phi_prev(i,j)

    end do

end do

call exchange_data(error_phi, nx, xi, xf)

end subroutine
