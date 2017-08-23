subroutine mg_residual(phi, ap, aw, ae, as, an, b, f, nx, ny, xi, xf, res)
use kind_parameters
implicit none

integer, intent(in):: nx, ny, xi, xf
real(kind=dp), intent(in), dimension(nx,xi-1:xf+1):: phi, ap, b, aw, ae, as, an, f
real(kind=dp), intent(out), dimension(nx,xi-1:xf+1):: res
!real(kind=dp), dimension(nx), intent(in):: dx
!real(kind=dp), dimension(ny), intent(in):: dy
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

        !res(i,j) = 1.0/(dx(i)*dy(j))*(ap(i,j)*phi(i,j) - (aw(i,j)*phi(i-1,j) + ae(i,j)*phi(i+1,j) + as(i,j)*phi(i,j-1) + an(i,j)*phi(i,j+1)) - b(i,j)*dx(i)*dy(j) - f(i,j)*dx(i)*dy(j))
        res(i,j) = ap(i,j)*phi(i,j) - (aw(i,j)*phi(i-1,j) + ae(i,j)*phi(i+1,j) + as(i,j)*phi(i,j-1) + an(i,j)*phi(i,j+1)) - b(i,j)- f(i,j)
        !res(i,j) = 1.0/(dx(i)*dy(j))*(ap(i,j)*phi(i,j) - (aw(i,j)*phi(i-1,j) + ae(i,j)*phi(i+1,j) + as(i,j)*phi(i,j-1) + an(i,j)*phi(i,j+1)) - b(i,j)*dx(i)*dy(j) - f(i,j)*dx(i)*dy(j))

    end do

end do

call exchange_data(res, nx, xi, xf)

end subroutine
