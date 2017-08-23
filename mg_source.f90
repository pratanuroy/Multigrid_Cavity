subroutine mg_source(phi, ap, aw, ae, as, an, b, res, nx, ny, dx, dy, source)
use kind_parameters
implicit none

integer:: i, j, nx, ny
real(kind=dp), intent(in), dimension(nx,ny):: phi, ap, b, res, aw, ae, as, an
real(kind=dp), intent(out), dimension(nx,ny):: source
real(kind=dp), dimension(nx), intent(in):: dx
real(kind=dp), dimension(ny), intent(in):: dy

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

do j = 2,ny-1

    do i = 2,nx-1

        source(i,j) = (ap(i,j)*phi(i,j) - (aw(i,j)*phi(i-1,j) + ae(i,j)*phi(i+1,j) + as(i,j)*phi(i,j-1) + an(i,j)*phi(i,j+1)) - b(i,j) - res(i,j))

    end do

end do


end subroutine
