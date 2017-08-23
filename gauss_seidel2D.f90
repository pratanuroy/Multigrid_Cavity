subroutine gauss_seidel2D(phi, ap, aw, ae, as, an, f, nx, ny, dx, dy, alpha)
use kind_parameters
implicit none

integer:: i, j, nx, ny
real(kind=dp), intent(inout), dimension(nx,ny):: phi
real(kind=dp), intent(in), dimension(nx,ny):: ap, aw, ae, as, an, f
real(kind=dp), dimension(nx):: dx
real(kind=dp), dimension(ny):: dy
real(kind=dp) :: alpha

! Apply Gauss Seidel algorithm to get an updated value of the variable u
! Input
! t      Matrix for variables
! f      Right Hand Side
! Output 
! t      Update variable

! Author: Pratanu Roy
! History:
! First Written: July 19, 2012

! Point by point method / Gauss-Siedel iteration method

do j = 2,ny-1

 do i = 2,nx-1
 
    phi(i,j) = phi(i,j) + alpha/ap(i,j)*(f(i,j)*dx(i)*dy(j) + aw(i,j)*phi(i-1,j) + ae(i,j)*phi(i+1,j) + as(i,j)*phi(i,j-1) + an(i,j)*phi(i,j+1) - ap(i,j)*phi(i,j))

 end do

end do

end subroutine