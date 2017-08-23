subroutine residual(u_local, res, aP_u, aW_u, aE_u, aS_u, aN_u, f_u, dx_mg, dy_mg, nx, ny, alpha)
use kind_parameters
implicit none 

integer:: i, j, nx, ny
real(kind=dp), intent(in), dimension(nx,ny):: u_local, aP_u, f_u, aW_u, aE_u, aS_u, aN_u
real(kind=dp), intent(out), dimension(nx,ny):: res
real(kind=dp), dimension(nx):: dx_mg
real(kind=dp), dimension(ny):: dy_mg
real(kind=dp):: alpha

! Calculate the residual r = f - Au, where u_local is the guessed value
! Input
! u_local      Matrix for variables
! f      Right Hand Side
! hx     grid spacing in x direction
! hy     grid spacing in y direction
! 
! Output 
! r      Residual

! Author: Pratanu Roy
! History:
! First Written: July 20, 2012

!call geometry_multigrid(dx_mg, dy_mg, delx_mg, dely_mg, nx, ny)

alpha = 1.0

do j = 2,ny-1
    
    do i = 2,nx-1
        
        res(i,j) = f_u(i,j) + 1.0/(dx_mg(i)*dy_mg(j))*(aW_u(i,j)*u_local(i-1,j) + aE_u(i,j)*u_local(i+1,j) + aS_u(i,j)*u_local(i,j-1) + aN_u(i,j)*u_local(i,j+1) - aP_u(i,j)*u_local(i,j))
        
    end do
    
end do

end subroutine