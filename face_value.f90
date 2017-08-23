subroutine face_value(phi, phi_uf, phi_vf, dx, dy, nx, ny, xi, xf)
use kind_parameters
implicit none

integer, intent(in):: nx, ny, xi, xf
real(kind=dp), intent(inout), dimension(nx,xi-1:xf+1):: phi
real(kind=dp), intent(out), dimension(nx-1,xi-1:xf+1):: phi_uf
real(kind=dp), intent(out), dimension(nx,xi-1:xf+1):: phi_vf
real(kind=dp), dimension(nx), intent(in):: dx
real(kind=dp), dimension(ny), intent(in):: dy
integer:: i, j
real(kind=dp) :: fe_area, fn_area
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

! Calculate boundary values by linear extrapolation
!
do j = xi-1,xf+1

    phi(1,j) = (15.0d0*phi(2,j) - 10.0d0*phi(3,j) + 3.0d0*phi(4,j))/8.0d0
    phi(nx,j) = (15.0d0*phi(nx-1,j) - 10.0d0*phi(nx-2,j) + 3.0d0*phi(nx-3,j))/8.0d0

end do

if (xi .eq. 2) then

do i = 2,nx-1

      phi(i,1) = (15.0d0*phi(i,2) - 10.0d0*phi(i,3) + 3.0d0*phi(i,4))/8.0d0
!      phi(i,1) = phi(i,2)

end do

else if (xf .eq. ny-1) then

do i = 2,nx-1

      phi(i,ny) = (15.0d0*phi(i,ny-1) - 10.0d0*phi(i,ny-2) + 3.0d0*phi(i,ny-3))/8.0d0
!     phi(i,ny) = phi(i,ny-1)

end do

end if

!phi(1,:) = phi(2,:)
!phi(nx,:) = phi(nx-1,:)
!phi(:,1) = phi(:,2)
!phi(:,ny) = phi(:,ny-1)

call exchange_data(phi, nx, xi, xf)

do j = xi,xf

     do i = 2,nx-1

            fe_area = dx(i+1)/(dx(i+1)+dx(i))
            phi_uf(i,j) = fe_area*phi(i,j) + (1.0d0-fe_area)*phi(i+1,j)

      end do

end do

phi_uf(1,:) = phi(1,:)
!phi_uf(1,:) = phi_uf(2,:)
!phi_uf(nx-1,:) = phi_uf(nx-2,:)

do j = xi,xf

    do i = 2,nx-1

            fn_area = dy(j+1)/(dy(j+1)+dy(j))
            phi_vf(i,j) = fn_area*phi(i,j) + (1.0d0-fn_area)*phi(i,j+1)

    end do

end do

if (xi .eq. 2) then

phi_vf(1:nx,1) = phi(1:nx,1)
!phi_vf(1:nx,1) = phi_vf(1:nx,2)

end if

if (xf .eq. ny-1) then

phi_vf(1:nx,ny-1) = phi(1:nx,ny)
!phi_vf(1:nx,ny-1) = phi_vf(1:nx,ny-2)

end if

call exchange_data(phi_uf, nx-1, xi, xf)
call exchange_data(phi_vf, nx, xi, xf)

end subroutine
