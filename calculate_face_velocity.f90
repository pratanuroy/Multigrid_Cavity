subroutine calculate_face_velocity(u, v, uf, vf, nx, ny, xi, xf)
use kind_parameters
implicit none

integer, intent(in):: nx, ny, xi, xf
real(kind=dp), intent(in), dimension(nx,xi-1:xf+1):: u, v
real(kind=dp), intent(out), dimension(nx-1,xi-1:xf+1):: uf
real(kind=dp), intent(out), dimension(nx,xi-1:xf+1):: vf
real(kind=dp):: phi
integer:: i, j

do j = xi,xf

    do i = 2,nx-2

        phi = (u(i,j)-u(i-1,j))/(u(i+1,j) - u(i-1,j))
        uf(i,j) = 0.50d0*(u(i,j)+u(i+1,j))*phi + 0.50d0*(3.0d0*u(i,j)-u(i-1,j))*(1.0d0-phi)

    end do

end do

do j = xi,xf

    do i = 2,nx-1

        phi = (v(i,j)-v(i,j-1))/(v(i,j+1) - v(i,j-1))
        vf(i,j) = 0.50d0*(v(i,j)+v(i,j+1))*phi + 0.50d0*(3.0d0*v(i,j)-v(i,j-1))*(1.0d0-phi)

    end do

end do

if (xf .eq. ny-1) then

    do i = 2,nx-1

        vf(i,ny-1) = v(i,ny)

    end do

end if

call exchange_data(uf, nx-1, xi, xf)
call exchange_data(vf, nx, xi, xf)

end subroutine
