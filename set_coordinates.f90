subroutine set_coordinates(delx, dely, xt, yt, nx, ny)
use kind_parameters
implicit none

integer:: nx, ny
real(kind=dp), intent(out), dimension(nx):: xt
real(kind=dp), intent(out), dimension(ny):: yt
real(kind=dp), dimension(nx-1):: delx
real(kind=dp), dimension(ny-1):: dely
integer:: it, jt

! locate the t position in x direction

xt(1) = 0.0

do it = 2,nx
    
    xt(it) = xt(it-1) + delx(it-1)
    
end do

! locate the t position in y direction

yt(1) = 0.0

do jt = 2,ny
    
    yt(jt) = yt(jt-1) + dely(jt-1)
    
end do

end subroutine set_coordinates
