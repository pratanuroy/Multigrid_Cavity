subroutine fine_to_coarse_res(r_h, nx, ny, xi_h, xf_h, xi_2h, xf_2h, r_2h)
use kind_parameters
implicit none

integer, intent(in):: nx, ny, xi_h, xf_h, xi_2h, xf_2h
real(kind=dp), intent(in), dimension(nx,xi_h-1:xf_h+1):: r_h
real(kind=dp), intent(out), dimension(nx/2+1,xi_2h-1:xf_2h+1):: r_2h
integer:: i,j, imax

imax = nx/2+1

! Apply full weighting to get interpolated values from fine to coarse grid
! Input
! r_h   Variable in fine grid
! f     Right Hand Side
! hx, hy     grid spacing
!
! Output
! r_2h     Variable in coarse grid

! Author: Pratanu Roy
! History:
! First Written: July 20, 2012
! Indices Rechecked on Jan 13, 2013

do j = xi_2h,xf_2h

    do i = 2,nx/2

        r_2h(i,j) = (r_h(2*i-2,2*j-2) + r_h(2*i-2,2*j-1) + r_h(2*i-1,2*j-2) + r_h(2*i-1,2*j-1))

     end do

end do

call exchange_data(r_2h, imax, xi_2h, xf_2h)

end subroutine
