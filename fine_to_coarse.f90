subroutine fine_to_coarse(r_h, nx, ny, xi_h, xf_h, xi_2h, xf_2h, r_2h)
use parallel_variables
use kind_parameters
implicit none

integer, intent(in):: nx, ny, xi_h, xf_h, xi_2h, xf_2h
real(kind=dp), intent(in), dimension(nx,xi_h-1:xf_h+1):: r_h
real(kind=dp), intent(out), dimension(nx/2+1,xi_2h-1:xf_2h+1):: r_2h
real(kind=dp):: w1, w2, w3, w4
integer:: imax, jmax
integer:: i, j

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

if (mod(nx,2) .eq. 0) then

	imax = nx/2+1

else

	imax = (nx+1)/2

end if

if (mod(ny,2) .eq. 0) then

	jmax = ny/2+1

else

	jmax = (ny+1)/2

end if

if (1 .eq. 1) then

do j = xi_2h,xf_2h

    do i = 2,imax-1

        r_2h(i,j) = 1.0d0/4.0d0*(r_h(2*i-2,2*j-2) + r_h(2*i-2,2*j-1) + r_h(2*i-1,2*j-2) + r_h(2*i-1,2*j-1))

     end do

end do

end if

if (1 .eq. 0) then

do j = xi_2h,xf_2h

    do i = 2,imax-1

        r_2h(i,j) = 1.0d0/4.0d0*3.0d0*(r_h(2*i-2,2*j-2) + r_h(2*i-2,2*j-1) + r_h(2*i-1,2*j-2) + r_h(2*i-1,2*j-1)) + 1.0d0/4.0d0*(r_2h(i-1,j) + r_2h(i+1,j) + r_2h(i,j-1) + r_2h(i,j+1))

     end do

end do

end if

if (xi_2h .eq. 2) then

w1 = 0.0d0
w2 = 2.0d0
w3 = 0.0d0
w4 = 2.0d0

do j = 1,1

 do i = 2,imax-1


    r_2h(i,j) = 1.0d0/4.0d0*(w2*r_h(2*i-2,2*j-1) + w4*r_h(2*i-1,2*j-1))

 end do

end do

else if (xf_2h .eq. jmax-1) then

w1 = 2.0d0
w2 = 0.0d0
w3 = 2.0d0
w4 = 0.0d0

do j = jmax,jmax

 do i = 2,imax-1


    r_2h(i,j) = 1.0d0/4.0d0*(w1*r_h(2*i-2,2*j-2) + w3*r_h(2*i-1,2*j-2))

 end do

end do

end if

w1 = 0.0d0
w2 = 0.0d0
w3 = 2.0d0
w4 = 2.0d0

do j = xi_2h,xf_2h

 do i = 1,1


    r_2h(i,j) = 1.0d0/4.0d0*(w3*r_h(2*i-1,2*j-2) + w4*r_h(2*i-1,2*j-1))

 end do

end do

w1 = 2.0d0
w2 = 2.0d0
w3 = 0.0d0
w4 = 0.0d0

do j = xi_2h,xf_2h

 do i = imax,imax

    r_2h(i,j) = 1.0d0/4.0d0*(w1*r_h(2*i-2,2*j-2) + w2*r_h(2*i-2,2*j-1))

 end do

end do

call exchange_data(r_2h, imax, xi_2h, xf_2h)

end subroutine
