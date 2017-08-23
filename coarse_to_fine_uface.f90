subroutine coarse_to_fine_uface(e_2h, nx, ny, e_h)
use kind_parameters
implicit none

integer, intent(in):: nx, ny
real(kind=dp), intent(in), dimension(nx,ny):: e_2h
real(kind=dp), intent(out), dimension(2*nx-1,2*ny-2):: e_h
integer:: i, j

!write(*,*) nx, ny

! Apply linear interpolation to get values from coarse to fine grid
! Input
! e_2h   Variable in coarse grid
! hx, hy     grid spacing
!
! Output
! e_h     Variable in fine grid

! Author: Pratanu Roy
! History:
! First Written: July 20, 2012

do j = 2,ny-2

    do i = 1,nx-1

        e_h(2*i-1,2*j-1) = (1.0/4.0)*(3.0*e_2h(i,j) + e_2h(i,j+1))
        e_h(2*i,2*j) = (1.0/8.0)*(3.0*e_2h(i,j+1) + 3.0*e_2h(i+1,j+1) + e_2h(i,j) + e_2h(i+1,j))
        e_h(2*i, 2*j-1) = (1.0/8.0)*(e_2h(i,j+1) + e_2h(i+1,j+1) + 3.0*e_2h(i,j) + 3.0*e_2h(i+1,j))
        e_h(2*i-1,2*j) = (1.0/4.0)*(e_2h(i,j) + 3.0*e_2h(i,j+1))

    end do

end do

j = 1

do i = 1,nx-1

	e_h(2*i-1,2*j) = (1.0/2.0)*(e_2h(i,j) + e_2h(i,j+1))
    e_h(2*i,2*j) = (1.0/4.0)*(e_2h(i,j) + e_2h(i,j+1) + e_2h(i+1,j) + e_2h(i+1,j+1))

end do

j = ny-1

do i = 1,nx-1

    e_h(2*i-1,2*j-1) = (1.0/2.0)*(e_2h(i,j) + e_2h(i,j+1))
    e_h(2*i,2*j-1) = (1.0/4.0)*(e_2h(i,j) + e_2h(i,j+1) + e_2h(i+1,j) + e_2h(i+1,j+1))

end do

end subroutine
