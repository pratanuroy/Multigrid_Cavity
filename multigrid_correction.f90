subroutine multigrid_correction(u_h, e_h, nx, ny)
!use variables
use kind_parameters
implicit none

integer:: i, j, nx, ny
real(kind=dp), intent(in), dimension(nx,ny):: e_h
real(kind=dp), intent(inout), dimension(nx,ny):: u_h
real(kind=dp):: alpha

! Correct the values of u_h using e_h
! Input
! u_h   Variable in grid h
! e_h   Estimated error in grid h
! h     grid spacing
! 
! Output 
! u_h   Updated value of u_h

! Author: Pratanu Roy
! History:
! First Written: July 20, 2012

do j = 2,ny-1
    
    do i = 2,nx-1
   
        u_h(i,j) = u_h(i,j) + e_h(i,j)
        
    end do
    
end do

end subroutine