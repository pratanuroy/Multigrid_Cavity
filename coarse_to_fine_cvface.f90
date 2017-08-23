subroutine coarse_to_fine_cvface(e_2h, e_h, nx, ny)
!use variables
use kind_parameters
implicit none 

integer:: i, j, if, jf, ic, jc, nx, ny
real(kind=dp), intent(in), dimension(nx,ny):: e_2h
real(kind=dp), intent(out), dimension(:,:):: e_h

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
    
    do i = 2,nx-2
   
        e_h(2*i-1,2*j-1) = (1.0/4.0)*(3.0*e_2h(i,j) + e_2h(i,j+1))
        e_h(2*i-1,2*j-1) = (1.0/16.0)*(9.0*e_2h(i+1,j) + 3.0*e_2h(i,j) + 3.0*e_2h(i+1,j+1) + e_2h(i,j+1))
        e_h(2*i-1,2*j) = (1.0/16.0)*(9.0*e_2h(i,j+1) + 3.0*e_2h(i,j) + 3.0*e_2h(i+1,j+1) + e_2h(i+1,j))
        e_h(2*i-1,2*j-1) = (1.0/16.0)*(9.0*e_2h(i,j) + 3.0*e_2h(i+1,j) + 3.0*e_2h(i,j+1) + e_2h(i+1,j+1))
            
    end do

end do

j = 1

    e_h(2,2*j) = (1.0/4.0)*(e_2h(1,j) + e_2h(1,j+1) + e_2h(2,j) + e_2h(2,j+1))
    
    do i = 2,nx-2
   
        
        e_h(2*i-1,2*j) = (1.0/8.0)*(3.0*e_2h(i,j) + 3.0*e_2h(i,j+1) + e_2h(i+1,j) + e_2h(i+1,j+1))
        e_h(2*i,2*j) = (1.0/8.0)*(e_2h(i,j) + e_2h(i,j+1) + 3.0*e_2h(i+1,j) + 3.0*e_2h(i+1,j+1))
            
    end do

    e_h(2*nx-3,2*j) = (1.0/4.0)*(e_2h(nx-1,j) + e_2h(nx-1,j+1) + e_2h(nx,j) + e_2h(nx,j+1))


do j = 2,ny-2

    i = 1 
        
        e_h(2*i,2*j-1) = (1.0/8.0)*(3.0*e_2h(i,j) + 3.0*e_2h(i+1,j) + e_2h(i,j+1) + e_2h(i+1,j+1))
        e_h(2*i,2*j) = (1.0/8.0)*(e_2h(i,j) + e_2h(i+1,j) + 3.0*e_2h(i,j+1) + 3.0*e_2h(i+1,j+1))

    i = nx-1 
        
        e_h(2*i-1,2*j-1) = (1.0/8.0)*(3.0*e_2h(i,j) + 3.0*e_2h(i+1,j) + e_2h(i,j+1) + e_2h(i+1,j+1))
        e_h(2*i-1,2*j) = (1.0/8.0)*(e_2h(i,j) + e_2h(i+1,j) + 3.0*e_2h(i,j+1) + 3.0*e_2h(i+1,j+1))

end do

j = ny-1
    
    e_h(2,2*j-1) = (1.0/4.0)*(e_2h(1,j) + e_2h(1,j+1) + e_2h(2,j) + e_2h(2,j+1))
    
    do i = 2,nx-2
        
        e_h(2*i-1,2*j-1) = (1.0/8.0)*(3.0*e_2h(i,j) + 3.0*e_2h(i,j+1) + e_2h(i+1,j) + e_2h(i+1,j+1))
        e_h(2*i,2*j-1) = (1.0/8.0)*(e_2h(i,j) + e_2h(i,j+1) + 3.0*e_2h(i+1,j) + 3.0*e_2h(i+1,j+1))
            
    end do
    
    e_h(2*nx-3,2*j-1) = (1.0/4.0)*(e_2h(nx-1,j) + e_2h(nx-1,j+1) + e_2h(nx,j) + e_2h(nx,j+1))  

end subroutine