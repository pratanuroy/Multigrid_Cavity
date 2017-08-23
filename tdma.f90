subroutine tdma(a, b, N)
use kind_parameters
implicit none

integer*4:: i, N
real(kind=dp), dimension(N,3):: a
real(kind=dp), intent(inout), dimension(N):: b 

a(1,3) = -a(1,3)/a(1,2)
b(1) = b(1)/a(1,2)

do i = 2,N
    
    a(i,3) = -a(i,3)/(a(i,2) + a(i,1)*a(i-1,3))
    b(i) = (b(i) - a(i,1)*b(i-1))/(a(i,2) + a(i,1)*a(i-1,3))
    
end do

do i = N-1,1,-1
    
    b(i) = a(i,3)*b(i+1) + b(i)
    
end do

end subroutine
