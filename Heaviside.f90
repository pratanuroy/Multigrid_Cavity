subroutine Heaviside(phi, result, eps_heaviside)
use parameters
use kind_parameters
implicit none

real(kind=dp), intent(in):: phi, eps_heaviside
real(kind=dp), intent(out):: result
!real(kind=dp), parameter:: pi = 3.1415926535897932384626433832795028841971693993751058209d0

if (phi .lt. -1.0d0*eps_heaviside) then

 result = 0.0d0

else if (phi .gt. eps_heaviside) then

 result = 1.0d0

else

 result = (phi + eps_heaviside)/(2.0d0*eps_heaviside) + (0.5d0/pi)*sin(pi*phi/eps_heaviside)

end if

end subroutine
