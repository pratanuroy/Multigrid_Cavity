real(kind=dp) function apfun(Peclet)

use kind_parameters
implicit none

! Calculating Ap using Power law scheme

real(kind=dp):: Peclet

Peclet = 1.0d0-0.1d0*dabs(Peclet)
!
apfun = max(0.0d0, (Peclet*Peclet*Peclet*Peclet*Peclet))
!apfun = max(0.0d0, 1.0d0 -0.50d0*dabs(Peclet))
end function
