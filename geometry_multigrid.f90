subroutine geometry_multigrid(dx_mg, dy_mg, delx_mg, dely_mg, imax, jmax)
use kind_parameters
implicit none

integer, intent(in):: imax, jmax
real(kind=dp), intent(inout), dimension(imax):: dx_mg
real(kind=dp), intent(inout), dimension(jmax)::  dy_mg
real(kind=dp), intent(inout), dimension(imax-1):: delx_mg
real(kind=dp), intent(inout), dimension(jmax-1):: dely_mg
real(kind=dp):: Lx, Ly, it, jt

Lx = 0.10d0
Ly = 0.10d0

! Calculate control volumes

dx_mg(1) = 0.0

do it=2,imax-1

 dx_mg(it) = Lx/(imax-2)

end do

dx_mg(imax) = 0.0

dy_mg(1) = 0.0

do jt = 2,jmax-1

 dy_mg(jt) = Ly/(jmax-2)

end do

dy_mg(jmax) = 0.0

! Calculate diffusion lengths

delx_mg(1) = 0.5*dx_mg(2)

do it = 2,imax-2

 delx_mg(it) = 0.5*(dx_mg(it)+dx_mg(it+1))

end do

delx_mg(imax-1) = 0.5*dx_mg(imax-1)

dely_mg(1) = 0.5*dy_mg(2)

do jt = 2,jmax-2

 dely_mg(jt) = 0.5*(dy_mg(jt) + dy_mg(jt+1))

end do

dely_mg(jmax-1) = 0.5*dy_mg(jmax-1)

end subroutine
