subroutine coordinates(xp, yp, dx, dy, imax, jmax)
use parameters
use kind_parameters
implicit none

integer, intent(in):: imax, jmax
real(kind=dp), intent(inout) :: xp(-1:imax+2), yp(-1:jmax+2)
real(kind=dp), intent(in) :: dx(imax), dy(jmax)
integer:: ip, jp


! locate the u position in x direction

!xu(1) = 0.0d0
!
!do iu = 2,iumax
!
!    xu(iu)=xu(iu-1) + delxu(iu-1)
!
!end do
!
!! locate the u position in y direction
!
!yu(1) = 0.0d0
!
!do ju = 2,jumax
!
!    yu(ju) = yu(ju-1) + delyu(ju-1)
!
!end do
!
!! locate the v position in x direction
!
!xv(1) = 0.0d0
!
!do iv = 2,ivmax
!
!  xv(iv) = xv(iv-1) + delxv(iv-1)
!
!end do
!
!! locate the v position in y direction
!
!do jv = 2,jvmax
!
!  yv(jv) = yv(jv-1) + delyv(jv-1)
!
!end do

! locate the p position in x direction

xp(1) = 0.0d0
xp(0) = xp(1) - 0.5*(dx(1) + dx(2))
xp(-1) = xp(0) - 0.5*(dx(2) + dx(3))

do ip = 2,imax

    xp(ip) = xp(ip-1) + 0.5*(dx(ip-1) + dx(ip))

end do

xp(imax+1) = xp(imax) + 0.5*(dx(imax-1) + dx(imax))
xp(imax+2) = xp(imax+1) + 0.5*(dx(imax-1) + dx(imax-2))

! locate the p position in y direction

yp(1) = 0.0d0
yp(0) = yp(1) - 0.5*(dy(1) + dy(2))
yp(-1) = yp(0) - 0.5*(dy(2) + dy(3))

do jp = 2,jmax

    yp(jp) = yp(jp-1) + 0.5*(dy(jp-1) + dy(jp))

end do

yp(imax+1) = yp(imax) + 0.5*(dy(jmax-1) + dy(jmax))
yp(imax+2) = yp(imax+1) + 0.5*(dy(jmax-1) + dy(jmax-2))

end subroutine
