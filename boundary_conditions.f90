subroutine boundary_conditions(u, v, uf, vf, imax, jmax, xi, xf)
use parameters
use parallel_variables
use kind_parameters
implicit none
include "mpif.h"

integer, intent(in):: imax, jmax, xi, xf
real(kind=dp), dimension(imax, xi-1:xf+1):: u, v
real(kind=dp), dimension(imax-1, xi-1:xf+1):: uf
real(kind=dp), dimension(imax, xi-1:xf+1):: vf
integer:: iu, ju, iv, jv
! Boundary Conditions for u-velocity

do ju = xi-1,xf+1

 u(1,ju) = 0.0d0
 u(imax,ju) = 0.0d0

end do

do ju = xi-1,xf+1

 uf(1,ju) = u(1,ju)
 uf(imax-1,ju) = u(imax,ju)

end do

if (xi-1 .eq. 1) then

do iu = 1,imax

 u(iu,xi-1) = 0.0d0

end do

do iu = 1,imax-1

 uf(iu,xi-1) = u(iu,xi-1)

end do

end if

if (xf+1 .eq. jmax) then

do iu = 1,imax

 u(iu,xf+1) = lid_velocity
 !u(iu,xf+1) = 0.0d0

end do

do iu = 1,imax-1

 uf(iu,xf+1) = u(iu,xf+1)

end do

end if

! Boundary Conditions for v-velocity

if (xi-1 .eq. 1) then

do iv = 1,imax

 v(iv,xi-1) = 0.0d0
 vf(iv,xi-1) = v(iv,xi-1)

end do

end if

if (xf+1 .eq. jmax) then

do iv = 1,imax

 v(iv,xf+1) = 0.0d0
 vf(iv,xf) = v(iv,xf+1)

end do

end if

do jv = xi-1,xf+1

 v(1,jv) = 0.0d0
 v(imax,jv) = 0.0d0

end do
!write(*,*) 'inside bc', taskid, xi ,xf, shape(vf)
do jv = xi-1, xf+1

 vf(1,jv) = v(1,jv)
 vf(imax,jv) = v(imax,jv)

end do

!! Boundary Conditions for temperature
!if (1 .eq. 0) then
!do jt=xi-1,xf+1
!
! t(1,jt) = t_west
! t(itmax,jt) = t_east
!
!end do
!
!if (xi-1 .eq. 1) then
!
!do it=1,itmax
!
! t(it,xi-1) = t_south
!
!end do
!
!end if
!
!if (xf+1 .eq. jtmax) then
!
!do it=1,itmax
!
! t(it,xf+1) = t_north
!
!end do
!
!end if
!
!end if


end subroutine
