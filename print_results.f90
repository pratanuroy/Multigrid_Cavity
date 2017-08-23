subroutine print_results(imax, jmax, xi, xf, u, v, p, uf, vf, xp, yp)
use parameters
use kind_parameters
implicit none

integer, intent(in):: imax, jmax, xi, xf
real(kind=dp), dimension(imax,xi-1:xf+1):: u, v, p
real(kind=dp), dimension(imax-1,xi-1:xf+1):: uf
real(kind=dp), dimension(imax,xi-1:xf+1):: vf
real(kind=dp), dimension(-1:imax+2):: xp
real(kind=dp), dimension(-1:jmax+2):: yp

!character*50:: string, filename
integer:: iu, ju, iv, jv, ip, jp
real(kind=dp):: u_cent_line, v_cent_line
!
!print *, ''
!print *, 'xp = ', xp
!
!print *, ''
!print *, 'yp = ', yp

!print *, imax
!print *, jmax

string = '(   (1x, E10.3))'

if (1 .eq. 1) then

! Write velocities in a txt file

write(string(2:4),'(i3)')imax
filename = 'uvelocity.txt'

open (unit = 11, file = filename)
!
!do ju=1,jmax
!
! write(11,string) u(1:imax,ju)
!
!end do

write(11,*) 'y_position',' ', 'uvelocity'

do ju = xi,xf
    u_cent_line = 0.50d0*(u(imax/2,ju)+u(imax/2+1,ju))/lid_velocity
    write(11,*) yp(ju)/Ly, u_cent_line
end do

! write(11,string) u(1:imax,ju)

!write(11,*) u

close(11)

!write(string(2:3),'(i3)')imax
filename = 'vvelocity.txt'

open (unit = 12, file = filename)

!do jv=1,jmax
!
! write(12,string) v(1:imax,jv)
!
!end do
write(12,*) 'x_position',' ', 'vvelocity'

do jv = xi,xf
do iv = 1,imax
   
    if(jv .eq. jmax/2) then 
    v_cent_line = 0.50d0*(v(iv,jmax/2)+v(iv,jmax/2+1))/lid_velocity
    write(12,*) xp(iv)/Lx, v_cent_line
    end if
end do
end do

!write(12,*) v

close(12)

!write(string(2:3),'(i3)')imax
filename = 'pressure.txt'

open (unit = 14, file = filename)

do jp = 1,jmax

 write(14,string) p(1:imax,jp)

end do

close(14)

filename = 'ufvelocity.txt'

open (unit = 15, file = filename)

do jp = 1,jmax

 write(15,string) uf(1:imax-1,jp)

end do

close(15)

!string = '(  (1x, e10.3))'
!write(string(3:3),'(i1)')(jmax-1)

filename = 'vfvelocity.txt'

open (unit = 16, file = filename)

do jp = 1,jmax-1

 write(16,string) vf(1:imax,jp)

end do

close(16)

end if

end subroutine
