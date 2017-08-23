subroutine pcoefficients(uf, vf, auf_bar, avf_bar, awp, aep, asp, anp, app, bp, dx, dy, imax, jmax, xi, xf)
use parallel_variables
use kind_parameters
implicit none

!! Input variables

integer, intent(in):: imax, jmax, xi, xf
real(kind=dp), dimension(imax, xi-1:xf+1), intent(inout):: awp, aep, asp, anp, app
real(kind=dp), dimension(imax, xi-1:xf+1), intent(inout):: bp
real(kind=dp), dimension(imax-1, xi-1:xf+1)::uf, auf_bar
real(kind=dp), dimension(imax, xi-1:xf+1)::vf, avf_bar
real(kind=dp), intent(in):: dx(imax), dy(jmax)

!! Other local variables

real(kind=dp):: totalmb
integer:: ip, jp

!   Velocities are specified along the boundaries

do jp = xi,xf

      do ip = 2,imax-1

            awp(ip,jp) = (dy(jp)/auf_bar(ip-1,jp))*dy(jp)
            aep(ip,jp) = (dy(jp)/auf_bar(ip,jp))*dy(jp)
            asp(ip,jp) = (dx(ip)/avf_bar(ip,jp-1))*dx(ip)
            anp(ip,jp) = (dx(ip)/avf_bar(ip,jp))*dx(ip)

      end do

end do

!  Linking P-coefficients are zeros

if (xi .eq. 2) then

do ip = 2,imax-1

      asp(ip,2) = 0.0d0

end do

end if

if (xf .eq. jmax-1) then

do ip = 2,imax-1

      anp(ip,jmax-1) = 0.0d0

end do

end if

do jp = xi,xf

     awp(2,jp) = 0.0d0
     aep(imax-1,jp) = 0.0d0

end do

do jp = xi,xf

     do ip = 2,imax-1

            app(ip,jp) = aep(ip,jp) + awp(ip,jp) + asp(ip,jp) + anp(ip,jp)

     end do

end do

totalmb = 0.0d0

do jp = xi,xf

     do ip = 2,imax-1

        !bp(ip,jp) = (uf(ip-1,jp)-uf(ip,jp))/dx(ip) + (vf(ip,jp-1)-vf(ip,jp))/dy(jp)

        bp(ip,jp) = (uf(ip-1,jp)-uf(ip,jp))*dy(jp) + (vf(ip,jp-1)-vf(ip,jp))*dx(ip)

        totalmb = totalmb + bp(ip,jp)

     end do

end do

!if (taskid .eq. 0) write(*,*) app
!write(*,*) ''

end subroutine
