subroutine momentum_interpolation_mg(u, v, p, su, sv, s_uf, s_vf, f_uf, f_vf, res_uf, res_vf, aw, ae, as, an, ap, uf, vf, dx, dy, imax, jmax, xi, xf)
use parameters
use kind_parameters
implicit none
!include "mpif.h"

!! Input variables

integer, intent(in):: imax, jmax, xi, xf
real(kind=dp), intent(in), dimension(imax,xi-1:xf+1):: u, v, p, aw, ae, as, an, ap, su, sv
real(kind=dp), intent(in), dimension(imax-1, xi-1:xf+1):: uf, s_uf, res_uf
real(kind=dp), intent(in), dimension(imax, xi-1:xf+1):: vf, s_vf, res_vf
real(kind=dp), intent(in) :: dx(imax), dy(jmax)
real(kind=dp), intent(inout), dimension(imax, xi-1:xf+1):: f_vf
real(kind=dp), intent(inout), dimension(imax-1, xi-1:xf+1):: f_uf
!real(kind=dp):: alphau, alphav

!! Other local variables
real(kind=dp), dimension(imax, xi-1:xf+1):: Hu, Hv
real(kind=dp), dimension(imax-1,xi-1:xf+1):: Huf_bar, auf_bar
real(kind=dp), dimension(imax,xi-1:xf+1):: Hvf_bar, avf_bar
real(kind=dp):: fe_area, fn_area
integer:: iu, ju, iv, jv

do ju = xi,xf

 do iu = 2, imax-1

   Hu(iu,ju) = aw(iu,ju)*u(iu-1,ju) + ae(iu,ju)*u(iu+1,ju) + as(iu,ju)*u(iu,ju-1) + an(iu,ju)*u(iu,ju+1)+ su(iu,ju)
 
   end do

end do

do jv = xi,xf

 do iv = 2, imax-1

   Hv(iv,jv) = aw(iv,jv)*v(iv-1,jv) + ae(iv,jv)*v(iv+1,jv) + as(iv,jv)*v(iv,jv-1) + an(iv,jv)*v(iv,jv+1)+ sv(iv,jv)

 end do

end do

call exchange_data(Hu, imax, xi, xf)
call exchange_data(Hv, imax, xi, xf)

! Internal u control volumes

do ju=xi,xf
      do iu=2,imax-2

            fe_area = dx(iu+1)/(dx(iu+1)+dx(iu))
            auf_bar(iu,ju) = 1.0d0/(fe_area/ap(iu,ju) + (1.0d0-fe_area)/ap(iu+1,ju))
            Huf_bar(iu,ju) = fe_area*Hu(iu,ju) + (1.0d0-fe_area)*Hu(iu+1,ju)
            f_uf(iu,ju) = auf_bar(iu,ju)*uf(iu,ju) - Huf_bar(iu,ju) - res_uf(iu,ju)

      end do

end do

! Internal v control volumes

do jv = xi,xf

      do iv = 2,imax-1

            fn_area = dy(jv+1)/(dy(jv+1)+dy(jv))
            avf_bar(iv,jv) = 1.0d0/(fn_area/ap(iv,jv) + (1.0d0-fn_area)/ap(iv,jv+1))
            Hvf_bar(iv,jv) = fn_area*Hv(iv,jv) + (1.0d0-fn_area)*Hv(iv,jv+1)
            f_vf(iv,jv) = avf_bar(iv,jv)*vf(iv,jv) - Hvf_bar(iv,jv) - res_vf(iv,jv)

      end do

end do

call exchange_data(f_uf, imax-1, xi, xf)
call exchange_data(f_vf, imax, xi, xf)

end subroutine

