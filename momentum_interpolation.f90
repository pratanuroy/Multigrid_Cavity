subroutine momentum_interpolation(u, v, p, aw, ae, as, an, ap, su, sv, f_uf, f_vf, s_uf, s_vf, res_uf, res_vf, auf_bar, avf_bar, uf, vf, dx, dy, imax, jmax, xi, xf)
use parameters
use kind_parameters
implicit none
include "mpif.h"

!! Input variables

integer, intent(in):: imax, jmax, xi, xf
real(kind=dp), intent(in), dimension(imax,xi-1:xf+1):: u, v, p, aw, ae, as, an, ap, su, sv
real(kind=dp), intent(inout), dimension(imax-1, xi-1:xf+1):: uf, auf_bar, f_uf, res_uf, s_uf
real(kind=dp), intent(inout), dimension(imax, xi-1:xf+1):: vf, avf_bar, f_vf, res_vf, s_vf
real(kind=dp):: dx(imax), dy(jmax)
!real(kind=dp):: alphau, alphav

!! Other local variables
real(kind=dp), dimension(imax, xi-1:xf+1):: Hu, Hv
real(kind=dp), dimension(imax-1,xi-1:xf+1):: Huf_bar
real(kind=dp), dimension(imax,xi-1:xf+1):: Hvf_bar
real(kind=dp):: fe_area, fn_area
integer:: iu, ju, iv, jv
!real(kind=dp):: temp_fu, temp_fv

do ju = xi,xf

 do iu = 2, imax-1

  Hu(iu,ju) = aw(iu,ju)*u(iu-1,ju) + ae(iu,ju)*u(iu+1,ju) + as(iu,ju)*u(iu,ju-1) + an(iu,ju)*u(iu,ju+1) + su(iu,ju)

 end do

end do


do jv = xi,xf

 do iv = 2, imax-1

  Hv(iv,jv) = aw(iv,jv)*v(iv-1,jv) + ae(iv,jv)*v(iv+1,jv) + as(iv,jv)*v(iv,jv-1) + an(iv,jv)*v(iv,jv+1) + sv(iv,jv)

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
            uf(iu,ju) = (1.0d0-alphau)*uf(iu,ju) + (alphau)/auf_bar(iu,ju)*(Huf_bar(iu,ju) + (p(iu,ju) - p(iu+1,ju))*dy(ju) + f_uf(iu,ju))

      end do

end do

do ju=xi,xf

      do iu=2,imax-2

            res_uf(iu,ju) = auf_bar(iu,ju)*uf(iu,ju) - Huf_bar(iu,ju) - (p(iu,ju) - p(iu+1,ju))*dy(ju) - f_uf(iu,ju)

      end do
end do


! Internal v control volumes

do jv = xi,xf

      do iv = 2,imax-1

            fn_area = dy(jv+1)/(dy(jv+1)+dy(jv))
            avf_bar(iv,jv) = 1.0d0/(fn_area/ap(iv,jv) + (1.0d0-fn_area)/ap(iv,jv+1))
            Hvf_bar(iv,jv) = fn_area*Hv(iv,jv) + (1.0d0-fn_area)*Hv(iv,jv+1)
            vf(iv,jv) = (1.0d0-alphav)*vf(iv,jv) + (alphav)/avf_bar(iv,jv)*(Hvf_bar(iv,jv) + (p(iv,jv) - p(iv,jv+1))*dx(iv)+ f_vf(iv,jv))
      end do

end do

if (xf .eq. jmax-1) then

      do iv = 2,imax-1

           vf(iv,jmax-1) = v(iv,jmax)

      end do

end if

!vf(1,xi-1:xf+1) = vf(2,xi-1:xf+1)
!vf(imax,xi-1:xf+1) = vf(imax-1,xi-1:xf+1)
!
!uf(1,xi-1:xf+1) = uf(2,xi-1:xf+1)
!uf(imax-1,xi-1:xf+1) = uf(imax-2,xi-1:xf+1)


do jv=xi,xf

      do iv=2,imax-1

            res_vf(iv,jv) = avf_bar(iv,jv)*vf(iv,jv) - Hvf_bar(iv,jv) - (p(iv,jv) - p(iv,jv+1))*dx(iv) - f_vf(iv,jv)

      end do
end do

call exchange_data(auf_bar, imax-1, xi, xf)
call exchange_data(avf_bar, imax, xi, xf)
call exchange_data(uf, imax-1, xi, xf)
call exchange_data(vf, imax, xi, xf)
call exchange_data(res_uf, imax-1, xi, xf)
call exchange_data(res_vf, imax, xi, xf)

end subroutine

