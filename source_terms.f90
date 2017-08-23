subroutine source_terms(u0, v0, uf0, vf0, ap_zero, f_uvel, f_vvel, su, sv, bu, bv,  s_uf, s_vf, rho_vf, p_uf, p_vf, dx, dy, imax, jmax, xi, xf)
use parameters
use kind_parameters
implicit none
include "mpif.h"

!! Input variables

integer, intent(in):: imax, jmax, xi, xf
real(kind=dp), dimension(imax,xi-1:xf+1):: u0, v0, f_uvel, f_vvel, bu, bv, ap_zero, su, sv
real(kind=dp), dimension(imax-1,xi-1:xf+1):: p_uf, s_uf, uf0
real(kind=dp), dimension(imax,xi-1:xf+1):: p_vf, rho_vf, s_vf, vf0
real(kind=dp), intent(in):: dx(imax), dy(jmax)

!! Other local variables
integer:: iu, ju, iv, jv
real(kind=dp):: temp_rho

!f_uvel(:,:) = 0.0d0
!f_vvel(:,:) = 0.0d0


! For both u and v control volumes
! For both u and v control volumes

do ju=xi,xf

      do iu=2,imax-1

            temp_rho = 0.50d0*(rho_vf(iu,ju-1) + rho_vf(iu,ju))
            su(iu,ju) = ap_zero(iu,ju)*(u0(iu,ju) + f_uvel(iu,ju)*time_step*A21)
            !bv(iu,ju) = ap_zero(iu,ju)*(v0(iu,ju) + f_vvel(iu,ju)*time_step*A21) - surface_tension*kappa(iu,ju)*nv_y(iu,ju)*dx(iu)*dy(ju)- rho(iu,ju)*gravity*dx(iu)*dy(ju)

            !bu(iu,ju) = ap_zero(iu,ju)*(u0(iu,ju) + f_uvel(iu,ju)*time_step*A21)
            !bv(iu,ju) = ap_zero(iu,ju)*(v0(iu,ju) + f_vvel(iu,ju)*time_step*A21) + rho(iu,ju)*gravity*beta*(t(iu,ju) - tambient)*dx(iu)*dy(ju)
            !sv(iu,ju) = ap_zero(iu,ju)*(v0(iu,ju) + f_vvel(iu,ju)*time_step*A21) - temp_rho*gravity*dx(iu)*dy(ju)
            !bv(iu,ju) = ap_zero(iu,ju)*(v0(iu,ju) + f_vvel(iu,ju)*time_step*A21)

            !su(iu,ju) = ap_zero(iu,ju)*(u0(iu,ju))
            bu(iu,ju) = (p_uf(iu-1,ju)-p_uf(iu,ju))*dy(ju) + su(iu,ju)

            !sv(iu,ju) = ap_zero(iu,ju)*(v0(iu,ju))  + rho_vf(iu,ju)*gravity*beta*(t(iu,ju) - tambient)*dx(iu)*dy(ju)
            sv(iu,ju) = ap_zero(iu,ju)*(v0(iu,ju)+ f_vvel(iu,ju)*time_step*A21)  - temp_rho*gravity*dx(iu)*dy(ju)
            !sv(iu,ju) = ap_zero(iu,ju)*(v0(iu,ju))
            !bv(iu,ju) = (p_vf(iu,ju-1)-p_vf(iu,ju))*dx(iu) + sv(iu,ju) - temp_rho*gravity*dx(iu)*dy(ju)
            bv(iu,ju) = (p_vf(iu,ju-1)-p_vf(iu,ju))*dx(iu) + sv(iu,ju)

       end do

end do

do ju = xi,xf

 do iu = 2,imax-2

   !s_uf(iu,ju) = 0.50d0*(ap_zero(iu,ju)*u0(iu,ju) + ap_zero(iu+1,ju)*u0(iu+1,ju))
   s_uf(iu,ju) = 0.50d0*(ap_zero(iu,ju) + ap_zero(iu+1,ju))*uf0(iu,ju)

 end do

end do

do jv = xi,xf

 do iv = 2,imax-1

   s_vf(iv,jv) = 0.50d0*(ap_zero(iv,jv)+ ap_zero(iv,jv+1))*vf0(iv,jv+1) - rho_vf(iv,jv)*gravity*dx(iv)*dy(jv)

 end do

end do

if (xf .eq. jmax-1) then

    s_vf(:,jmax-1) = sv(:,jmax)

end if

call exchange_data(su, imax, xi, xf)
call exchange_data(sv, imax, xi, xf)
call exchange_data(s_uf, imax-1, xi, xf)
call exchange_data(s_vf, imax, xi, xf)

end subroutine

