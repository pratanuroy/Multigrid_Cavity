subroutine velocity_coefficients(aw, ae, as, an, ap, ap_zero, d_u, d_v, uf, vf, rho_uf, mu_uf, rho_vf, mu_vf, gammat_uf, gammat_vf, dx, dy, delx, dely, imax, jmax, xi, xf)
use parameters
use parallel_variables
use kind_parameters
implicit none
include "mpif.h"

!! Input variables

integer, intent(in):: imax, jmax, xi, xf
real(kind=dp), dimension(imax,xi-1:xf+1):: aw, ae, as, an, ap, d_u, d_v, ap_zero
real(kind=dp), dimension(imax-1,xi-1:xf+1):: rho_uf, mu_uf, uf, gammat_uf
real(kind=dp), dimension(imax,xi-1:xf+1):: rho_vf, mu_vf, vf, gammat_vf
real(kind=dp), intent(in):: dx(imax), delx(-2:imax+2), dy(jmax), dely(-2:jmax+2)

!! Other local variables
real(kind=dp):: Fn, Fs, Fw, Fe, Dn, Ds, Dw, De, Pn, Ps, Pw, Pe
real(kind=dp):: Ap_Pn, Ap_Ps, Ap_Pw, Ap_Pe
real(kind=dp):: apfun
integer:: iu, ju
real(kind=dp):: ap_sum

! For both u and v control volumes

do ju = xi,xf

      do iu = 2,imax-1

            Fw = rho_uf(iu-1,ju)*uf(iu-1,ju)*dy(ju)
            Dw = mu_uf(iu-1,ju)*dy(ju)/delx(iu-1)
            Pw = Fw/Dw
            Ap_Pw = apfun(Pw)
            aw(iu,ju) = Dw*Ap_Pw + max(Fw,0.0d0)

            Fe = rho_uf(iu,ju)*uf(iu,ju)*dy(ju)
            De = mu_uf(iu,ju)*dy(ju)/delx(iu)
            Pe = Fe/De
            Ap_Pe = apfun(Pe)
            ae(iu,ju) = De*Ap_Pe + max(-Fe,0.0d0)

            Fs = rho_vf(iu,ju-1)*vf(iu,ju-1)*dx(iu)
            Ds = mu_vf(iu,ju-1)*dx(iu)/dely(ju-1)
            Ps = Fs/Ds
            Ap_Ps = apfun(Ps)
            as(iu,ju) = Ds*Ap_Ps + max(Fs,0.0d0)

            Fn = rho_vf(iu,ju)*vf(iu,ju)*dx(iu)
            Dn = mu_vf(iu,ju)*dx(iu)/dely(ju)
            Pn = Fn/Dn
            Ap_Pn = apfun(Pn)
            an(iu,ju) = Dn*Ap_Pn + max(-Fn,0.0d0)

       end do

end do


do ju = xi,xf

      do iu = 2,imax-1

            ap_sum = aw(iu,ju)+ae(iu,ju)+as(iu,ju)+an(iu,ju)
            ap(iu,ju) = ap_sum + ap_zero(iu,ju)
            d_u(iu,ju) = dy(ju)/ap(iu,ju)
            d_v(iu,ju) = dx(iu)/ap(iu,ju)

      end do

end do


tlocal1 = MPI_WTime()

call exchange_data(ap, imax, xi, xf)

tlocal2 = MPI_WTime()

communication_time = communication_time + (tlocal2 - tlocal1)

end subroutine

