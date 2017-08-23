subroutine unsteady_functions(ap_zero, u, v, f_uvel, f_vvel, aw, ae, as, an, p_uf, p_vf, su, sv, fu, fv, rho, dx, dy, imax, jmax, xi, xf)
use parameters
use kind_parameters
implicit none

integer, intent(in):: imax, jmax, xi, xf
real(kind=dp), intent(in), dimension(imax, xi-1:xf+1):: u, v, rho, aw, ae, as, an
real(kind=dp), intent(inout), dimension(imax, xi-1:xf+1):: su, sv, ap_zero, fu, fv
real(kind=dp), intent(in) :: dx(imax), dy(jmax)
real(kind=dp), dimension(imax, xi-1:xf+1):: f_uvel, f_vvel
real(kind=dp), dimension(imax-1, xi-1:xf+1):: p_uf
real(kind=dp), dimension(imax, xi-1:xf+1):: p_vf
integer:: iu, ju, iv, jv

do ju = xi,xf

    do iu = 2,imax-1

    f_uvel(iu,ju) = (-(aw(iu,ju)+ae(iu,ju)+as(iu,ju)+an(iu,ju))*u(iu,ju)+aw(iu,ju)*u(iu-1,ju)+ae(iu,ju)*u(iu+1,ju)+as(iu,ju)*u(iu,ju-1)+an(iu,ju)*u(iu,ju+1)+dy(ju)*(p_uf(iu-1,ju)-p_uf(iu,ju)) + fu(iu,ju))/(rho(iu,ju)*dy(ju)*dx(iu))

    end do

end do

do jv = xi,xf

    do iv = 2,imax-1

     f_vvel(iv,jv) = (-(aw(iv,jv)+ae(iv,jv)+as(iv,jv)+an(iv,jv))*v(iv,jv)+aw(iv,jv)*v(iv-1,jv)+ae(iv,jv)*v(iv+1,jv)+an(iv,jv)*v(iv,jv+1)+as(iv,jv)*v(iv,jv-1)+dx(iv)*(p_vf(iv,jv-1)-p_vf(iv,jv)) + fv(iv,jv))/(rho(iv,jv)*dx(iv)*dy(jv))

    end do

end do

call exchange_data(f_uvel, imax, xi, xf)
call exchange_data(f_vvel, imax, xi, xf)

do ju = xi,xf

 do iu = 2,imax-1

  ap_zero(iu,ju) = rho(iu,ju)*dx(iu)*dy(ju)/(time_step*A22)

 end do

end do

call exchange_data(ap_zero, imax, xi, xf)

end subroutine
