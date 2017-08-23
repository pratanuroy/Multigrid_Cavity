subroutine calculate_face_velocity_interp(u, v, uf, vf, p, p_uf, p_vf, s_uf, s_vf, ap, dx, dy, nx, ny, xi, xf)
use kind_parameters
implicit none

integer, intent(in):: nx, ny, xi, xf
real(kind=dp), intent(in), dimension(nx,xi-1:xf+1):: u, v, ap, p
real(kind=dp), intent(in):: dx(nx), dy(ny)
real(kind=dp), intent(in), dimension(nx-1,xi-1:xf+1):: p_uf, s_uf
real(kind=dp), intent(in), dimension(nx,xi-1:xf+1):: p_vf, s_vf
real(kind=dp), intent(inout), dimension(nx-1,xi-1:xf+1):: uf
real(kind=dp), intent(inout), dimension(nx,xi-1:xf+1):: vf

!! Other local variables
real(kind=dp), dimension(nx-1, xi-1:xf+1):: auf_bar
real(kind=dp), dimension(nx, xi-1:xf+1):: avf_bar
real(kind=dp):: fe_area, fn_area
integer:: iu, ju, iv, jv

do ju=xi,xf
      do iu=2,nx-2

            fe_area = dx(iu+1)/(dx(iu+1)+dx(iu))
            auf_bar(iu,ju) = 1.0d0/(fe_area/ap(iu,ju) + (1.0d0-fe_area)/ap(iu+1,ju))

      end do

end do

do ju = xi,xf

    do iu = 2,nx-2

        uf(iu,ju) = 0.50d0*(u(iu,ju) + u(iu+1,ju)) + 1.0d0/auf_bar(iu,ju)*((p(iu,ju) - p(iu+1,ju))*dy(ju) + s_uf(iu,ju)) - 0.50d0*(1.0d0/ap(iu,ju)*(p_uf(iu-1,ju)-p_uf(iu,ju))*dy(ju) + 1.0d0/ap(iu+1,ju)*(p_uf(iu,ju)-p_uf(iu+1,ju))*dy(ju))
!        uf(iu,ju) = 0.50d0*(u(iu,ju) + u(iu+1,ju)) + 1.0d0/auf_bar(iu,ju)*((p(iu,ju) - p(iu+1,ju))*dy(ju)) - 0.50d0*(1.0d0/ap(iu,ju)*(p_uf(iu-1,ju)-p_uf(iu,ju))*dy(ju) + 1.0d0/ap(iu+1,ju)*(p_uf(iu,ju)-p_uf(iu+1,ju))*dy(ju))
!
    end do

end do

do jv = xi,xf

      do iv = 2,nx-1

            fn_area = dy(jv+1)/(dy(jv+1)+dy(jv))
            avf_bar(iv,jv) = 1.0d0/(fn_area/ap(iv,jv) + (1.0d0-fn_area)/ap(iv,jv+1))

      end do

end do

do jv = xi,xf

    do iv = 2,nx-1

        vf(iv,jv) = 0.50d0*(v(iv,jv) + v(iv,jv+1)) + 1.0d0/avf_bar(iv,jv)*((p(iv,jv) - p(iv,jv+1))*dx(iv) + s_vf(iv,jv)) - 0.50d0*(1.0d0/ap(iv,jv)*(p_vf(iv,jv-1)-p_vf(iv,jv))*dx(iv) + 1.0d0/ap(iv,jv+1)*(p_vf(iv,jv)-p_vf(iv,jv+1))*dx(iv))
!        vf(iv,jv) = 0.50d0*(v(iv,jv) + v(iv,jv+1)) + 1.0d0/avf_bar(iv,jv)*((p(iv,jv) - p(iv,jv+1))*dx(iv)) - 0.50d0*(1.0d0/ap(iv,jv)*(p_vf(iv,jv-1)-p_vf(iv,jv))*dx(iv) + 1.0d0/ap(iv,jv+1)*(p_vf(iv,jv)-p_vf(iv,jv+1))*dx(iv))

    end do

end do

if (xf .eq. ny-1) then

    do iv = 2,nx-1

        vf(iv,ny-1) = v(iv,ny)

    end do

end if

call exchange_data(uf, nx-1, xi, xf)
call exchange_data(vf, nx, xi, xf)

end subroutine
