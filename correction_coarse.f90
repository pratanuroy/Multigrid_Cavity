subroutine correction_coarse(u, v, uf, vf, pc, p, auf_bar, avf_bar, p_uf, p_vf, d_u, d_v, dx, dy, imax, jmax, alphap)
!use variables
use kind_parameters
implicit none

!! Input variables

integer, intent(in):: imax, jmax
real(kind=dp), dimension(imax, jmax):: u, v, pc, p, d_u, d_v
real(kind=dp), dimension(imax-1, jmax):: uf, auf_bar, p_uf, pc_uf
real(kind=dp), dimension(imax, jmax-1):: vf, avf_bar, p_vf, pc_vf
real(kind=dp):: dx(imax), dy(jmax)
real(kind=dp):: alphap

!! Other local variables
integer:: iu, ju, ip, jp, iv, jv
real(kind=dp):: fe_area, fn_area
!real(kind=dp):: sum_num, sum_den
!real(kind=dp):: massflux_in, massflux_out, massflux_out_nom, massflux_out_den

!write(*,*) 'alphap = ', alphap

! Correct face velocities

do ju = 2,jmax-1

      do iu = 2,imax-2

            uf(iu,ju) = uf(iu,ju) + dy(ju)/auf_bar(iu,ju)*(pc(iu,ju)-pc(iu+1,ju))

      end do

end do

do jv = 2,jmax-2

      do iv=2,imax-1

            vf(iv,jv) = vf(iv,jv) + dx(iv)/avf_bar(iv,jv)*(pc(iv,jv)-pc(iv,jv+1))

      end do

end do

do jp = 2,jmax-1

     do ip = 2,imax-1

            p(ip,jp) = p(ip,jp) + alphap*pc(ip,jp)

     end do

end do

do jp=2,jmax-1

!    p(1,jp) = (15.0*p(2,jp) - 10.0*p(3,jp) + 3.0*p(4,jp))/8.0d0
!    p(imax,jp) = (15.0*p(imax-1,jp) - 10.0*p(imax-2,jp) + 3.0*p(imax-3,jp))/8.0d0
    !p(1,jp) = p(2,jp) + delx(1)*(p(2,jp) - p(3,jp))/delx(2)
    !p(imax,jp) = p(imax-1,jp) + delx(imax-1)*(p(imax-1,jp) - p(imax-2,jp))/delx(imax-2)
    p(1,jp) = p(2,jp)
    p(imax,jp) = p(imax-1,jp)

    !pc(1,jp) = (15.0*pc(2,jp) - 10.0*pc(3,jp) + 3.0*pc(4,jp))/8.0d0
    !pc(imax,jp) = (15.0*pc(imax-1,jp) - 10.0*pc(imax-2,jp) + 3.0*pc(imax-3,jp))/8.0d0

    !pc(1,jp) = pc(2,jp) + delx(1)*(pc(2,jp) - pc(3,jp))/delx(2)
    !pc(1,jp) = ((xp(1) - xp(2))*pc(3,jp) - (xp(1) - xp(3))*pc(2,jp))/(xp(3) - xp(2))
    !pc(imax,jp) = pc(imax-1,jp) + delx(imax-1)*(pc(imax-1,jp) - pc(imax-2,jp))/delx(imax-2)
    pc(1,jp) = pc(2,jp)
    pc(imax,jp) = pc(imax-1,jp)

end do

do ip=2,imax-1

!      p(ip,1) = (15.0*p(ip,2) - 10.0*p(ip,3) + 3.0*p(ip,4))/8.0d0
!      p(ip,jmax) = (15.0*p(ip,jmax-1) - 10.0*p(ip,jmax-2) + 3.0*p(ip,jmax-3))/8.0d0
     !p(ip,1) = p(ip,2) + dely(2)*(p(ip,2) - p(ip,3))/dely(2)
     !p(ip,1) = ((yp(1) - yp(2))*p(ip,3) - (yp(1) - yp(3))*p(ip,2))/(yp(3) - yp(2))
     !p(ip,jmax) = ((yp(jmax) - yp(jmax-1))*p(ip,jmax-2) - (yp(jmax) - yp(jmax-2))*p(ip,jmax-1))/(yp(jmax-2) - yp(jmax-1))
     !p(ip,jmax) = p(ip,jmax-1) + dely(jmax-1)*(p(ip,jmax-1) - p(ip,jmax-2))/dely(jmax-2)
     p(ip,1) = p(ip,2)
     p(ip,jmax) = p(ip,jmax-1)

     !pc(ip,1) = (15.0*pc(ip,2) - 10.0*pc(ip,3) + 3.0*pc(ip,4))/8.0d0
     !pc(ip,jmax) = (15.0*pc(ip,jmax-1) - 10.0*pc(ip,jmax-2) + 3.0*pc(ip,jmax-3))/8.0d0
     !pc(ip,1) = ((yp(1) - yp(2))*pc(ip,3) - (yp(1) - yp(3))*pc(ip,2))/(yp(3) - yp(2))
     !pc(ip,jmax) = ((yp(jmax) - yp(jmax-1))*pc(ip,jmax-2) - (yp(jmax) - yp(jmax-2))*pc(ip,jmax-1))/(yp(jmax-1) - yp(jmax-2))

     !pc(ip,1) = pc(ip,2) + dely(1)*(pc(ip,2) - pc(ip,3))/dely(2)
     !pc(ip,jmax) = pc(ip,jmax-1) + dely(jmax-1)*(pc(ip,jmax-1) - pc(ip,jmax-2))/dely(jmax-2)
     pc(ip,1) = pc(ip,2)
     pc(ip,jmax) = pc(ip,jmax-1)

end do


! Calculate pressure and pressure correction values at interface

do jp = 2,jmax-1

     do ip = 2,imax-1

            fe_area = dx(ip+1)/(dx(ip+1)+dx(ip))
            p_uf(ip,jp) = fe_area*p(ip,jp) + (1.0d0-fe_area)*p(ip+1,jp)
            pc_uf(ip,jp) = fe_area*pc(ip,jp) + (1.0d0-fe_area)*pc(ip+1,jp)

      end do

end do

do jp = 2,jmax-1

    do ip = 2,imax-1

            fn_area = dy(jp+1)/(dy(jp+1)+dy(jp))
            p_vf(ip,jp) = fn_area*p(ip,jp) + (1.0d0-fn_area)*p(ip,jp+1)
            pc_vf(ip,jp) = fn_area*pc(ip,jp) + (1.0d0-fn_area)*pc(ip,jp+1)

    end do

end do

!pc_uf(imax-1,:) = pc_uf(imax-2,:)
!pc_uf(1,:) = pc_uf(2,:)
!pc_vf(:,jmax-1) = pc_vf(:,jmax-2)
!pc_vf(:,1) = pc_vf(:,2)

!pc_uf(imax-1,:) = pc(imax,:)
!pc_uf(1,:) = pc(1,:)
!pc_uf(1,:) = pc_uf(2,:)
!pc_uf(imax-1,:) = pc_uf(imax-2,:)
!pc_uf(1,:) = pc_uf(2,:) + dx(2)*(pc_uf(2,:) - pc_uf(3,:))/dx(3)
!p(imax,jp) = p(imax-1,jp) + delx(imax-1)*(p(imax-1,jp) - p(imax-2,jp))/delx(imax-2)
!pc_vf(:,jmax-1) = pc(:,jmax)
pc_vf(:,1) = pc(:,1)
p_vf(:,1) = p(:,1)
!pc_vf(:,1) = pc_vf(:,2)
!pc_vf(:,jmax-1) = pc_vf(:,jmax-2)
!pc_vf(:,1) = pc_vf(:,2) + dy(2)*(pc_vf(:,2) - pc_vf(:,3))/dy(3)

!p_uf(imax-1,:) = p(imax,:)
p_uf(1,:) = p(1,:)
pc_uf(1,:) = pc(1,:)
!p_uf(1,:) = p_uf(2,:) + dx(2)*(p_uf(2,:) - p_uf(3,:))/dx(3)
!p_vf(:,jmax-1) = p(:,jmax)
!p_vf(:,1) = p(:,1)
!p_vf(:,1) = p_vf(:,2) + dy(2)*(p_vf(:,2) - p_vf(:,3))/dy(3)


! Correct nodal velocities using pressure correction at interfaces

do ju = 2,jmax-1

      do iu = 2,imax-1

            u(iu,ju) = u(iu,ju) + d_u(iu,ju)*(pc_uf(iu-1,ju)-pc_uf(iu,ju))

      end do

end do

do jv = 2,jmax-1

      do iv=2,imax-1

	     v(iv,jv) = v(iv,jv) + d_v(iv,jv)*(pc_vf(iv,jv-1)-pc_vf(iv,jv))

      end do

end do

!u(1,:) = u(2,:)
!u(iumax,:) = u(iumax-1,:)
!
!uf(1,:) = uf(2,:)
!uf(iumax-1,:) = uf(iumax-2,:)

!uf(1,:) = u(1,:)
!uf(iumax-1,:) = u(iumax,:)

!v(1,:) = v(2,jv)
!v(ivmax,:) = v(ivmax-1,jv)

!vf(1,1:jvmax-1) = v(1,1:jvmax-1)
!vf(ivmax,1:jvmax-1) = v(ivmax,1:jvmax-1)

!vf(1,1:jvmax-1) = vf(2,1:jvmax-1)
!vf(ivmax,1:jvmax-1) = vf(ivmax-1,1:jvmax-1)

end subroutine

