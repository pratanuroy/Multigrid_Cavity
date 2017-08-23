subroutine correction(u, v, uf, vf, pc, p, auf_bar, avf_bar, p_uf, p_vf, d_u, d_v, dx, dy, imax, jmax, xi, xf, alphap)
        !use variables
        use parallel_variables
use kind_parameters
        implicit none

        !! Input variables

        integer, intent(in):: imax, jmax, xi, xf
        real(kind=dp), dimension(imax, xi-1:xf+1):: u, v, pc, p, d_u, d_v
        real(kind=dp), dimension(imax-1, xi-1:xf+1):: uf, auf_bar, p_uf, pc_uf
        real(kind=dp), dimension(imax, xi-1:xf+1):: vf, avf_bar, p_vf, pc_vf
        real(kind=dp):: dx(imax), dy(jmax)
        real(kind=dp):: alphap

        !! Other local variables
        integer:: iu, ju, ip, jp, iv, jv
        real(kind=dp):: fe_area, fn_area

        !write(*,*) 'alphap = ', alphap

        ! Calculate boundary pressure correction values

        if (xi .eq. 2) then

                do ip=2,imax-1

                !     pc(ip,1) = (15.0*pc(ip,2) - 10.0*pc(ip,3) + 3.0*pc(ip,4))/8.0d0
                pc(ip,1) = pc(ip,2)

                end do

        end if

        if (xf .eq. jmax-1) then

                do ip=2,imax-1

                !       pc(ip,jmax) = (15.0*pc(ip,jmax-1) - 10.0*pc(ip,jmax-2) + 3.0*pc(ip,jmax-3))/8.0d0
                pc(ip,jmax) = pc(ip,jmax-1)

                end do

        end if

        do jp=xi,xf

        !    pc(1,jp) = (15.0*pc(2,jp) - 10.0*pc(3,jp) + 3.0*pc(4,jp))/8.0d0
        !    pc(imax,jp) = (15.0*pc(imax-1,jp) - 10.0*pc(imax-2,jp) + 3.0*pc(imax-3,jp))/8.0d0
        pc(1,jp) = pc(2,jp)
        pc(imax,jp) = pc(imax-1,jp)

        end do

        call exchange_data(pc, imax, xi, xf)

        ! Correct face velocities

        do ju = xi,xf

        do iu = 2,imax-2

        uf(iu,ju) = uf(iu,ju) + dy(ju)/auf_bar(iu,ju)*(pc(iu,ju)-pc(iu+1,ju))

        end do

        end do

        do jv = xi,xf

        do iv=2,imax-1

        vf(iv,jv) = vf(iv,jv) + dx(iv)/avf_bar(iv,jv)*(pc(iv,jv)-pc(iv,jv+1))

        end do

        end do

        if (xf .eq. jmax-1) then

                do iv = 2,imax-1

                vf(iv,jmax-1) = v(iv,jmax)

                end do

        end if

        ! Correct nodal pressures

        do jp = xi,xf

        do ip = 2,imax-1

        p(ip,jp) = p(ip,jp) + alphap*pc(ip,jp)

        end do

        end do


        if (xi .eq. 2) then

                do ip=2,imax-1

                p(ip,1) = (15.0*p(ip,2) - 10.0*p(ip,3) + 3.0*p(ip,4))/8.0d0
                !     p(ip,1) = p(ip,2)

                end do

        end if

        if (xf .eq. jmax-1) then

                do ip=2,imax-1

                p(ip,jmax) = (15.0*p(ip,jmax-1) - 10.0*p(ip,jmax-2) + 3.0*p(ip,jmax-3))/8.0d0
                !       p(ip,jmax) = p(ip,jmax-1)

                end do

        end if

        do jp=xi,xf

        p(1,jp) = (15.0*p(2,jp) - 10.0*p(3,jp) + 3.0*p(4,jp))/8.0d0
        p(imax,jp) = (15.0*p(imax-1,jp) - 10.0*p(imax-2,jp) + 3.0*p(imax-3,jp))/8.0d0
        !    p(1,jp) = p(2,jp)
        !    p(imax,jp) = p(imax-1,jp)

        end do

        call exchange_data(p, imax, xi, xf)

        ! Calculate pressure and pressure correction values at interface

        do jp = xi,xf

        do ip = 2,imax-1

        fe_area = dx(ip+1)/(dx(ip+1)+dx(ip))
        p_uf(ip,jp) = fe_area*p(ip,jp) + (1.0d0-fe_area)*p(ip+1,jp)
        pc_uf(ip,jp) = fe_area*pc(ip,jp) + (1.0d0-fe_area)*pc(ip+1,jp)

        end do

        end do

        do jp =xi,xf

        do ip = 2,imax-1

        fn_area = dy(jp+1)/(dy(jp+1)+dy(jp))
        p_vf(ip,jp) = fn_area*p(ip,jp) + (1.0d0-fn_area)*p(ip,jp+1)
        pc_vf(ip,jp) = fn_area*pc(ip,jp) + (1.0d0-fn_area)*pc(ip,jp+1)

        end do

        end do

        pc_uf(1,xi-1:xf+1) = pc(1,xi-1:xf+1)
        p_uf(1,xi-1:xf+1) = p(1,xi-1:xf+1)

        pc_uf(imax-1,xi-1:xf+1) = pc(imax,xi-1:xf+1)
        p_uf(imax-1,xi-1:xf+1) = p(imax,xi-1:xf+1)

        !pc_uf(1,xi-1:xf+1) = pc_uf(2,xi-1:xf+1)
        !p_uf(1,xi-1:xf+1) = p_uf(2,xi-1:xf+1)

        !pc_uf(imax-1,xf-1:xf+1) = pc_uf(imax-2,xf-1:xf+1)
        !p_uf(imax-1,xf-1:xf+1) = p_uf(imax-2,xf-1:xf+1)

        if (xi .eq. 2) then

                pc_vf(1:imax,1) = pc(1:imax,1)
                p_vf(1:imax,1) = p(1:imax,1)

                !pc_vf(1:imax,1) = pc_vf(2:imax,1)
                !p_vf(1:imax,1) = p_vf(2:imax,1)

        end if

        if (xf .eq. jmax-1) then

                pc_vf(1:imax,jmax-1) = pc(1:imax,jmax)
                p_vf(1:imax,jmax-1) = p(1:imax,jmax)

                !pc_vf(1:imax,jmax-1) = pc_vf(1:imax,jmax-2)
                !p_vf(1:imax,jmax-1) = p_vf(1:imax,jmax-2)

        end if

        call exchange_data(pc_uf, imax-1, xi, xf)
        call exchange_data(pc_vf, imax, xi, xf)
        call exchange_data(p_uf, imax-1, xi, xf)
        call exchange_data(p_vf, imax, xi, xf)

        ! Correct nodal velocities using pressure correction at interfaces

        do ju = xi,xf

        do iu = 2,imax-1

        u(iu,ju) = u(iu,ju) + d_u(iu,ju)*(pc_uf(iu-1,ju)-pc_uf(iu,ju))

        end do

        end do

        do jv = xi,xf

        do iv=2,imax-1

        v(iv,jv) = v(iv,jv) + d_v(iv,jv)*(pc_vf(iv,jv-1)-pc_vf(iv,jv))

        end do

        end do

        uf(1,xi-1:xf+1) = u(1,xi-1:xf+1)
        uf(imax-1,xi-1:xf+1) = u(imax,xi-1:xf+1)

        vf(1,xi-1:xf+1) = v(1,xi-1:xf+1)
        vf(imax,xi-1:xf+1) = v(imax,xi-1:xf+1)

        call exchange_data(uf, imax-1, xi, xf)
        call exchange_data(vf, imax, xi, xf)
        call exchange_data(u, imax, xi, xf)
        call exchange_data(v, imax, xi, xf)

end subroutine

