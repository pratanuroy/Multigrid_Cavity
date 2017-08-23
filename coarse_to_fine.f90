subroutine coarse_to_fine(e_2h, nx, ny, xi_2h, xf_2h, xi_h, xf_h, e_h)
use parallel_variables
use kind_parameters
implicit none
!include 'mpif.h'

integer, intent(in):: nx, ny, xi_h, xf_h, xi_2h, xf_2h
real(kind=dp), intent(in), dimension(nx,xi_2h-1:xf_2h+1):: e_2h
real(kind=dp), intent(out), dimension(nx*2-2,xi_h-1:xf_h+1):: e_h
!character*50:: filename, string
integer:: i, j, io, jo, ie, je

! Apply linear interpolation to get values from coarse to fine grid
! Input
! e_2h   Variable in coarse grid
! hx, hy     grid spacing
!
! Output
! e_h     Variable in fine grid

! Author: Pratanu Roy
! History:
! First Written: July 20, 2012
! Indices Rechecked on Jan 13, 2013

do j = xi_h,xf_h

    do i = 2,2*nx-3

    io = (i+1)/2
    jo = (j+1)/2

    ie = i/2
    je = j/2

    if (i .eq. 2) then

        if (j .eq. 2) then

        e_h(i,j) =  (1.0/4.0)*(e_2h(ie,je) + e_2h(ie,je+1) + e_2h(ie+1,je) + e_2h(ie+1,je+1))

        else if (j .eq. 2*ny-3) then

        e_h(i,j) = (1.0/4.0)*(e_2h(ie,jo) + e_2h(ie,jo+1) + e_2h(ie+1,jo) + e_2h(ie+1,jo+1))

        else

            if(mod(j,2) .eq. 0) then

            e_h(i,j) = (1.0/8.0)*(3.0*e_2h(ie,je+1) + 3.0*e_2h(ie+1,je+1) + e_2h(ie,je) + e_2h(ie+1,je))

            else

            e_h(i,j) = (1.0/8.0)*(e_2h(ie,jo+1) + e_2h(ie+1,jo+1) + 3.0*e_2h(ie,jo) + 3.0*e_2h(ie+1,jo))

            end if

        end if

   else if (i .eq. 2*nx-3) then

        if (j .eq. 2) then

        e_h(i,j) =  (1.0/4.0)*(e_2h(io,je) + e_2h(io,je+1) + e_2h(io+1,je) + e_2h(io+1,jo+1))

        else if (j .eq. 2*ny-3) then

        e_h(i,j) = (1.0/4.0)*(e_2h(io,jo) + e_2h(io,jo+1) + e_2h(io+1,jo) + e_2h(io+1,jo+1))

        else

            if(mod(j,2) .eq. 0) then

            e_h(i,j) = (1.0/8.0)*(3.0*e_2h(io,je+1) + 3.0*e_2h(io+1,je+1) + e_2h(io,je) + e_2h(io+1,je))

            else

            e_h(i,j) = (1.0/8.0)*(e_2h(io,jo+1) + e_2h(io+1,jo+1) + 3.0*e_2h(io,jo) + 3.0*e_2h(io+1,jo))

            end if

        end if

    else

        if (j .eq. 2) then

            if(mod(i,2) .eq. 0) then

            e_h(i,j) = (1.0/8.0)*(3.0*e_2h(ie+1,je) + 3.0*e_2h(ie+1,je+1) + e_2h(ie,je) + e_2h(ie,je+1))

            else

            e_h(i,j) = (1.0/8.0)*(e_2h(io+1,je) + e_2h(io+1,je+1) + 3.0*e_2h(io,je) + 3.0*e_2h(io,je+1))

            end if

        else if (j .eq. 2*ny-3) then

            if(mod(i,2) .eq. 0) then

            e_h(i,j) = (1.0/8.0)*(3.0*e_2h(ie+1,jo+1) + 3.0*e_2h(ie+1,jo) + e_2h(ie,jo) + e_2h(ie,jo+1))

            else

            e_h(i,j) = (1.0/8.0)*(3.0*e_2h(io,jo+1) + e_2h(io+1,jo+1) + 3.0*e_2h(io,jo) + e_2h(io+1,jo))

            end if

        else

            if (mod(i,2) .ne. 0 .and. mod(j,2) .ne. 0) e_h(i,j) = (1.0/16.0)*(9.0*e_2h(io,jo) + 3.0*e_2h(io,jo+1) + 3.0*e_2h(io+1,jo) + e_2h(io+1,jo+1))
            if (mod(i,2) .eq. 0 .and. mod(j,2) .eq. 0) e_h(i,j) = (1.0/16.0)*(e_2h(ie,je) + 3.0*e_2h(ie,je+1) + 3.0*e_2h(ie+1,je) + 9.0*e_2h(ie+1,je+1))
            if (mod(i,2) .ne. 0 .and. mod(j,2) .eq. 0) e_h(i,j) = (1.0/16.0)*(3.0*e_2h(io,je) + 3.0*e_2h(io+1,je+1) + e_2h(io+1,je) + 9.0*e_2h(io,je+1))
            if (mod(i,2) .eq. 0 .and. mod(j,2) .ne. 0) e_h(i,j) = (1.0/16.0)*(3.0*e_2h(ie,jo) + 3.0*e_2h(ie+1,jo+1) + 9.0*e_2h(ie+1,jo) + e_2h(ie,jo+1))

        end if

    end if

    end do

end do

call exchange_data(e_h, nx*2-2, xi_h, xf_h)

end subroutine
