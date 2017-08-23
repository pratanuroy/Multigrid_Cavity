subroutine coarse_to_fine_ls(e_2h, nx, ny, xi_2h, xf_2h, xi_h, xf_h, e_h, delx, dely)
!USE variables
use kind_parameters
implicit none
include 'mpif.h'

integer:: i, j, nx, ny, io, jo, ie, je
integer:: xi_h, xf_h, xi_2h, xf_2h
real(kind=dp), intent(in), dimension(-1:nx+2,xi_2h-3:xf_2h+3):: e_2h
real(kind=dp), intent(out), dimension(-1:nx*2,xi_h-3:xf_h+3):: e_h
real(kind=dp), dimension(-2:2*nx):: delx
real(kind=dp), dimension(-2:2*ny):: dely

!write(*,*) nx, xi_2h, xf_2h

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


!do j = xi_2h,xf_2h
!
!    do i = 2,nx-2
!
!    if (j /= ny-1) then
!
!        e_h(2*i,2*j) = (1.0/16.0)*(9.0*e_2h(i+1,j+1) + 3.0*e_2h(i,j+1) + 3.0*e_2h(i+1,j) + e_2h(i,j))
!        e_h(2*i,2*j-1) = (1.0/16.0)*(9.0*e_2h(i+1,j) + 3.0*e_2h(i,j) + 3.0*e_2h(i+1,j+1) + e_2h(i,j+1))
!        e_h(2*i-1,2*j) = (1.0/16.0)*(9.0*e_2h(i,j+1) + 3.0*e_2h(i,j) + 3.0*e_2h(i+1,j+1) + e_2h(i+1,j))
!        e_h(2*i-1,2*j-1) = (1.0/16.0)*(9.0*e_2h(i,j) + 3.0*e_2h(i+1,j) + 3.0*e_2h(i,j+1) + e_2h(i+1,j+1))
!
!    end if
!
!    end do
!
!end do

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


!  e_h(1,xi_h-3:xf_h+3) = e_h(2,xi_h-3:xf_h+3) + delx(1)*(e_h(2,xi_h-3:xf_h+3) - e_h(3,xi_h-3:xf_h+3))/delx(2)
!  e_h(0,xi_h-3:xf_h+3) = e_h(1,xi_h-3:xf_h+3) + delx(0)*(e_h(1,xi_h-3:xf_h+3) - e_h(2,xi_h-3:xf_h+3))/delx(1)
!  e_h(-1,xi_h-3:xf_h+3) = e_h(0,xi_h-3:xf_h+3)+ delx(-1)*(e_h(0,xi_h-3:xf_h+3) - e_h(1,xi_h-3:xf_h+3))/delx(0)
!
!  e_h(nx,xi_h-3:xf_h+3) = e_h(nx-1,xi_h-3:xf_h+3) + delx(nx-1)*(e_h(nx-1,xi_h-3:xf_h+3) - e_h(nx-2,xi_h-3:xf_h+3))/delx(nx-2)
!  e_h(nx+1,xi_h-3:xf_h+3) = e_h(nx,xi_h-3:xf_h+3)+ delx(nx)*(e_h(nx,xi_h-3:xf_h+3) - e_h(nx-1,xi_h-3:xf_h+3))/delx(nx-1)
!  e_h(nx+2,xi_h-3:xf_h+3) = e_h(nx+1,xi_h-3:xf_h+3)+ delx(nx+1)*(e_h(nx+1,xi_h-3:xf_h+3) - e_h(nx,xi_h-3:xf_h+3))/delx(nx)
!
!if (xf_h .eq. ny-1) then
!
!  e_h(:,ny) = e_h(:,ny-1) + dely(ny-1)*(e_h(:,ny-1) - e_h(:,ny-2))/dely(ny-2)
!  e_h(:,ny+1) = e_h(:,ny) + dely(ny)*(e_h(:,ny) - e_h(:,ny-1))/dely(ny-1)
!  e_h(:,ny+2) = e_h(:,ny+1)+ dely(ny+1)*(e_h(:,ny+1) - e_h(:,ny))/dely(ny)
!
!end if
!
!if (xi_h .eq. 2) then
!
!  e_h(:,1) = e_h(:,2)+ dely(1)*(e_h(:,2) - e_h(:,3))/dely(2)
!  e_h(:,0) = e_h(:,1)+ dely(0)*(e_h(:,1) - e_h(:,2))/dely(1)
!  e_h(:,-1) = e_h(:,0)+ dely(-1)*(e_h(:,0) - e_h(:,1))/dely(0)
!
!end if
!
call exchange_lsdata(e_h, nx*2+2, xi_h, xf_h)
!if (xi_2h .eq. 2) then
!
!j = 1
!
!    e_h(2,2*j) = (1.0/4.0)*(e_2h(1,j) + e_2h(1,j+1) + e_2h(2,j) + e_2h(2,j+1))
!
!    do i = 2,nx-2
!
!        e_h(2*i-1,2*j) = (1.0/8.0)*(3.0*e_2h(i,j) + 3.0*e_2h(i,j+1) + e_2h(i+1,j) + e_2h(i+1,j+1))
!        e_h(2*i,2*j) = (1.0/8.0)*(e_2h(i,j) + e_2h(i,j+1) + 3.0*e_2h(i+1,j) + 3.0*e_2h(i+1,j+1))
!
!    end do
!
!    e_h(2*nx-3,2*j) = (1.0/4.0)*(e_2h(nx-1,j) + e_2h(nx-1,j+1) + e_2h(nx,j) + e_2h(nx,j+1))
!
!
!end if
!
!do j = xi_2h,xf_2h
!
!if( j /= ny-1) then
!
!    i = 1
!
!        e_h(2*i,2*j-1) = (1.0/8.0)*(3.0*e_2h(i,j) + 3.0*e_2h(i+1,j) + e_2h(i,j+1) + e_2h(i+1,j+1))
!        e_h(2*i,2*j) = (1.0/8.0)*(e_2h(i,j) + e_2h(i+1,j) + 3.0*e_2h(i,j+1) + 3.0*e_2h(i+1,j+1))
!
!    i = nx-1
!
!        e_h(2*i-1,2*j-1) = (1.0/8.0)*(3.0*e_2h(i,j) + 3.0*e_2h(i+1,j) + e_2h(i,j+1) + e_2h(i+1,j+1))
!        e_h(2*i-1,2*j) = (1.0/8.0)*(e_2h(i,j) + e_2h(i+1,j) + 3.0*e_2h(i,j+1) + 3.0*e_2h(i+1,j+1))
!
!end if
!
!end do
!
!
!if (xf_2h .eq. ny-1) then
!
!j = ny-1
!
!    e_h(2,2*j-1) = (1.0/4.0)*(e_2h(1,j) + e_2h(1,j+1) + e_2h(2,j) + e_2h(2,j+1))
!
!    do i = 2,nx-2
!
!        e_h(2*i-1,2*j-1) = (1.0/8.0)*(3.0*e_2h(i,j) + 3.0*e_2h(i,j+1) + e_2h(i+1,j) + e_2h(i+1,j+1))
!        e_h(2*i,2*j-1) = (1.0/8.0)*(e_2h(i,j) + e_2h(i,j+1) + 3.0*e_2h(i+1,j) + 3.0*e_2h(i+1,j+1))
!
!    end do
!
!    e_h(2*nx-3,2*j-1) = (1.0/4.0)*(e_2h(nx-1,j) + e_2h(nx-1,j+1) + e_2h(nx,j) + e_2h(nx,j+1))
!
!end if

!string = '(   (1x, e10.3))'
!
!write(string(2:4),'(i3)')(nx*2-2)
!
!if (taskid .eq. 0) then
!
!filename = 'e_h0_prev.txt'
!
!open (unit = 13, file = filename)
!
!do jt=xi_h-1,xf_h+1
!
! write(13,string) e_h(1:nx*2-2,jt)
!
!end do
!
!write(13,*) 'xi_h = ',xi_h, 'xf_h= ',xf_h
!write(13,*) 'xi2 = ',xi_2h, 'xf2= ',xf_2h
!
!close(13)
!
!else
!
!filename = 'e_h1_prev.txt'
!
!open (unit = 13, file = filename)
!
!do jt=xi_h-1,xf_h+1
!
! write(13,string) e_h(1:nx*2-2,jt)
!
!end do
!
!write(13,*) 'xi_h = ',xi_h, 'xf_h= ',xf_h
!write(13,*) 'xi2 = ',xi_2h, 'xf2= ',xf_2h
!
!close(13)
!
!end if
!
!call exchange_data(e_h, nx*2-2, xi_h, xf_h)
!
!string = '(   (1x, e10.3))'
!
!write(string(2:4),'(i3)')(nx*2-2)
!
!if (taskid .eq. 0) then
!
!filename = 'e_h0.txt'
!
!open (unit = 13, file = filename)
!
!do jt=xi_h-1,xf_h+1
!
! write(13,string) e_h(1:nx*2-2,jt)
!
!end do
!
!close(13)
!
!else
!
!filename = 'e_h1.txt'
!
!open (unit = 13, file = filename)
!
!do jt=xi_h-1,xf_h+1
!
! write(13,string) e_h(1:nx*2-2,jt)
!
!end do
!
!close(13)
!
!end if

!stop

end subroutine
