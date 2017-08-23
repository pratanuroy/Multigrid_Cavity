subroutine fine_to_coarse_ls(r_h, nx, ny, xi_h, xf_h, xi_2h, xf_2h, r_2h, delx, dely)
use kind_parameters
implicit none

integer:: i, j, nx, ny
integer:: xi_2h, xf_2h, xi_h, xf_h
real(kind=dp), intent(in), dimension(-1:nx+2,xi_h-3:xf_h+3):: r_h
real(kind=dp), intent(out), dimension(-1:nx/2+3,xi_2h-3:xf_2h+3):: r_2h
real(kind=dp), dimension(-2:nx/2+3):: delx
real(kind=dp), dimension(-2:ny/2+3):: dely
integer:: imax, jmax

! Apply full weighting to get interpolated values from fine to coarse grid
! Input
! r_h   Variable in fine grid
! f     Right Hand Side
! hx, hy     grid spacing
!
! Output
! r_2h     Variable in coarse grid

! Author: Pratanu Roy
! History:
! First Written: July 20, 2012
! Indices Rechecked on Jan 13, 2013

imax = nx/2
jmax = ny/2

do j = xi_2h, xf_2h
    do i = 2,nx/2

        r_2h(i,j) = 1.0/4.0*(r_h(2*i-2,2*j-2) + r_h(2*i-2,2*j-1) + r_h(2*i-1,2*j-2) + r_h(2*i-1,2*j-1))

     end do

end do

!  r_2h(1,xi_2h-3:xf_2h+3) = r_2h(2,xi_2h-3:xf_2h+3) + delx(1)*(r_2h(2,xi_2h-3:xf_2h+3) - r_2h(3,xi_2h-3:xf_2h+3))/delx(2)
!  r_2h(0,xi_2h-3:xf_2h+3) = r_2h(1,xi_2h-3:xf_2h+3) + delx(0)*(r_2h(1,xi_2h-3:xf_2h+3) - r_2h(2,xi_2h-3:xf_2h+3))/delx(1)
!  r_2h(-1,xi_2h-3:xf_2h+3) = r_2h(0,xi_2h-3:xf_2h+3)+ delx(-1)*(r_2h(0,xi_2h-3:xf_2h+3) - r_2h(1,xi_2h-3:xf_2h+3))/delx(0)
!
!  r_2h(imax,xi_2h-3:xf_2h+3) = r_2h(imax-1,xi_2h-3:xf_2h+3) + delx(imax-1)*(r_2h(imax-1,xi_2h-3:xf_2h+3) - r_2h(imax-2,xi_2h-3:xf_2h+3))/delx(imax-2)
!  r_2h(imax+1,xi_2h-3:xf_2h+3) = r_2h(imax,xi_2h-3:xf_2h+3)+ delx(imax)*(r_2h(imax,xi_2h-3:xf_2h+3) - r_2h(imax-1,xi_2h-3:xf_2h+3))/delx(imax-1)
!  r_2h(imax+2,xi_2h-3:xf_2h+3) = r_2h(imax+1,xi_2h-3:xf_2h+3)+ delx(imax+1)*(r_2h(imax+1,xi_2h-3:xf_2h+3) - r_2h(imax,xi_2h-3:xf_2h+3))/delx(imax)
!
!if (xf_2h .eq. jmax-1) then
!
!  r_2h(:,jmax) = r_2h(:,jmax-1) + dely(jmax-1)*(r_2h(:,jmax-1) - r_2h(:,jmax-2))/dely(jmax-2)
!  r_2h(:,jmax+1) = r_2h(:,jmax) + dely(jmax)*(r_2h(:,jmax) - r_2h(:,jmax-1))/dely(jmax-1)
!  r_2h(:,jmax+2) = r_2h(:,jmax+1)+ dely(jmax+1)*(r_2h(:,jmax+1) - r_2h(:,jmax))/dely(jmax)
!
!end if
!
!if (xi_2h .eq. 2) then
!
!  r_2h(:,1) = r_2h(:,2)+ dely(1)*(r_2h(:,2) - r_2h(:,3))/dely(2)
!  r_2h(:,0) = r_2h(:,1)+ dely(0)*(r_2h(:,1) - r_2h(:,2))/dely(1)
!  r_2h(:,-1) = r_2h(:,0)+ dely(-1)*(r_2h(:,0) - r_2h(:,1))/dely(0)
!
!end if

call exchange_lsdata(r_2h, (nx/2)+5, xi_2h, xf_2h)

end subroutine
