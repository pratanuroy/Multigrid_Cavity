subroutine setup_ls_func(T, T_exact, grad_phi, delx, dely, dx, dy, xp, yp, nx, ny, xi, xf)
use parallel_variables
use parameters
use kind_parameters
implicit none
include "mpif.h"

! Input variables

integer, intent(in):: nx, ny, xi, xf
real(kind=dp), intent(inout), dimension(-1:nx+2, xi-3:xf+3):: T
real(kind=dp), intent(inout), dimension(1:nx, xi-1:xf+1):: T_exact, grad_phi
real(kind=dp), intent(in):: delx(-2:nx+2), dely(-2:ny+2), dx(1:nx), dy(1:ny)
real(kind=dp), intent(in):: xp(-1:nx+2), yp(-1:ny+2)
! Local variables

integer:: ip, jp, it, jt, i, j
real(kind=dp):: center_x1, center_x2, center_y1, center_y2, radius1, radius2, slot_width, slot_length, temp1
real(kind=dp):: local_vol_exact

! Initialize the level-set function, T(x,y)

!radius = 0.20*pi
!center_x = 0.50*pi
!center_y = 0.20*(pi+1)
!L_perimeter = 2.0*pi*radius
!vol_exact = pi*radius*radius
!L_perimeter = 4*(2.0*radius)
!vol_exact = 4.0*radius*radius

radius = 1.2d-2
center_x = Lx/2.0d0
center_y = 6.0d-2
!radius = 1.0d-2
!center_x = Lx/2.0d0
!center_y = 1.2d-1

! Zalesak's disk
!radius = 15.0d0
!slot_width = 5.0d0
!slot_length = 25.0d0
!center_x = 50.0d0
!center_y = 75.0d0


! For merging of two bubble

!center_x1 = 0.5d0
!center_y1 = 0.35d0
!radius1 = 0.1d0

!center_x2 = 0.5d0
!center_y2 = 0.65d0
!radius2 = 0.15d0

!radius = radius1

L_perimeter = 2.0*pi*radius
vol_exact = pi*radius*radius

!inlet_velocity = reynolds_number*mu(1,1)/(rho(1,1)*hydraulic_dia)
v_ref = dsqrt(gravity*radius)
time_ref = v_ref/radius

rho_ref = rho_l
mu_ref = mu_l

reynolds_number = rho_ref*v_ref*radius/mu_ref
bond_number = rho_ref*gravity*radius*radius/surface_tension

write(*,*) 'Reference Velocity = ', v_ref
write(*,*) 'gravity = ', gravity
write(*,*) 'epsilon = ', epsilon

do jt = xi-3,xf+3

 do it = -1,nx+2

	! For Bubble drop
    !T(it,jt) = ((xp(it) - center_x)*(xp(it) - center_x) + (yp(jt) - center_y)*(yp(jt) - center_y)) - radius*radius
    T(it,jt) = dsqrt((xp(it) - center_x)*(xp(it) - center_x) + (yp(jt) - center_y)*(yp(jt) - center_y)) - radius
    !T(it,jt) = min(yp(jt)-0.05d0,dsqrt((xp(it) - center_x)*(xp(it) - center_x) + (yp(jt) - center_y)*(yp(jt) - center_y)) - radius)

    ! Merging of two bubbles
    !T(it,jt) = min(dsqrt((xp(it) - center_x1)**2 + (yp(jt) -center_y1)**2) - radius1, dsqrt((xp(it) - center_x2)*(xp(it) - center_x2) + (yp(jt) - center_y2)*(yp(jt) - center_y2)) - radius2)

    !For Rayleigh Taylor instability
    !T(it,jt) = yp(jt) -(2.0d0 + 0.05d0*cos(2.0d0*pi*xp(it)/1.0d0))

    !For dam breaking
    !T(it,jt) = max(yp(jt)-1.0d0,xp(it)-1.0d0)

    !For capillary rise
    !T(it,jt) = yp(jt) - 7.0832d-4

    ! Zalesak's disk
!    temp1 = max(abs(yp(jt)-(75.0d0-2.50d0))-25.0d0/2.0d0, abs(xp(it)-50.0d0)-5.0d0/2.0d0)
!    T(it,jt) = min(radius - dsqrt((xp(it) - center_x)*(xp(it) - center_x) + (yp(jt) - center_y)*(yp(jt) - center_y)), temp1)

 end do

end do

!! Boundary conditions for level-set function
!
!  T(1,xi-3:xf+3) = T(2,xi-3:xf+3) + delx(1)*(T(2,xi-3:xf+3) - T(3,xi-3:xf+3))/delx(2)
!  T(0,xi-3:xf+3) = T(1,xi-3:xf+3) + delx(0)*(T(1,xi-3:xf+3) - T(2,xi-3:xf+3))/delx(1)
!  T(-1,xi-3:xf+3) = T(0,xi-3:xf+3)+ delx(-1)*(T(0,xi-3:xf+3) - T(1,xi-3:xf+3))/delx(0)
!
!  T(nx,xi-3:xf+3) = T(nx-1,xi-3:xf+3) + delx(nx-1)*(T(nx-1,xi-3:xf+3) - T(nx-2,xi-3:xf+3))/delx(nx-2)
!  T(nx+1,xi-3:xf+3) = T(nx,xi-3:xf+3)+ delx(nx)*(T(nx,xi-3:xf+3) - T(nx-1,xi-3:xf+3))/delx(nx-1)
!  T(nx+2,xi-3:xf+3) = T(nx+1,xi-3:xf+3)+ delx(nx+1)*(T(nx+1,xi-3:xf+3) - T(nx,xi-3:xf+3))/delx(nx)
!
!if (xf .eq. ny-1) then
!
!  T(:,ny) = T(:,ny-1) + dely(ny-1)*(T(:,ny-1) - T(:,ny-2))/dely(ny-2)
!  T(:,ny+1) = T(:,ny) + dely(ny)*(T(:,ny) - T(:,ny-1))/dely(ny-1)
!  T(:,ny+2) = T(:,ny+1)+ dely(ny+1)*(T(:,ny+1) - T(:,ny))/dely(ny)
!
!end if
!
!if (xi .eq. 2) then
!
!  T(:,1) = T(:,2)+ dely(1)*(T(:,2) - T(:,3))/dely(2)
!  T(:,0) = T(:,1)+ dely(0)*(T(:,1) - T(:,2))/dely(1)
!  T(:,-1) = T(:,0)+ dely(-1)*(T(:,0) - T(:,1))/dely(0)
!
!end if

tlocal1 = MPI_Wtime()

call exchange_lsdata(T, nx+4, xi, xf)

tlocal2 = MPI_Wtime()

communication_time = communication_time + (tlocal2 - tlocal1)

!call kappa_func1

! Calculate the exact area/ volume
vol_exact = 0.0d0
local_vol_exact = 0.0d0
T_exact = T
do jp = xi,xf

 do ip = 1,nx

     call Heaviside(T_exact(ip,jp), H_phi_exact, epsilon)
     local_vol_exact = local_vol_exact + (1.0d0 - H_phi_exact)*dx(ip)*dy(jp)

 end do

end do

call MPI_Allreduce(local_vol_exact, vol_exact, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

grad_phi(:,:) = 0.0d0
grad_phi(2:nx-1,xi:xf) = 1.0d0
!error_check(:,:) = 0.0d0
end subroutine
