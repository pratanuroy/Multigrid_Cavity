subroutine line_gs_solver(phi, ap_phi, aw_phi, ae_phi, as_phi, an_phi, b_phi, f_phi, nx, ny, alpha, xi, xf)
use kind_parameters
implicit none

integer, intent(in):: nx, ny, xi, xf
real(kind=dp), intent(inout), dimension(nx,xi-1:xf+1):: phi
real(kind=dp), intent(in), dimension(nx,xi-1:xf+1):: ap_phi, aw_phi, ae_phi, as_phi, an_phi, b_phi, f_phi
real(kind=dp), intent(in):: alpha

integer:: i, j, N
real(kind=dp):: alpha_inv
real(kind=dp), allocatable:: a(:,:), b(:)
!real(kind=dp):: dx(nx), dy(ny), delx(-2:nx+2), dely(-2:ny+2)

alpha_inv = 1.0d0/alpha

!	Sweep from bottom to top
!	Implicit in X

N = nx-2
allocate(a(N,3), b(N))

a(1,1) = 0.0d0
a(nx-2,3) = 0.0d0

do j = xi,xf

        a(1,2) = ap_phi(2,j)*alpha_inv
        a(1,3) = -ae_phi(2,j)
        b(1) = as_phi(2,j)*phi(2,j-1) + an_phi(2,j)*phi(2,j+1) + phi(2,j)*ap_phi(2,j)*(alpha_inv-1.0d0) + aw_phi(2,j)*phi(1,j) + b_phi(2,j) + f_phi(2,j)

        do i=3,nx-2

            a(i-1,1) = -aw_phi(i,j)
            a(i-1,2) = ap_phi(i,j)*alpha_inv
            a(i-1,3) = -ae_phi(i,j)
            b(i-1) = an_phi(i,j)*phi(i,j+1) + as_phi(i,j)*phi(i,j-1) + ap_phi(i,j)*(alpha_inv-1.0d0)*phi(i,j) + b_phi(i,j) + f_phi(i,j)

        end do

        a(nx-2,1) = -aw_phi(nx-1,j)
        a(nx-2,2) = ap_phi(nx-1,j)*alpha_inv
        b(nx-2) = as_phi(nx-1,j)*phi(nx-1,j-1) + an_phi(nx-1,j)*phi(nx-1,j+1) + ap_phi(nx-1,j)*phi(nx-1,j)*(alpha_inv-1.0d0) + ae_phi(nx-1,j)*phi(nx,j) + b_phi(nx-1,j) + f_phi(nx-1,j)

        call tdma(a, b, N)

        do i = 1,N

            phi(i+1,j) = b(i)

        end do

end do

call exchange_data(phi, nx, xi, xf)

!	Sweep from top to bottom
!	Implicit in X


do j = xf,xi,-1

        a(1,2) = ap_phi(2,j)*alpha_inv
        a(1,3) = -ae_phi(2,j)
        b(1) = as_phi(2,j)*phi(2,j-1) + an_phi(2,j)*phi(2,j+1) + phi(2,j)*ap_phi(2,j)*(alpha_inv-1.0d0) + aw_phi(2,j)*phi(1,j) + b_phi(2,j)+ f_phi(2,j)

        do i=3,nx-2

            a(i-1,1) = -aw_phi(i,j)
            a(i-1,2) = ap_phi(i,j)*alpha_inv
            a(i-1,3) = -ae_phi(i,j)
            b(i-1) = an_phi(i,j)*phi(i,j+1) + as_phi(i,j)*phi(i,j-1) + ap_phi(i,j)*(alpha_inv-1.0d0)*phi(i,j) + b_phi(i,j) + f_phi(i,j)

        end do

        a(nx-2,1) = -aw_phi(nx-1,j)
        a(nx-2,2) = ap_phi(nx-1,j)*alpha_inv
        b(nx-2) = as_phi(nx-1,j)*phi(nx-1,j-1) + an_phi(nx-1,j)*phi(nx-1,j+1) + ap_phi(nx-1,j)*phi(nx-1,j)*(alpha_inv-1.0d0) + ae_phi(nx-1,j)*phi(nx,j) + b_phi(nx-1,j) + f_phi(nx-1,j)

        call tdma(a, b, N)

        do i = 1,N

            phi(i+1,j) = b(i)

        end do

end do


deallocate(a,b)

call exchange_data(phi, nx, xi, xf)

end subroutine
