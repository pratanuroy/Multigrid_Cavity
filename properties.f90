subroutine properties(T, rho, mu, rho_uf, rho_vf, mu_uf, mu_vf, kcond, cp, gammat, imax, jmax, xi, xf)
use parameters
use kind_parameters
implicit none

integer, intent(in):: imax, jmax, xi, xf
real(kind=dp), dimension(imax,xi-1:xf+1):: rho, mu, kcond, cp, gammat
real(kind=dp), dimension(imax-1,xi-1:xf+1):: rho_uf, mu_uf
real(kind=dp), dimension(imax,xi-1:xf+1):: rho_vf, mu_vf
real(kind=dp), intent(in), dimension(-1:imax+2, xi-3:xf+3):: T
integer:: ip, jp
!real(kind=dp):: temp_phi

! Air/ Water

if (prop_flag .eq. 1) then

! Water Properties

do jp=xi-1,xf+1

 do ip=1,imax

  rho(ip,jp) = 997.0d0                       ! density in Kg/m**3
  mu(ip,jp) = 8.71d-4                      ! dynamic viscosity in N. S/m**2
  kcond(ip,jp) = 0.563d0                     ! thermal conductivity in W/m k
  cp(ip,jp) = 4179.0d0                      ! specific heat at constant pressure, J/Kg K
  gammat(ip,jp) = kcond(ip,jp)/cp(ip,jp) 	! gamma temp.(thermal diffusivity), Kg/m. s

  end do

end do

else if (prop_flag .eq. 2) then

! Air Properties

do jp = xi-1,xf+1

 do ip = 1,imax

  rho(ip,jp) = 1.177d0    			! density in Kg/m**3
  mu(ip,jp) = 1.843d-5    			! dynamic viscosity in N. S/m**2
  kcond(ip,jp) = 0.0242d0			!thermal conductivity in w/m k
  cp(ip,jp) = 1006.43d0    			! specific heat at constant pressure, J/Kg K
  gammat(ip,jp) = kcond(ip,jp)/cp(ip,jp) 	! gamma temp.(thermal diffusivity), Kg/m. s

 end do

end do

end if

rho_uf(1:imax-1,xi-1:xf+1) = rho(1:imax-1,xi-1:xf+1)
mu_uf(1:imax-1,xi-1:xf+1) = mu(1:imax-1,xi-1:xf+1)

rho_vf(1:imax,xi-1:xf+1) = rho(1:imax,xi-1:xf+1)
mu_vf(1:imax,xi-1:xf+1) = mu(1:imax,xi-1:xf+1)

end subroutine
