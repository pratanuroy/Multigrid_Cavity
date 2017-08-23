subroutine initial_conditions(qlevel)
use global_variables
use multigrid_levels
use parallel_variables
use parameters
use kind_parameters
implicit none

integer, intent(in):: qlevel
integer:: ilevel
integer:: nx, ny, xi, xf


do ilevel = 1, qlevel

!write(*,*) 'ilevel =', ilevel

 nx = mg_level(ilevel)%nx
 ny = mg_level(ilevel)%ny

 xi = mg_level(ilevel)%xi
 xf = mg_level(ilevel)%xf

!! Velocity variables

 u  => mg_level(ilevel)%u
 v  => mg_level(ilevel)%v

 uf  => mg_level(ilevel)%uf
 vf  => mg_level(ilevel)%vf

 aw => mg_level(ilevel)%aw
 ae => mg_level(ilevel)%ae
 as => mg_level(ilevel)%as
 an => mg_level(ilevel)%an
 ap => mg_level(ilevel)%ap

 bu => mg_level(ilevel)%bu
 bv => mg_level(ilevel)%bv

 d_u => mg_level(ilevel)%d_u
 d_v => mg_level(ilevel)%d_v

 fu => mg_level(ilevel)%fu
 fv => mg_level(ilevel)%fv

 uf  => mg_level(ilevel)%uf
 vf  => mg_level(ilevel)%vf

 auf_bar  => mg_level(ilevel)%auf_bar
 avf_bar  => mg_level(ilevel)%avf_bar

 f_uf  => mg_level(ilevel)%f_uf
 f_vf  => mg_level(ilevel)%f_vf

 res_uf  => mg_level(ilevel)%res_uf
 res_vf  => mg_level(ilevel)%res_vf

 resu => mg_level(ilevel)%resu
 resv => mg_level(ilevel)%resv

 su => mg_level(ilevel)%su
 sv => mg_level(ilevel)%sv

 s_uf  => mg_level(ilevel)%s_uf
 s_vf  => mg_level(ilevel)%s_vf

!! Pressure variables

 p  => mg_level(ilevel)%p
 pc  => mg_level(ilevel)%pc
 pc_uf  => mg_level(ilevel)%pc_uf
 pc_vf  => mg_level(ilevel)%pc_vf
 p_uf  => mg_level(ilevel)%p_uf
 p_vf  => mg_level(ilevel)%p_vf

 awp => mg_level(ilevel)%awp
 aep => mg_level(ilevel)%aep
 asp => mg_level(ilevel)%asp
 anp => mg_level(ilevel)%anp
 app => mg_level(ilevel)%app
 fp  => mg_level(ilevel)%fp
 bp  => mg_level(ilevel)%bp

! Temperature Variables

! t => mg_level(ilevel)%t
! awt => mg_level(ilevel)%awt
! aet => mg_level(ilevel)%aet
! ast => mg_level(ilevel)%ast
! ant => mg_level(ilevel)%ant
! apt => mg_level(ilevel)%apt
! ft  => mg_level(ilevel)%ft
! bt  => mg_level(ilevel)%bt

!! Geometry variables

 dx => mg_level(ilevel)%dx
 dy => mg_level(ilevel)%dy

 delx => mg_level(ilevel)%delx
 dely => mg_level(ilevel)%dely

 xp => mg_level(ilevel)%xp
 yp => mg_level(ilevel)%yp

!! Properties variables

 rho => mg_level(ilevel)%rho
 mu => mg_level(ilevel)%mu
 kcond => mg_level(ilevel)%kcond
 cp => mg_level(ilevel)%cp
 gammat => mg_level(ilevel)%gammat

 rho_uf => mg_level(ilevel)%rho_uf
 rho_vf => mg_level(ilevel)%rho_vf
 mu_uf => mg_level(ilevel)%mu_uf
 mu_vf => mg_level(ilevel)%mu_vf
 gammat_uf => mg_level(ilevel)%gammat_uf
 gammat_vf => mg_level(ilevel)%gammat_vf

 !! Unsteady variables

 ap_zero => mg_level(ilevel)%ap_zero
! apt_zero => mg_level(ilevel)%apt_zero
 u0 => mg_level(ilevel)%u0
 v0 => mg_level(ilevel)%v0
 f_uvel => mg_level(ilevel)%f_uvel
 f_vvel => mg_level(ilevel)%f_vvel
! t0 => mg_level(ilevel)%t0
 uf0 => mg_level(ilevel)%uf0
 vf0 => mg_level(ilevel)%vf0


 !! Level Set variables
  T => mg_level(ilevel)%T

! initial conditions

 u(:,:) = 0.0d0
 v(:,:) = 0.0d0
 uf(:,:) = 0.0d0
 vf(:,:) = 0.0d0
 p(:,:) = 0.0d0
 pc(:,:) = 0.0d0

 fu(:,:) = 0.0d0
 fv(:,:) = 0.0d0

 p_uf(:,:) = 0.0d0
 p_vf(:,:) = 0.0d0

 pc_uf(:,:) = 0.0d0
 pc_vf(:,:) = 0.0d0

 fp(:,:) = 0.0d0
 bp(:,:) = 0.0d0

f_uf(:,:) = 0.0d0
f_vf(:,:) = 0.0d0

res_uf(:,:) = 0.0d0
res_vf(:,:) = 0.0d0

resu(:,:) = 0.0d0
resv(:,:) = 0.0d0

su(:,:) = 0.0d0
sv(:,:) = 0.0d0

s_uf(:,:) = 0.0d0
s_vf(:,:) = 0.0d0

f_uvel(:,:) = 0.0d0
f_vvel(:,:) = 0.0d0

!write(*,*) 'info: ', taskid,ilevel, xi, xf
!call print_results(nx, ny, u, v, p, uf, vf)

!gravity = 9.80d0

 call coordinates(xp, yp, dx, dy, nx, ny)

! call setup_ls_func(T, T_exact, grad_phi, delx, dely, dx, dy, xp, yp, nx, ny, xi, xf)
 call properties(T, rho, mu, rho_uf, rho_vf, mu_uf, mu_vf, kcond, cp, gammat, nx, ny, xi, xf)

! Calculate necessary flow parameters

 hydraulic_dia = 4.0*Lx*Ly/(2.0*(Lx+Ly))
!hydraulic_dia = 2.0*Ly

lid_velocity = reynolds_number*mu(1,xi)/(rho(1,xi)*hydraulic_dia)
v_ref = lid_velocity
rho_ref = rho(1,xi)
!v_reference = inlet_velocity

!inlet_velocity = sqrt(gravity*2.0*radius)
!Rayleigh = rho(1,xi)*gravity*beta*(t_west-t_east)*hydraulic_dia**3/(mu(1,xi)*alpha)

!if (taskid .eq. root+1) then
!write(*,*) 'Rayleigh Number = ', Rayleigh
!end if
!write(*,*) 'mu = ', mu(1,xi)
!write(*,*) 'rho = ', rho(1,xi)

 call boundary_conditions(u, v, uf, vf, nx, ny, xi, xf)

!call print_results(nx, ny, u, v, p, uf, vf)

 mg_level(ilevel)%uprev = mg_level(ilevel)%u
 mg_level(ilevel)%vprev = mg_level(ilevel)%v
! mg_level(ilevel)%ufprev = mg_level(ilevel)%uf
! mg_level(ilevel)%vfprev = mg_level(ilevel)%vf

 mg_level(ilevel)%u0 = mg_level(ilevel)%u
 mg_level(ilevel)%v0 = mg_level(ilevel)%v

if (ilevel .eq. num_level) then

allocate(uprint(1:nx,1:ny))
allocate(vprint(1:nx,1:ny))
allocate(pprint(1:nx,1:ny))
!allocate(tprint(1:nx,1:ny))
allocate(tprint(1:nx+4,1:ny+4))
allocate(rho_print(1:nx,1:ny))
allocate(uexact(1:nx,1:ny))
allocate(vexact(1:nx,1:ny))
allocate(u_error(1:nx,1:ny))
allocate(v_error(1:nx,1:ny))

uprint(1,1:ny) = 0.0d0
uprint(nx,1:ny) = 0.0d0
uprint(1:nx,1) = 0.0d0
uprint(1:nx,ny) = lid_velocity

vprint(1:nx,1) = 0.0d0
vprint(1:nx,ny) = 0.0d0
vprint(1,1:ny) = 0.0d0
vprint(nx,1:ny) = 0.0d0

!tprint(1:itmax,1) = t_south
!tprint(1:itmax,jtmax) = t_north
!tprint(1,1:jtmax) = t_west
!tprint(itmax,1:jtmax) = t_east

u0 = u
v0 = v
!t0 = t

end if

 end do

! Set number of iteration at each MG level

allocate(iter_num(qlevel))

iter_num(qlevel) = 2
!iter_num(qlevel-1) = 3
!iter_num(qlevel-2) = 6
!iter_num(qlevel-3) = 7

 do ilevel = qlevel,2,-1

!iter_num(ilevel-1) = iter_num(ilevel) + 2
iter_num(ilevel-1) = 2

end do

!iter_num(qlevel) = 1
iter_num(1) = 40

if (taskid .eq. root) write(*,*) 'iter_num', iter_num

end subroutine initial_conditions
