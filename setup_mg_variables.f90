subroutine setup_mg_variables(imax, jmax, q_level)
use multigrid_levels

! input argument variables

integer, intent(in):: q_level, imax, jmax

! local variables

integer:: i, nx, ny, xi, xf

nx = imax
ny = jmax

! main body of the subroutine

!nx = 2**q_level+2
!ny = 2**q_level+2
num_level = q_level

allocate(mg_level(num_level))

do i = num_level, 1, -1

 !write(*,*) i, nx, ny

 mg_level(i)%nx = nx
 mg_level(i)%ny = ny

 call map_task(i)

 xi = mg_level(i)%xi
 xf = mg_level(i)%xf

 write(*,*) 'ilevel = ', i, nx, ny, xi, xf

! Property variables

 allocate(mg_level(i)%rho(nx,xi-1:xf+1))
 allocate(mg_level(i)%mu(nx,xi-1:xf+1))
 allocate(mg_level(i)%kcond(nx,xi-1:xf+1))
 allocate(mg_level(i)%cp(nx,xi-1:xf+1))
 allocate(mg_level(i)%gammat(nx,xi-1:xf+1))

 allocate(mg_level(i)%rho_uf(nx-1,xi-1:xf+1))
 allocate(mg_level(i)%mu_uf(nx-1,xi-1:xf+1))
 allocate(mg_level(i)%rho_vf(nx,xi-1:xf+1))
 allocate(mg_level(i)%mu_vf(nx,xi-1:xf+1))
 allocate(mg_level(i)%gammat_uf(nx,xi-1:xf+1))
 allocate(mg_level(i)%gammat_vf(nx,xi-1:xf+1))

! Velocity variables

 allocate(mg_level(i)%u(nx,xi-1:xf+1))
 allocate(mg_level(i)%v(nx,xi-1:xf+1))

 allocate(mg_level(i)%ap(nx,xi-1:xf+1))
 allocate(mg_level(i)%aw(nx,xi-1:xf+1))
 allocate(mg_level(i)%ae(nx,xi-1:xf+1))
 allocate(mg_level(i)%as(nx,xi-1:xf+1))
 allocate(mg_level(i)%an(nx,xi-1:xf+1))

 allocate(mg_level(i)%bu(nx,xi-1:xf+1))
 allocate(mg_level(i)%bv(nx,xi-1:xf+1))

 allocate(mg_level(i)%d_u(nx,xi-1:xf+1))
 allocate(mg_level(i)%d_v(nx,xi-1:xf+1))

 allocate(mg_level(i)%fu(nx,xi-1:xf+1))
 allocate(mg_level(i)%fv(nx,xi-1:xf+1))

 allocate(mg_level(i)%su(nx,xi-1:xf+1))
 allocate(mg_level(i)%sv(nx,xi-1:xf+1))

  allocate(mg_level(i)%s_uf(nx-1,xi-1:xf+1))
 allocate(mg_level(i)%s_vf(nx,xi-1:xf+1))

 allocate(mg_level(i)%f_uf(nx-1,xi-1:xf+1))
 allocate(mg_level(i)%f_vf(nx,xi-1:xf+1))

 allocate(mg_level(i)%uf(nx-1,xi-1:xf+1))
 allocate(mg_level(i)%vf(nx,xi-1:xf+1))

 allocate(mg_level(i)%auf_bar(nx-1,xi-1:xf+1))
 allocate(mg_level(i)%avf_bar(nx,xi-1:xf+1))

! Pressure correction variables

 allocate(mg_level(i)%p(nx,xi-1:xf+1))
 allocate(mg_level(i)%pc(nx,xi-1:xf+1))
 allocate(mg_level(i)%app(nx,xi-1:xf+1))
 allocate(mg_level(i)%awp(nx,xi-1:xf+1))
 allocate(mg_level(i)%aep(nx,xi-1:xf+1))
 allocate(mg_level(i)%asp(nx,xi-1:xf+1))
 allocate(mg_level(i)%anp(nx,xi-1:xf+1))
 allocate(mg_level(i)%bp(nx,xi-1:xf+1))
 allocate(mg_level(i)%fp(nx,xi-1:xf+1))

 allocate(mg_level(i)%pc_uf(nx-1,xi-1:xf+1))
 allocate(mg_level(i)%p_uf(nx-1,xi-1:xf+1))
 allocate(mg_level(i)%pc_vf(nx,xi-1:xf+1))
 allocate(mg_level(i)%p_vf(nx,xi-1:xf+1))

! Residuals

 allocate(mg_level(i)%resu(nx,xi-1:xf+1))
 allocate(mg_level(i)%resv(nx,xi-1:xf+1))

 allocate(mg_level(i)%res_uf(nx-1,xi-1:xf+1))
 allocate(mg_level(i)%res_vf(nx,xi-1:xf+1))

! For storing previous mg iteration

 allocate(mg_level(i)%uprev(nx,xi-1:xf+1))
 allocate(mg_level(i)%vprev(nx,xi-1:xf+1))

! Temperature variables

! allocate(mg_level(i)%t(nx,ny))
! allocate(mg_level(i)%apt(nx,ny))
! allocate(mg_level(i)%awt(nx,ny))
! allocate(mg_level(i)%aet(nx,ny))
! allocate(mg_level(i)%ast(nx,ny))
! allocate(mg_level(i)%ant(nx,ny))
! allocate(mg_level(i)%bt(nx,ny))
! allocate(mg_level(i)%ft(nx,ny))


! Level set variables

allocate(mg_level(i)%T(-1:nx+2,xi-3:xf+3))
allocate(mg_level(i)%Tprev(-1:nx+2,xi-3:xf+3))
allocate(mg_level(i)%T_exact(1:nx,xi-1:xf+1))

! Geometry variables

 allocate(mg_level(i)%dx(nx))
 allocate(mg_level(i)%dy(ny))
 allocate(mg_level(i)%delx(-2:nx+2))
 allocate(mg_level(i)%dely(-2:ny+2))

 allocate(mg_level(i)%xp(-1:nx+2))
 allocate(mg_level(i)%yp(-1:ny+2))

 allocate(mg_level(i)%error_u(nx,xi-1:xf+1))
 allocate(mg_level(i)%error_v(nx,xi-1:xf+1))

 ! Unsteady variables

 allocate(mg_level(i)%ap_zero(nx,xi-1:xf+1))
 allocate(mg_level(i)%u0(nx,xi-1:xf+1))
 allocate(mg_level(i)%v0(nx,xi-1:xf+1))
 allocate(mg_level(i)%uf0(nx-1,xi-1:xf+1))
 allocate(mg_level(i)%vf0(nx,xi-1:xf+1))
 allocate(mg_level(i)%f_uvel(nx,xi-1:xf+1))
 allocate(mg_level(i)%f_vvel(nx,xi-1:xf+1))

 ! Divide the CV-s into half

 nx = nx/2+1
 ny = ny/2+1

 !write(*,*) shape(mg_level(i)%phi)

end do

end subroutine setup_mg_variables
