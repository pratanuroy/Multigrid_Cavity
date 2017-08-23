subroutine deallocate_mg_variables()
use multigrid_levels

! To do: Make deallocation more robust and concise

! input argument variables

!integer, intent(in):: num_level

! local variables

integer:: i

! main body of the subroutine

do i = num_level, 1, -1

! Property variables

 deallocate(mg_level(i)%rho)
 deallocate(mg_level(i)%mu)
 deallocate(mg_level(i)%kcond)
 deallocate(mg_level(i)%cp)
 deallocate(mg_level(i)%gammat)

 deallocate(mg_level(i)%rho_uf)
 deallocate(mg_level(i)%mu_uf)
 deallocate(mg_level(i)%rho_vf)
 deallocate(mg_level(i)%mu_vf)
 deallocate(mg_level(i)%gammat_uf)
 deallocate(mg_level(i)%gammat_vf)

! Velocity variables

 deallocate(mg_level(i)%u)
 deallocate(mg_level(i)%v)

 deallocate(mg_level(i)%ap)
 deallocate(mg_level(i)%aw)
 deallocate(mg_level(i)%ae)
 deallocate(mg_level(i)%as)
 deallocate(mg_level(i)%an)

 deallocate(mg_level(i)%bu)
 deallocate(mg_level(i)%bv)

 deallocate(mg_level(i)%d_u)
 deallocate(mg_level(i)%d_v)

 deallocate(mg_level(i)%fu)
 deallocate(mg_level(i)%fv)

 deallocate(mg_level(i)%su)
 deallocate(mg_level(i)%sv)

  deallocate(mg_level(i)%s_uf)
 deallocate(mg_level(i)%s_vf)

 deallocate(mg_level(i)%f_uf)
 deallocate(mg_level(i)%f_vf)

 deallocate(mg_level(i)%uf)
 deallocate(mg_level(i)%vf)

 deallocate(mg_level(i)%auf_bar)
 deallocate(mg_level(i)%avf_bar)

! Pressure correction variables

 deallocate(mg_level(i)%p)
 deallocate(mg_level(i)%pc)
 deallocate(mg_level(i)%app)
 deallocate(mg_level(i)%awp)
 deallocate(mg_level(i)%aep)
 deallocate(mg_level(i)%asp)
 deallocate(mg_level(i)%anp)
 deallocate(mg_level(i)%bp)
 deallocate(mg_level(i)%fp)

 deallocate(mg_level(i)%pc_uf)
 deallocate(mg_level(i)%p_uf)
 deallocate(mg_level(i)%pc_vf)
 deallocate(mg_level(i)%p_vf)

! Residuals

 deallocate(mg_level(i)%resu)
 deallocate(mg_level(i)%resv)

 deallocate(mg_level(i)%res_uf)
 deallocate(mg_level(i)%res_vf)

! For storing previous mg iteration

 deallocate(mg_level(i)%uprev)
 deallocate(mg_level(i)%vprev)

! Temperature variables

! deallocate(mg_level(i)%t)
! deallocate(mg_level(i)%apt)
! deallocate(mg_level(i)%awt)
! deallocate(mg_level(i)%aet)
! deallocate(mg_level(i)%ast)
! deallocate(mg_level(i)%ant)
! deallocate(mg_level(i)%bt)
! deallocate(mg_level(i)%ft)

! Level set variables

deallocate(mg_level(i)%T)

! Geometry variables

 deallocate(mg_level(i)%dx)
 deallocate(mg_level(i)%dy)
 deallocate(mg_level(i)%delx)
 deallocate(mg_level(i)%dely)

 deallocate(mg_level(i)%xp)
 deallocate(mg_level(i)%yp)

 deallocate(mg_level(i)%error_u)
 deallocate(mg_level(i)%error_v)

 ! Unsteady variables

 deallocate(mg_level(i)%ap_zero)
 deallocate(mg_level(i)%u0)
 deallocate(mg_level(i)%v0)
 deallocate(mg_level(i)%uf0)
 deallocate(mg_level(i)%vf0)
 deallocate(mg_level(i)%f_uvel)
 deallocate(mg_level(i)%f_vvel)

 !write(*,*) shape(mg_level(i)%phi)

end do

deallocate(mg_level)

end subroutine deallocate_mg_variables
