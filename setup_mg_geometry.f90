subroutine setup_mg_geometry
 use global_variables
 use multigrid_levels
 use parameters
use kind_parameters
 implicit none

 integer :: imax, jmax
 integer :: ilevel
 integer :: it, jt

 do ilevel = num_level,1, -1

 imax = mg_level(ilevel)%nx
 jmax = mg_level(ilevel)%ny

! write(*,*) 'I am here'
! write(*,*) ilevel, imax, jmax

 ! Calculate control volumes

 mg_level(ilevel)%dx(1) = 0.0d0

 do it = 2,imax-1

  mg_level(ilevel)%dx(it) = Lx/(imax-2)

 end do

 mg_level(ilevel)%dx(imax) = 0.0d0

 mg_level(ilevel)%dy(1) = 0.0d0

 do jt = 2,jmax-1

  mg_level(ilevel)%dy(jt) = Ly/(jmax-2)

 end do

 mg_level(ilevel)%dy(jmax) = 0.0d0

 ! Calculate diffusion lengths

 mg_level(ilevel)%delx(1) = 0.5d0*mg_level(ilevel)%dx(2)

 do it = 2,imax-2

  mg_level(ilevel)%delx(it) = 0.5d0*(mg_level(ilevel)%dx(it)+mg_level(ilevel)%dx(it+1))

 end do

 mg_level(ilevel)%delx(imax-1) = 0.5d0*mg_level(ilevel)%dx(imax-1)

 mg_level(ilevel)%delx(0) =  mg_level(ilevel)%delx(1)
 mg_level(ilevel)%delx(-1) = mg_level(ilevel)%delx(2)
 mg_level(ilevel)%delx(-2) = mg_level(ilevel)%delx(3)
 mg_level(ilevel)%delx(imax) = mg_level(ilevel)%delx(imax-1)
 mg_level(ilevel)%delx(imax+1) = mg_level(ilevel)%delx(imax-2)
 mg_level(ilevel)%delx(imax+2) = mg_level(ilevel)%delx(imax-3)

 mg_level(ilevel)%dely(1) = 0.5d0*mg_level(ilevel)%dy(2)

 do jt = 2,jmax-2

  mg_level(ilevel)%dely(jt) = 0.5d0*(mg_level(ilevel)%dy(jt) + mg_level(ilevel)%dy(jt+1))

 end do

 mg_level(ilevel)%dely(jmax-1) = 0.5d0*mg_level(ilevel)%dy(jmax-1)

 mg_level(ilevel)%dely(0) =  mg_level(ilevel)%dely(1)
 mg_level(ilevel)%dely(-1) = mg_level(ilevel)%dely(2)
 mg_level(ilevel)%dely(-2) = mg_level(ilevel)%dely(3)
 mg_level(ilevel)%dely(jmax) = mg_level(ilevel)%dely(jmax-1)
 mg_level(ilevel)%dely(jmax+1) = mg_level(ilevel)%dely(jmax-2)
 mg_level(ilevel)%dely(jmax+2) = mg_level(ilevel)%dely(jmax-3)

! write(*,*) ''
! write(*,*) mg_level(ilevel)%dx
! write(*,*) ''
! write(*,*) mg_level(ilevel)%dy
! write(*,*) ''
! write(*,*) mg_level(ilevel)%delx
! write(*,*) ''
! write(*,*) mg_level(ilevel)%dely
! write(*,*) ''

 end do

end subroutine setup_mg_geometry
