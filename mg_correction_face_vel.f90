subroutine mg_correction_face_vel(uf, vf, pc, auf_bar, avf_bar, dx, dy, imax, jmax)
!use variables
use kind_parameters
implicit none

!! Input variables

integer, intent(in):: imax, jmax
real(kind=dp), dimension(imax, jmax):: pc
real(kind=dp), dimension(imax-1, jmax):: uf, auf_bar
real(kind=dp), dimension(imax, jmax-1):: vf, avf_bar
real(kind=dp):: dx(imax), dy(jmax)

!! Other local variables
integer:: iu, ju, iv, jv

! Correct face velocities

do ju = 2,jmax-1

      do iu = 2,imax-2

            uf(iu,ju) = uf(iu,ju) + dy(ju)/auf_bar(iu,ju)*(pc(iu,ju)-pc(iu+1,ju))

      end do

end do

do jv = 2,jmax-2

      do iv=2,imax-1

            vf(iv,jv) = vf(iv,jv) + dx(iv)/avf_bar(iv,jv)*(pc(iv,jv)-pc(iv,jv+1))

      end do

end do

end subroutine
