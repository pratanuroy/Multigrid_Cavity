subroutine calc_residualt(phi, atp, atw, ate, ats, atn, bt, imax, jmax, dx, dy, residualt)
use kind_parameters
implicit none

integer :: it, jt, imax, jmax
real(kind=dp), intent(inout) :: residualt
real(kind=dp), dimension(imax,jmax), intent(in) :: atp, atw, ate, ats, atn, bt
real(kind=dp), dimension(imax,jmax), intent(in) :: phi
real(kind=dp), dimension(imax), intent(in) :: dx
real(kind=dp), dimension(jmax), intent(in) :: dy
real(kind=dp) :: tempsum1

tempsum1 = 0.0d0


!write(*,*) imax, jmax
!write(*,*) size(phi)

do jt = 2,jmax-1

 do it = 2,imax-1

  tempsum1 = tempsum1 + abs(atp(it,jt)*phi(it,jt) - atw(it,jt)*phi(it-1,jt) - ate(it,jt)*phi(it+1,jt) - ats(it,jt)*phi(it,jt-1) - atn(it,jt)*phi(it,jt+1) - bt(it,jt)*dx(it)*dy(jt))

  !if(bt(it,jt) /= 0.0d0) write(*,*) 'bt not zero', it, jt
 end do

end do

residualt = tempsum1

end subroutine
