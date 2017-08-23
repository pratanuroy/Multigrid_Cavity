subroutine convergencet(t, atw, ate, ats, atn, atp, bt, itmax, jtmax, residualt)
!use variables
use kind_parameters
implicit none

integer, intent(in):: itmax, jtmax
real(kind=dp), dimension(itmax,jtmax):: t, atw, ate, ats, atn, atp, bt
real(kind=dp):: temp1, temp2, sum_nom, sum_den
integer:: it, jt
real(kind=dp), intent(out):: residualt

sum_nom = 0.0d0
sum_den = 0.0d0

do jt=2,jtmax-1

       do it=2,itmax-1

              temp1 = atp(it,jt)*t(it,jt)
              temp2 = atw(it,jt)*t(it-1,jt) + ate(it,jt)*t(it+1,jt) + ats(it,jt)*t(it,jt-1) + atn(it,jt)*t(it,jt+1)
              sum_nom = sum_nom + dabs(temp1 - temp2 - bt(it,jt))
              sum_den = sum_den + dabs(temp1)

       end do

end do
	
residualt = sum_nom/sum_den

end subroutine
