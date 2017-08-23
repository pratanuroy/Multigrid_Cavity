subroutine tcoefficients(t, awt, aet, ast, ant, apt, bt, uf, vf, rho_uf, rho_vf, gammat_uf, gammat_vf, dx, dy, delx, dely, imax, jmax)
!use variables
use kind_parameters
implicit none

integer, intent(in):: imax, jmax
real(kind=dp), intent(in), dimension(imax,jmax):: t, apt, awt, aet, ast, ant, bt
real(kind=dp), intent(in), dimension(imax-1, jmax):: uf, rho_uf, gammat_uf
real(kind=dp), intent(in), dimension(imax, jmax-1):: vf, rho_vf, gammat_vf
real(kind=dp), intent(in) :: dx(imax), dy(jmax), delx(imax-1), dely(jmax-1)

integer:: it, jt
real(kind=dp):: Fn, Fs, Fw, Fe, Dn, Ds, Dw, De, Pn, Ps, Pw, Pe
real(kind=dp):: Ap_Pn, Ap_Ps, Ap_Pw, Ap_Pe
real(kind=dp):: apfun

do jt=2,jmax-1

    do it=2,imax-1

        Fw = rho_uf(it-1,jt)*dy(jt)*uf(it-1,jt)
        Dw = gammat_uf(it-1,jt)*dy(jt)/delx(it-1)
        Pw = Fw/Dw
        Ap_Pw = apfun(Pw)
        awt(it,jt) = Dw*Ap_Pw + max(Fw,0.0)

        Fe = rho_uf(it,jt)*dy(jt)*uf(it,jt)
        De = gammat_uf(it,jt)*dy(jt)/delx(it)
        Pe = Fe/De
        Ap_Pe = apfun(Pe)
        aet(it,jt) = De*Ap_Pe+ max(-Fe,0.0)

        Fs = rho_vf(it,jt-1)*dx(it)*vf(it,jt-1)
        Ds = gammat_vf(it,jt-1)*dx(it)/dely(jt-1)
        Ps = Fs/Ds
        Ap_Ps = apfun(Ps)
        ast(it,jt) = Ds*Ap_Ps+max(Fs,0.0)

        Fn = rho_vf(it,jt)*dx(it)*vf(it,jt)
        Dn = gammat_vf(it,jt)*dx(it)/dely(jt)
        Pn = Fn/Dn
        Ap_Pn = apfun(Pn)
        ant(it,jt) = Dn*Ap_Pn+max(-Fn,0.0)

       !atp_zero(it,jt)=rho(it,jt)*dx(it)*dy(jt)/(time_step*A22)
       !bt(it,jt) = atp_zero(it,jt)*(t0(it,jt)+f_temp(it,jt)*(time_step*A21))
        bt(it,jt) = 0.0d0

    end do

end do

! locally parabolic condition

!do jt=2,jmax-1

!    ate(imax-1,jt) = 0.0

!end do

do jt=2,jmax-1

    do it=2,imax-1

!        atp(it,jt) = atw(it,jt) + ate(it,jt) + atn(it,jt) + ats(it,jt) + atp_zero(it,jt)
         atp(it,jt) = awt(it,jt) + aet(it,jt) + ant(it,jt) + ast(it,jt)

    end do

end do

end subroutine
