subroutine face_velocity(u, v, p, uw, ue, vs, vn, pw, pe, ps, pn, imax, jmax)
use kind_parameters
implicit none

integer:: i, j, imax, jmax
real(kind=dp), dimension(imax, jmax):: u,v, p, uw, ue, vs, vn, pw, pe, ps, pn

uw(2,:) = u(1,:)
ue(2,:) = 0.5*(u(2,:)+u(3,:))
pw(2,:) = p(1,:)
pe(2,:) = 0.5*(p(2,:)+p(3,:))

do j = 1,jmax

 do i = 3,imax-1

  uw(i,j) = 0.5*(u(i-1,j) + u(i,j)) 
  ue(i,j) = 0.5*(u(i,j) + u(i+1,j))
  pw(i,j) = 0.5*(p(i-1,j) + p(i,j)) 
  pe(i,j) = 0.5*(p(i,j) + p(i+1,j))

 end do

end do

ue(imax-1,:) = u(imax,:)
pe(imax-1,:) = p(imax,:)

vs(:,2) = v(:,1)
vn(:,2) = 0.5*(v(:,2) + v(:,3)) 

ps(:,2) = p(:,1)
pn(:,2) = 0.5*(p(:,2) + p(:,3)) 

do j = 3,jmax-1

 do i = 1,imax

  vs(i,j) = 0.5*(v(i,j-1) + v(i,j)) 
  vn(i,j) = 0.5*(v(i,j) + v(i,j+1))
  ps(i,j) = 0.5*(p(i,j-1) + p(i,j)) 
  pn(i,j) = 0.5*(p(i,j) + p(i,j+1))

 end do

end do

vn(:,jmax-1) = v(:,jmax) 
pn(:,jmax-1) = p(:,jmax) 

end subroutine
