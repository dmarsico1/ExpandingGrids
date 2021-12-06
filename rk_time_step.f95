module rk_time_step

use grid
use advect
use therm
use diffuse
use Boundary
use poisson
use matrix
use therm

contains

subroutine rk_stepper(U,W,Bu,Bs,P,M,N,dt,hx,hz,grid_num,fu,pxu)


implicit none

integer :: M,N,grid_num
real*8, dimension(M+2*NGz,N+2*NGx), intent(inout) :: U,Bu,Bs
real*8, dimension(M+2,N+2), intent(inout) :: P
real*8, dimension(M+2*NGz-1,N+2*NGx), intent(inout) :: W
real*8, dimension(N) :: fu,pxu
real*8, dimension(M,N) :: Px
real*8, dimension(N) :: fl
real*8, dimension(N) :: fu2
real*8, dimension(M-1,N) :: Pz
real*8 :: U1(M+2*NGz,N+2*NGx,3),W1(M+2*NGz-1,N+2*NGx,3),Bu1(M+2*NGz,N+2*NGx,3),Bs1(M+2*NGz,N+2*NGx,3)
real*8 :: dt,hx,hz,Ubar
integer :: i,j

do i=1,3
	U1(:,:,i)=U
	W1(:,:,i)=W
	Bu1(:,:,i)=Bu
	Bs1(:,:,i)=Bs
enddo


if(grid_num==1)then
	U1(NGz,NGx+1:NGx+N,1)=U1(NGz+1,NGx+1:NGx+N,1)
	U1(1,NGx+1:NGx+N,1)=U1(NGz+2,NGx+1:NGx+N,1)
endif



call vel_rk_step(U1(NGz+1:NGz+M,NGx+1:NGx+N,2),U1(:,:,1),U1(:,:,1),U1(:,:,1),&
				W1(NGz+1:NGz+M-1,NGx+1:NGx+N,2),W1(:,:,1),W1(:,:,1),W1(:,:,1),&
				Bu1(:,:,1),Bs1(:,:,1),P,M,N,dt,hx,hz,grid_num)
if(grid_num==1)then
W1(NGz+M,NGx+1:NGx+N,2)=0.0
else

do i=1,N
	if(Bs1(NGz+M,NGx+i,1) < qvs(grid_num)%p_var(M))then
		W1(NGz+M,NGx+i,2)=W1(NGz+M,NGx+i,1)+g*(0.5/th0)*dt*(Bu1(NGz+M+1,NGx+i,1)+Bu1(NGz+M,NGx+i,1))+&
								g*(0.5)*(Rvd-(Lv/(cp*th0)))*dt*(Bs1(NGz+M+1,NGx+i,1)+Bs1(NGz+M,NGx+i,1))
	else
		W1(NGz+M,NGx+i,2)=W1(NGz+M,NGx+i,1)+g*(0.5/th0)*dt*(Bu1(NGz+M+1,NGx+i,1)+Bu1(NGz+M,NGx+i,1))+&
								+g*dt*(Rvd-(Lv/(cp*th0))+1)*(qvs(grid_num)%p_var(M))&
								-g*(0.5)*dt*(Bs1(NGz+M+1,NGx+i,1)+Bs1(NGz+M,NGx+i,1))

	endif
enddo

do j=1,N
	Ubar=(1.0/4.0)*(U1(NGz+M,NGx+j,1)+U1(NGz+M,NGx+j+1,1)+U1(NGz+M+1,NGx+j,1)+U1(NGz+M+1,NGx+j+1,1))
	W1(NGz+M,NGx+j,2)=W1(NGz+M,NGx+j,2)-dt*Ubar*(W1(NGz+M,NGx+j+1,1)-W1(NGz+M,NGx+j-1,1))/(2.0*hx) &

						-dt*mu1*(1.0/hx**4)*(W1(NGz+M,NGx+j-2,1)-4.0*W1(NGz+M,NGx+j-1,1)+6.0*W1(NGz+M,NGx+j,1)&
						-4.0*W1(NGz+M,NGx+j+1,1)+W1(NGz+M,NGx+j+2,1))&
						-dt*W1(NGz+M,NGx+j,1)*(3.0*W1(NGz+M,NGx+j,1)-4.0*W1(NGz+M-1,NGx+j,1)+W1(NGz+M-2,NGx+j,1))/(2.0*hz)&
						+mu2*(dt/hz**2)*(W1(NGz+M+1,NGx+j,1)-2.0*W1(NGz+M,NGx+j,1)+W1(NGz+M-1,NGx+j,1))

enddo

endif



fl=0.0


fu2=0.0*fu2/hz

W1(NGz,NGx+1:NGx+N,2)=0.0
W1(1,NGx+1:NGx+N,2)=-W1(3,NGx+1:NGx+N,1)






call PeriodicBoundary(U1(:,:,2),W1(:,:,2),Bu1(:,:,2),Bs1(:,:,2),P,M,N)


if(grid_num==1)then
call PoissonSolve(U1(NGz+1:NGz+M,NGx+1:NGx+N+1,2),W1(NGz:NGz+M,NGx+1:NGx+N,2),U1(2:NGz+M+1,1:NGx+N+1,2),&
				W1(1:NGz+M+1,2:NGx+N+1,2),P,N,M,dt,dble(10.0**(-4)),hx,hz,grid_num,fu2,&
				fl)
else

call PoissonSolve(U1(NGz+1:NGz+M,NGx+1:NGx+N+1,2),W1(NGz:NGz+M,NGx+1:NGx+N,2),U1(2:NGz+M+1,1:NGx+N+1,2),&
				W1(1:NGz+M+1,2:NGx+N+1,2),P,N,M,dt,dble(10.0**(-4)),hx,hz,grid_num,fu,&
				fl)
endif

call GradP(Px,Pz,P(2:M+1,2:N+1),hx,hz,M,N)


U1(NGz+1:NGz+M,NGx+1:NGx+N,2)=U1(NGz+1:NGz+M,NGx+1:NGx+N,2)-dt*Px
W1(NGz+1:NGz+M-1,NGx+1:NGx+N,2)=W1(NGz+1:NGz+M-1,NGx+1:NGx+N,2)-dt*Pz

if(grid_num > 1)then
W1(NGz+M,NGx+1:NGx+N,2)=W1(NGz+M,NGx+1:NGx+N,2)-dt*hz*fu
endif



call B_rk_time_step(Bu1(NGz+1:NGz+M,NGx+1:NGx+N,2),Bu1(:,:,1),Bu1(:,:,1),Bu1(:,:,1),&
					Bs1(NGz+1:NGz+M,NGx+1:NGx+N,2),Bs1(:,:,1),Bs1(:,:,1),Bs1(:,:,1),&
					U1(:,:,2),W1(:,:,2),M,N,dt,hx,hz,grid_num)




call PeriodicBoundary(U1(:,:,2),W1(:,:,2),Bu1(:,:,2),Bs1(:,:,2),P,M,N)



U(NGz+1:NGz+M,NGx+1:NGx+N)=U1(NGz+1:NGz+M,NGx+1:NGx+N,2)
W(NGz+1:NGz+M,NGx+1:NGx+N)=W1(NGz+1:NGz+M,NGx+1:NGx+N,2)

Bu(NGz+1:NGz+M,NGx+1:NGx+N)=Bu1(NGz+1:NGz+M,NGx+1:NGx+N,2)
Bs(NGz+1:NGz+M,NGx+1:NGx+N)=Bs1(NGz+1:NGz+M,NGx+1:NGx+N,2)

end subroutine rk_stepper


subroutine vel_rk_step(U,U1,U2,U3,W,W1,W2,W3,Bu1,Bs1,P,M,N,dt,hx,hz,grid_num)

implicit none
real*8, dimension(M,N) :: U
real*8, dimension(M-1,N) :: W
real*8, dimension(M+2,N+2) :: P
real*8, dimension(M+4,N+4) :: U1,U2,U3
real*8, dimension(M+3,N+4) :: W1,W2,W3
real*8, dimension(M+4,N+4) :: Bu1,Bs1
real*8, dimension(M,N) :: Ua,Ud,Px,ls_fu,Uhd,U_spg,W_spg
real*8, dimension(M-1,N) :: f
real*8, dimension(M-1,N) :: Wa,Wd,Pz,Bavg,Whd
integer :: M,N
integer :: i,j,grid_num
real*8 :: hx,hz,dt

call advectVel(Ua,Wa,U2,W2,M,N,hx,hz)
call Normal_viscVel(Ud,U3(NGz:M+NGz+1,NGx+1:N+NGx),Wd,W3(NGz:M+Ngz,NGx+1:NGx+N),M,N,hz)
call force(f,Bu1(3:M+NGz,3:NGx+N),Bs1(3:M+NGz,3:NGx+N),qvs(grid_num)%p_var(:),M,N)
call ls_force_vel(ls_fu,U1(NGz+1:M+NGz,NGx+1:N+NGx),M,N)
call Horiz_HyperDiff(U1(NGz+1:NGz+M,:),W1(NGz+1:NGz+M-1,:),Uhd,Whd,M,N,hx)
call sponge_layer(U,U_spg,M,N,grid_num)
call sponge_layer(W,W_spg,M,N,grid_num)

!average B to w grid
do i=NGz+1,NGz+M-1
	do j=NGx+1,NGz+N
		Bavg(i-NGz,j-NGx)=g*(0.5/th0)*(Bu1(i+1,j)+Bu1(i,j))
	enddo
enddo

U=U1(NGz+1:NGz+M,NGx+1:NGx+N)-dt*Ua+dt*ls_fu+dt*Uhd+dt*Ud+dt*U_spg
W=W1(NGz+1:NGz+M-1,NGx+1:NGx+N)+dt*f-dt*Wa+dt*Whd+dt*Wd+dt*W_spg

end subroutine vel_rk_step

subroutine B_rk_time_step(Bu,Bu1,Bu2,Bu3,Bs,Bs1,Bs2,Bs3,U1,W1,M,N,dt,hx,hz,grid_num)

implicit none
real*8, dimension(M,N) :: Bu,Bs
real*8, dimension(M+4,N+4) :: Bu1,Bs1,Bu2,Bs2,Bu3,Bs3
real*8, dimension(M+4,N+4) :: U1
real*8, dimension(M+3,N+4) :: W1
real*8, dimension(M,N) :: Bua,Bud,Bsa,Bsd,Wavg,Uavg,ls_fbu,ls_fbs,Buhd,Bshd,Bu_spg,Bs_spg
real*8, dimension(M+4,N+4) :: ql
real*8, dimension(M,N) :: qlz
integer :: M,N,grid_num
real*8 :: hx,hz,dt
integer :: i,j

call advectB(Bua,Bsa,Bu2,Bs2,U1,W1,M,N,hx,hz)

call Normal_viscB(Bud,Bu3(Ngz:NGz+M+1,NGx:NGx+N+1),Bsd,Bs3(Ngz:NGz+M+1,NGx:NGx+N+1),M,N,hz)
call ls_force_therm(ls_fBu,ls_fBs,Bu1(NGz+1:NGz+M,NGx+1:NGx+N),Bs1(NGz+1:NGz+M,NGx+1:NGx+N),M,N,grid_num)
call Horiz_HyperDiffB(Bu1(NGz+1:NGz+M,:),Bs1(NGz+1:NGz+M,:),Buhd,Bshd,M,N,hx)
call sponge_layer(Bu,Bu_spg,M,N,grid_num)
call sponge_layer(Bs,Bs_spg,M,N,grid_num)

do i=NGz+1,M+NGz
	do j=NGx+1,N+NGx
		Wavg(i-2,j-2)=0.5*(W1(i-1,j)+W1(i,j))
	enddo
enddo

call calc_ql_grid_k(ql(3:M+2,3:N+2),Bs,qvs(grid_num)%p_var(:),M,N,grid_num)

Bu=Bu1(NGz+1:NGz+M,NGx+1:NGx+N)-dt*Bua+dt*Buhd+dt*Bud+dt*ls_fbu+dt*Bu_spg
Bs=Bs1(NGz+1:NGz+M,NGx+1:NGx+N)-dt*Bsa+dt*Bshd+dt*Bsd+dt*ls_fbs+dt*Bs_spg

do i=1,M
	do j=1,N
		Bu(i,j)=Bu(i,j)-dt*Wavg(i,j)*dtheta_e_vec(grid_num)%p_var(i)
		Bs(i,j)=Bs(i,j)-dt*Wavg(i,j)*dqt_vec(grid_num)%p_var(i)
	enddo
enddo

ql(1,3:N+2)=ql(4,3:N+2)
ql(2,3:N+2)=ql(3,3:N+2)
if(grid_num==1)then
ql(M+3,3:N+2)=ql(M+2,3:N+2)
ql(M+4,3:N+2)=ql(M+1,3:N+2)
else

do i=1,N
	ql(M+3,2+i)=max(Bs1(NGz+M+1,NGx+i)-qvs_func(z_grid_ext(grid_num)%p_var(M)+hz),0.0)
	ql(M+4,2+i)=max(Bs1(NGz+M+2,Ngx+i)-qvs_func(z_grid_ext(grid_num)%p_var(M)+2.0*hz),0.0)
enddo

endif

ql(:,1)=ql(:,N+1)
ql(:,2)=ql(:,N+2)
ql(:,N+3)=ql(:,3)
ql(:,N+4)=ql(:,4)

do i=3,M+2
	do j=3,N+2
		qlz(i-2,j-2)=(1.0/(12.0*hz))*(-ql(i+2,j)+8.0*(ql(i+1,j)-ql(i-1,j))+ql(i-2,j)) &
					+(1.0/(12.0*hz))*sign(dble(1.0),Wavg(i-2,j-2))*(ql(i+2,j)-4.0*ql(i+1,j)+6.0*ql(i,j)-4.0*ql(i-1,j)+ql(i-2,j))

	enddo
enddo


Bu=Bu-VT*dt*(Lv/cp)*qlz+dt*sf(grid_num)%p_var(:,:,1)
Bs=Bs+VT*dt*qlz


end subroutine B_rk_time_step


subroutine vel_rk_prelim_step(U,U1,U2,W,W1,W2,Bu1,Bs1,P,M,N,dt,hx,hz)

implicit none
real*8, dimension(M,N) :: U
real*8, dimension(M-1,N) :: W
real*8, dimension(M+2,N+2) :: P
real*8, dimension(M+4,N+4) :: U1,U2
real*8, dimension(M+3,N+4) :: W1,W2
real*8, dimension(:,:) :: Bu1,Bs1
real*8, dimension(M,N) :: Ua,Px,f
real*8, dimension(M-1,N) :: Wa,Pz,Bavg
integer :: M,N
integer :: i,j
real*8 :: hx,hz,dt

call advectVel(Ua,Wa,U2,W2,M,N,hx,hz)
call GradP(Px,Pz,P(2:M+1,2:N+1),hx,hz,M,N)


do i=1,M-1
	do j=1,N
		Bavg(i,j)=g*(0.5/th0)*(Bu1(i+1,j)+Bu1(i,j))+g*(0.5)*(Rvd-(Lv/(cp*th0)))*(Bs1(i+1,j)+Bs1(i,j))
	enddo
enddo

U=U1(NGz+1:NGz+M,NGx+1:NGx+N)-dt*Ua
W=W1(NGz+1:NGz+M-1,NGx+1:NGx+N)-dt*Wa+dt*Bavg

end subroutine vel_rk_prelim_step



end module rk_time_step


