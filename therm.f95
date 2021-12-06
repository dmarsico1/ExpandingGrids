module therm

use grid
use pointer_type

type(single_var_pointer), dimension(:), allocatable :: qvs_tot

contains

subroutine init_qvs

implicit none

integer :: i

allocate(qvs_tot(N_grids))

do i=1,N_grids
	allocate(qvs_tot(i)%p_var(z_dim_vec_ext(i)))
enddo


do i=1,N_grids
	call init_qvs_grid_k(i,1)
enddo

call calc_qt_vec
call calc_dqt_vec

do i=1,N_grids
	call init_qvs_grid_k(i,2)
enddo

end subroutine init_qvs

subroutine init_qvs_grid_k(i,m)

implicit none

integer :: i,j,m

	if(m==1)then
	call calc_qvs(qvs_tot(i)%p_var(:),z_grid_ext(i)%p_var(:),z_dim_vec_ext(i))


	else
	call calc_qvs_pert(qvs(i)%p_var(:),qvs_tot(i)%p_var(:),z_grid_ext(i)%p_var(:),z_dim_vec_ext(i),i)
	endif

end subroutine init_qvs_grid_k


subroutine calc_qvs(qvs,zt,M)

!calculate the total saturation mixing ratio

implicit none

real*8 :: p0=10**5
real*8 :: e_s0=3500.0
integer, intent(in) :: M
real*8, dimension(M) :: qvs,zt
real*8, dimension(M) :: theta_tilde,p_tilde,T_tilde
integer :: i

do i=1,M
	if(zt(i)<BL_top)then
		theta_tilde(i)=th0
	else
		theta_tilde(i)=zt(i)*dtheta+th0
	endif
enddo

do i=1,M
	p_tilde(i)=p0*(1-(g/(cp*dtheta))*log(1+(dtheta/th0)*(zt(i))))**(cp/Rd)
enddo

do i=1,M
	T_tilde(i)=theta_tilde(i)*(p_tilde(i)/p0)**(Rd/cp)
enddo

do i=1,M
	qvs(i)=(Rd/Rv)*(e_s0/p_tilde(i))*exp(-(Lv/Rv)*((1.0/T_tilde(i))-(1.0/th0)))
enddo


end subroutine calc_qvs

subroutine calc_qvs_pert(qvs_pert,qvs_tot,zt,M,k)

!calculate the perturbation saturation mixing ratio by subtracting off the backround total water

implicit none

integer, intent(in) :: M
real*8, dimension(M) :: qvs_pert,qvs_tot,zt
integer :: i,k

	do i=1,M
		qvs_pert(i)=qvs_tot(i)-qt_vec(k)%p_var(i)
	enddo

end subroutine calc_qvs_pert

function VT_func(z)

implicit none

real*8, intent(in) :: z
real*8 :: rho_0,rho
real*8 :: gam=-0.3654
real*8 :: f0=14.34
real*8 :: p0=10.0**5
real*8 :: p_tilde,T_tilde,theta_tilde
integer :: i
real*8 :: VT_func

rho_0=p0/(Rd*th0)

if(z<BL_top)then
	theta_tilde=th0
else
	theta_tilde=z*dtheta+th0
endif

p_tilde=p0*(1-(g/(cp*dtheta))*log(1+(dtheta/th0)*(z)))**(cp/Rd)

T_tilde=theta_tilde*(p_tilde/p0)**(Rd/cp)

rho=p_tilde/(Rd*T_tilde)

VT_func=f0*sqrt(rho_0)*rho**gam

end function VT_func


function qvs_func(z)

implicit none

real*8 :: z
real*8 :: p0=10**5
real*8 :: e_s0=3500.0
real*8 :: theta_tilde,p_tilde,T_tilde,qvs,qvs_pert
real*8 :: qvs_func

if(z<BL_top)then
	theta_tilde=th0
else
	theta_tilde=z*dtheta+th0
endif


p_tilde=p0*(1-(g/(cp*dtheta))*log(1+(dtheta/th0)*(z)))**(cp/Rd)

T_tilde=theta_tilde*(p_tilde/p0)**(Rd/cp)

qvs=(Rd/Rv)*(e_s0/p_tilde)*exp(-(Lv/Rv)*((1.0/T_tilde)-(1.0/th0)))

qvs_func=qvs-0.8*qvs


end function

subroutine init_sponge

implicit none

integer :: k

do k=1,N_grids
	call sponge_grid_k(z_dim_vec_ext(k),x_dim_vec(k),k)
enddo

end subroutine init_sponge

subroutine sponge_grid_k(M,N,grid_num)

implicit none

integer :: i,j,grid_num,M,N

do i=1,M
	if(z_grid_ext(grid_num)%p_var(i)<N_vert*hz*spge_h)then
		spge(grid_num)%p_var(i,:,1)=0.0
	else
	do j=1,N
		spge(grid_num)%p_var(i,j,1)=(1.0/1800.0)*(sin((0.5*4.0*atan(1.0)/(N_vert*hz-spge_h*N_vert*hz))*(z_grid_ext(grid_num)%p_var(i)&
										-N_vert*hz*spge_h))**2)
	enddo
	endif
enddo


end subroutine sponge_grid_k

subroutine sponge_layer(u,t,M,N,grid_num)

implicit none

real*8, dimension(M,N), intent(out) :: t
real*8, dimension(M,N), intent(in) :: u
integer :: i,j,M,N,grid_num

do i=1,M
	do j=1,N
		t(i,j)=-spge(grid_num)%p_var(i,j,1)*u(i,j)
	enddo
enddo

end subroutine sponge_layer

subroutine force(f,Bu,Bs,qvs,M,N)

implicit none

real*8, intent(out), dimension(M-1,N) :: f
real*8, dimension(M,N) :: Bu,Bs
real*8, dimension(M) :: qvs
integer, intent(in) :: M,N
integer :: i,j
real*8 :: row_sum

do i=1,M-1
	do j=1,N
		if(Bs(i,j)<qvs(i))then
			f(i,j)=0.5*g*((1.0/th0)*(Bu(i,j)+Bu(i+1,j))+(Rvd-(Lv/(cp*th0)))*(Bs(i,j)+Bs(i+1,j)))
		else
			f(i,j)=0.5*g*((1.0/th0)*(Bu(i,j)+Bu(i+1,j))+(Rvd-(Lv/(cp*th0))+1)*(qvs(i)+qvs(i+1))-(Bs(i,j)+Bs(i+1,j)))
		endif
	enddo
enddo

do i=1,M-1
	row_sum=sum(f(i,:))/N
		do j=1,N
			f(i,j)=f(i,j)-row_sum
		enddo
enddo



end subroutine force

subroutine calc_dthetae

!calculate derivative of background equivalent potential temperature

implicit none

dthetae=dtheta+(Lv/cp)*dqt

end subroutine calc_dthetae


subroutine calc_dtheta_e_vec

implicit none

integer :: i,j

do i=1,N_grids
	do j=1,z_dim_vec_ext(i)
		if(z_grid_ext(i)%p_var(j)<=BL_top)then
			dtheta_e_vec(i)%p_var(j)=(Lv/cp)*dqt_vec(i)%p_var(j)
		else
			dtheta_e_vec(i)%p_var(j)=dtheta+(Lv/cp)*dqt_vec(i)%p_var(j)
		endif
	enddo
enddo

end subroutine calc_dtheta_e_vec


subroutine calc_qt_vec

implicit none

integer :: i,j,r

if(qt_tilde_type==1)then
do i=1,N_grids
	do j=1,z_dim_vec_ext(i)
			qt_vec(i)%p_var(j)=qt0*exp(3.0*(dqt/qt0)*(1.0/nls*hz)*(z_grid_ext(i)%p_var(j)-BL_top))
			qt_vec(i)%p_var(j)=0.87*qvs_tot(i)%p_var(j)
	enddo
enddo


elseif(qt_tilde_type==0)then
	do i=1,N_grids
		do j=1,z_dim_vec_ext(i)
			if(z_grid_ext(i)%p_var(j)<=BL_top)then
				qt_vec(i)%p_var(j)=qt0
			else
				qt_vec(i)%p_var(j)=qt0+dqt*(z_grid_ext(i)%p_var(j)-BL_top)
			endif
		enddo
	enddo
endif


end subroutine calc_qt_vec


subroutine calc_dqt_vec

implicit none

integer :: i,j,r

if(qt_tilde_type==1)then
do i=1,N_grids
	do j=1,z_dim_vec_ext(i)

			if(j==1)then
				dqt_vec(i)%p_var(j)=(qt_vec(i)%p_var(j+1)-qt_vec(i)%p_var(j))/(hz/(2.0**(i-1)))
			elseif(j==z_dim_vec_ext(i))then
				dqt_vec(i)%p_var(j)=(3.0*qt_vec(i)%p_var(j)-4.0*qt_vec(i)%p_var(j-1)+qt_vec(i)%p_var(j-2))/(2.0*hz/(2.0**(i-1)))
			else
				dqt_vec(i)%p_var(j)=(qt_vec(i)%p_var(j+1)-qt_vec(i)%p_var(j-1))/(2.0*hz/(2.0**(i-1)))
			endif
	enddo
enddo

elseif(qt_tilde_type==0)then
	do i=1,N_grids
		do j=1,z_dim_vec_ext(i)
			if(z_grid_ext(i)%p_var(j)<=BL_top)then
				dqt_vec(i)%p_var(j)=0.0
			else
				dqt_vec(i)%p_var(j)=dqt
			endif
		enddo
	enddo
endif

end subroutine calc_dqt_vec

subroutine calc_ql

implicit none

integer :: i


do i=1,N_grids
	call calc_ql_grid_k(ql(i)%p_var(:,:,1),var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),4),&
						qvs(i)%p_var(z_dim_vec_ext(i)-z_dim_vec(i)+1:z_dim_vec_ext(i)),z_dim_vec(i),x_dim_vec(i),i)
enddo

end subroutine calc_ql

subroutine init_surf_flux

integer :: k

do k=1,N_grids
	call surf_flux(k,sf(k)%p_var(:,:,1),z_dim_vec_ext(k),x_dim_vec(k))
enddo

end subroutine init_surf_flux

subroutine surf_flux(grid_num,sf,M,N)

implicit none

integer, intent(in) :: grid_num,M,N
real*8, intent(out), dimension(M,N) :: sf
integer :: i,j

do i=1,M
	do j=1,N
		if(z_grid_ext(grid_num)%p_var(i)<301)then
			sf(i,j)=(1.0/1800.0)*(sin((0.5*4.0*atan(1.0)/(301))*(z_grid_ext(grid_num)%p_var(i)&
										-301))**2)
		else
			sf(i,j)=0.0
		endif
	enddo
enddo

end subroutine surf_flux

subroutine calc_ql_grid_k(ql,qt,qvs,M,N,k)

implicit none

real*8, dimension(M,N) :: qt,ql
real*8, dimension(M) :: qvs
integer :: M,N
integer :: k
integer :: i,j

do i=1,M
	do j=1,N
		ql(i,j)=max(qt(i,j)-qvs(i),0.0)
	enddo
enddo

end subroutine calc_ql_grid_k


subroutine ls_force_vel(ls_fu,U,M,N)

implicit none

integer :: i,M,N
real*8, dimension(M,N) :: ls_fu,U

do i=1,M
	ls_fu(i,:)=-(sum(U(i,:))/N)/relax_u
enddo

end subroutine ls_force_vel

subroutine ls_force_therm(ls_fBu,ls_fBs,Bu,Bs,M,N,grid_num)

implicit none

integer :: i,M,N,j,grid_num
real*8, dimension(M,N) :: ls_fBu,ls_fBs,Bu,Bs


do i=1,M
	do j=1,N
		ls_fbs(i,j)=max_moist*(1.0/therm_relax)*(1.0/(max_moist_z*N_vert*hz))*z_grid_ext(grid_num)%p_var(i)*&
					exp(-(1.0/2.0)*((1.0/(max_moist_z*N_vert*hz))**2)*z_grid_ext(grid_num)%p_var(i)**2)
	enddo
enddo


do i=1,M
	do j=1,N
		ls_fbu(i,j)=max_t*(1.0/therm_relax)*(1.0/(max_t_z*N_vert*hz))*z_grid_ext(grid_num)%p_var(i)*&
					exp(-(1.0/2.0)*((1.0/(max_t_z*N_vert*hz))**2)*z_grid_ext(grid_num)%p_var(i)**2)&
					+(Lv/cp)*ls_fbs(i,j)
	enddo
enddo



end subroutine ls_force_therm

end module therm

