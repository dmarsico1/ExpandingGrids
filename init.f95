module init

use grid
use therm

contains

subroutine init_var(first_step)

implicit none

integer :: i
logical :: first_step
if(first_step)then
	do i=1,N_grids
		call init_var_k(i,first_step,var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),1),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),2),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),3),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),4),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),5))
	enddo
else
	do i=1,N_grids
		call init_var_k(i,first_step,var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),1),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),2),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),3),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),4),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),5))
	enddo
endif

end subroutine init_var

subroutine init_var_k(k,init,U1,U2,U3,U4,U5)

implicit none

integer, intent(in) :: k
integer ::i,j
real*8, dimension(x_dim_vec(k)) :: xm,xt
real*8, dimension(z_dim_vec(k)) :: zm,zt
logical :: init
real*8, dimension(:,:), intent(inout) :: U1,U2,U3,U4,U5

xm=x_grid(k)%p_grid(1,:)
xt=x_grid(k)%p_grid(2,:)
zm=z_grid(k)%p_grid(1,:)
zt=z_grid(k)%p_grid(2,:)

if(init)then
	do i=1,z_dim_vec(k)
		do j=1,x_dim_vec(k)
			U1(i,j)=0.0
			U2(i,j)=0.0
			if(zt(i)<301)then
				call random_number(U3(i,j))
				U3(i,j)=0.2*(2.0*U3(i,j)-1.0)
			else
				U3(i,j)=0.0
			endif
			U4(i,j)=0.0
			U5(i,j)=0.0
		enddo
	enddo
else
	do i=1,z_dim_vec(k)
		do j=1,x_dim_vec(k)
			U1(i,j)=0.0
			U2(i,j)=0.0
			if(zt(i)<301)then
				call random_number(U3(i,j))
				U3(i,j)=0.2*(2.0*U3(i,j)-1.0)
			else
				U3(i,j)=0.0
			endif
			U4(i,j)=0.0
		enddo
	enddo
endif

end subroutine init_var_k

end module init
