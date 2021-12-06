module step

use grid
use interp
use Boundary
use average
use therm
use netcdf_data
use rk_time_step

integer :: rk=1
character(len=100) :: temp_str1,temp_str2,temp_str3,temp_str4
character(len=100) :: cwd
integer :: r,k,q,i,m
real*8, allocatable, dimension(:,:) :: U_T,W_T,Bu_T,Bs_T,P_T

real*8, allocatable, dimension(:,:) :: BC_u,BC_w,BC_bu,BC_bs
real*8, allocatable, dimension(:) :: Bc_p
real*8, allocatable, dimension(:) :: Bx1,PTBx,PTBx2

contains

subroutine time_stepper(start)

implicit none

logical :: start
real*8 :: time=0


q=r


i=1



m=1
if(start)then
	call time_step(var,dt/1000.0)
else
	do while(time<Tend)
		if((write_freq_vec(m)-write_freq_dt <= time) .and. (time <= write_freq_vec(m)+write_freq_dt))then
			call write_full_grid_vars(m)
			m=m+1
			!write(*,*) m
		endif
		write(*,*) "Time =", time
		call time_step(var,dt)
		call calc_ql
		time=time+dt
		rk=1
	enddo
endif




end subroutine time_stepper

subroutine time_step(var,dt)

implicit none

type(var_pointer), dimension(:) :: var
real*8 :: dt
integer, dimension(N_grids) :: step_vec,stop_vec
integer :: i,k

step_vec=0

do i=1,N_grids
	stop_vec(i)=rat**(i-1)
enddo

if(N_grids==1)then
	allocate(Bc_p(x_dim_Vec(1)),Ptbx2(x_dim_vec(1)))
	call Lowerboundary(var(1)%p_var(:,:,1),var(1)%p_var(:,:,2),var(1)%p_var(:,:,3),&
						var(1)%p_var(:,:,4),var(1)%p_var(NGz+1,NGx+1:NGx+x_dim_vec(1),1),&
						var(1)%p_var(NGz+1,NGx+1:NGx+x_dim_vec(1),3),&
						var(1)%p_var(NGz+1,NGx+1:NGx+x_dim_vec(1),4),&
						hx,hz,z_dim_vec(1),x_dim_vec(1))
	var(1)%p_var(2,:,2)=0.0
	var(1)%p_var(1,:,2)=0.0
	var(1)%p_var(2*NGz+z_dim_vec(1)-1,:,3)=var(1)%p_var(NGz+z_dim_vec(1),:,3)
	var(1)%p_var(2*NGz+z_dim_vec(1),:,3)=var(1)%p_var(NGz+z_dim_vec(1)-1,:,3)
	var(1)%p_var(2*NGz+z_dim_vec(1)-1,:,4)=var(1)%p_var(NGz+z_dim_vec(1),:,4)
	var(1)%p_var(2*NGz+z_dim_vec(1),:,4)=var(1)%p_var(NGz+z_dim_vec(1)-1,:,4)
	var(1)%p_var(2*NGz+z_dim_vec(1)-1,:,1)=var(1)%p_var(NGz+z_dim_vec(1),:,1)
	var(1)%p_var(2*NGz+z_dim_vec(1),:,1)=var(1)%p_var(NGz+z_dim_vec(1)-1,:,1)
	var(1)%p_var(2*NGz+z_dim_vec(1)-1,:,2)=0.0

	call PeriodicBoundary(var(1)%p_var(:,:,1),var(1)%p_var(:,:,2),var(1)%p_var(:,:,3),var(1)%p_var(:,:,4),&
						var(1)%p_var(2:z_dim_vec(1)+NGz+1,NGx:NGx+x_dim_vec(1)+1,5),z_dim_vec(1),x_dim_vec(1))
	call rk_stepper(var(1)%p_var(:,:,1),var(1)%p_var(1:z_dim_vec(1)+2*NGz-1,:,2),&
					var(1)%p_var(:,:,3),var(1)%p_var(:,:,4),var(1)%p_var(2:z_dim_vec(1)+NGz+1,NGx:NGx+x_dim_vec(1)+1,5),&
					z_dim_vec(1),x_dim_vec(1),dt,hx,hz,1,0.0*Bc_p/(hz/(rat**(1-1))),Ptbx2)
	deallocate(Bc_p,Ptbx2)
else

do while(step_vec(N_grids)<rat**(N_grids-1))
	call recursive_step(var,step_vec,stop_vec)
enddo

endif

end subroutine time_step

recursive subroutine recursive_step(var,step_vec,stop_vec)

implicit none

type(var_pointer), dimension(:) :: var
integer, dimension(N_grids),intent(out) :: step_vec
integer, dimension(N_grids), intent(in) :: stop_vec
integer :: i=1

if(i==N_grids)then
	allocate(BC_p(x_dim_vec(i)),ptbx2(x_dim_vec(i)))
	
	call Lowerboundary(var(i)%p_var(:,:,1),var(i)%p_var(:,:,2),var(i)%p_var(:,:,3),&
						var(i)%p_var(:,:,4),var(i)%p_var(NGz+1,NGx+1:NGx+x_dim_vec(i),1),&
						var(i)%p_var(NGz+1,NGx+1:NGx+x_dim_vec(i),3),&
						var(i)%p_var(NGz+1,NGx+1:NGx+x_dim_vec(i),4),&
						hx/(rat**(i-1)),hz/(rat**(i-1)),z_dim_vec(i),x_dim_vec(i))
	

		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,3),&
													var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,3),&
													var(i)%p_var(NGz+z_dim_vec(i)+1,NGx+1:x_dim_vec(i)+NGx,3),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
													z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(2,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-2)))

		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,3),&
													var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,3),&
													var(i)%p_var(2*NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,3),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
													z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(2,z_dim_vec(i))+2.0*hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-2)))
		
		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,4),&
													var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,4),&
													var(i)%p_var(NGz+z_dim_vec(i)+1,NGx+1:x_dim_vec(i)+NGx,4),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
													z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(2,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-2)))

		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,4),&
													var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,4),&
													var(i)%p_var(2*NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,4),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
													z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(2,z_dim_vec(i))+2.0*hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-2)))


	call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,2),&
													var(i)%p_var(NGz+z_dim_vec(i)-2:NGz+z_dim_vec(i)-1,NGx+1:x_dim_vec(i)+NGx,2),&
													var(i)%p_var(NGz+z_dim_vec(i)+1,NGx+1:x_dim_vec(i)+NGx,2),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(1,1),z_grid(i-1)%p_grid(1,2),&
													z_grid(i)%p_grid(1,z_dim_vec(i)-2),z_grid(i)%p_grid(1,z_dim_vec(i)-1),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(1,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-2)))



		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,1),&
											var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,1),&
											var(i)%p_var(NGz+z_dim_vec(i)+1,NGx+1:x_dim_vec(i)+NGx,1),&
											x_grid(i-1)%p_grid(1,:),x_grid(i)%p_grid(1,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
											z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(1,:),&
											z_grid(i)%p_grid(2,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz)

		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,1),&
											var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,1),&
											var(i)%p_var(NGz+z_dim_vec(i)+2,NGx+1:x_dim_vec(i)+NGx,1),&
											x_grid(i-1)%p_grid(1,:),x_grid(i)%p_grid(1,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
											z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(1,:),&
											z_grid(i)%p_grid(2,z_dim_vec(i))+2.0*hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz)
	


	call NeumannBoundary(Bc_p,0.5*(var(i-1)%p_var(3:6,NGx+1:x_dim_vec(i-1)+NGx,5)+BC_PT(i-1)%p_var(:,:,1)),&
						x_dim_vec(i-1),x_grid(i-1)%p_grid(2,:),&
						x_grid(i)%p_grid(2,:),hx/rat**(i-2),hz/(rat**(i-2)))

	

	call PeriodicBoundary(var(i)%p_var(:,:,1),var(i)%p_var(:,:,2),var(i)%p_var(:,:,3),var(i)%p_var(:,:,4),&
						var(i)%p_var(2:z_dim_vec(i)+NGz+1,NGx:NGx+x_dim_vec(i)+1,5),z_dim_vec(i),x_dim_vec(i))
	!call rk time step
	call rk_stepper(var(i)%p_var(:,:,1),var(i)%p_var(1:z_dim_vec(i)+2*NGz-1,:,2),&
					var(i)%p_var(:,:,3),var(i)%p_var(:,:,4),var(i)%p_var(2:z_dim_vec(i)+NGz+1,NGx:NGx+x_dim_vec(i)+1,5),&
					z_dim_vec(i),x_dim_vec(i),dt/(rat**(i-1)),hx/(rat**(i-1)),hz/(rat**(i-1)),i,Bc_p/(hz/(rat**(i-1))),Ptbx2)

	step_vec(N_grids)=step_vec(N_grids)+1
	

	Bc(i-1)%p_var(:,:,1)=0.5*(Bc(i-1)%p_var(:,:,1)+&
							var(i-1)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i-1),1))
	Bc(i-1)%p_var(:,:,2)=0.5*(Bc(i-1)%p_var(:,:,2)+&
							var(i-1)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i-1),2))
	Bc(i-1)%p_var(:,:,3)=0.5*(Bc(i-1)%p_var(:,:,3)+&
							var(i-1)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i-1),3))
	Bc(i-1)%p_var(:,:,4)=0.5*(Bc(i-1)%p_var(:,:,4)+&
							var(i-1)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i-1),4))
		
	call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,2),&
													var(i)%p_var(NGz+z_dim_vec(i)-2:NGz+z_dim_vec(i)-1,NGx+1:x_dim_vec(i)+NGx,2),&
													var(i)%p_var(NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,2),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(1,1),z_grid(i-1)%p_grid(1,2),&
													z_grid(i)%p_grid(1,z_dim_vec(i)-2),z_grid(i)%p_grid(1,z_dim_vec(i)-1),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(1,z_dim_vec(i)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-1)))




		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,3),&
													var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,3),&
													var(i)%p_var(NGz+z_dim_vec(i)+1,NGx+1:x_dim_vec(i)+NGx,3),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
													z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(2,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-2)))

		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,3),&
													var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,3),&
													var(i)%p_var(2*NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,3),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
													z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(2,z_dim_vec(i))+2.0*hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-2)))
		
		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,4),&
													var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,4),&
													var(i)%p_var(NGz+z_dim_vec(i)+1,NGx+1:x_dim_vec(i)+NGx,4),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
													z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(2,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-2)))

		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,4),&
													var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,4),&
													var(i)%p_var(2*NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,4),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
													z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(2,z_dim_vec(i))+2.0*hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-2)))

	call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,2),&
													var(i)%p_var(NGz+z_dim_vec(i)-2:NGz+z_dim_vec(i)-1,NGx+1:x_dim_vec(i)+NGx,2),&
													var(i)%p_var(NGz+z_dim_vec(i)+1,NGx+1:x_dim_vec(i)+NGx,2),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(1,1),z_grid(i-1)%p_grid(1,2),&
													z_grid(i)%p_grid(1,z_dim_vec(i)-2),z_grid(i)%p_grid(1,z_dim_vec(i)-1),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(1,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-2)))
	
		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,1),&
											var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,1),&
											var(i)%p_var(NGz+z_dim_vec(i)+1,NGx+1:x_dim_vec(i)+NGx,1),&
											x_grid(i-1)%p_grid(1,:),x_grid(i)%p_grid(1,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
											z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(1,:),&
											z_grid(i)%p_grid(2,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz)

		call CubicInterpT_coarse_fine(Bc(i-1)%p_var(:,:,1),&
											var(i)%p_var(NGz+z_dim_vec(i)-1:NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,1),&
											var(i)%p_var(NGz+z_dim_vec(i)+2,NGx+1:x_dim_vec(i)+NGx,1),&
											x_grid(i-1)%p_grid(1,:),x_grid(i)%p_grid(1,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
											z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(1,:),&
											z_grid(i)%p_grid(2,z_dim_vec(i))+2.0*hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz)

	



	
	call Lowerboundary(var(i)%p_var(:,:,1),var(i)%p_var(:,:,2),var(i)%p_var(:,:,3),&
						var(i)%p_var(:,:,4),var(i)%p_var(NGz+1,NGx+1:NGx+x_dim_vec(i),1),&
						var(i)%p_var(NGz+1,NGx+1:NGx+x_dim_vec(i),3),&
						var(i)%p_var(NGz+1,NGx+1:NGx+x_dim_vec(i),4),&
						hx/(rat**(i-1)),hz/(rat**(i-1)),z_dim_vec(i),x_dim_vec(i))

	
	call NeumannBoundary(Bc_p,var(i-1)%p_var(3:6,NGx+1:x_dim_vec(i-1)+NGx,5),&
						x_dim_vec(i-1),x_grid(i-1)%p_grid(2,:),&
						x_grid(i)%p_grid(2,:),hx/rat**(i-2),hz/(rat**(i-2)))


	call PeriodicBoundary(var(i)%p_var(:,:,1),var(i)%p_var(:,:,2),var(i)%p_var(:,:,3),var(i)%p_var(:,:,4),&
						var(i)%p_var(2:z_dim_vec(i)+NGz+1,NGx:NGx+x_dim_vec(i)+1,5),z_dim_vec(i),x_dim_vec(i))

	call rk_stepper(var(i)%p_var(:,:,1),var(i)%p_var(1:z_dim_vec(i)+2*NGz-1,:,2),&
					var(i)%p_var(:,:,3),var(i)%p_var(:,:,4),var(i)%p_var(2:z_dim_vec(i)+NGz+1,NGx:NGx+x_dim_vec(i)+1,5),&
					z_dim_vec(i),x_dim_vec(i),dt/(rat**(i-1)),hx/(rat**(i-1)),hz/(rat**(i-1)),i,Bc_p/(hz/(rat**(i-1))),&
					var(i-1)%p_var(NGz+z_dim_vec(i-1),NGx+1:NGx+x_dim_vec(i-1),2))

	step_vec(N_grids)=step_vec(N_grids)+1


	call CubicInterpT_coarse_fine(var(i-1)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i-1),2),&
													var(i)%p_var(NGz+z_dim_vec(i)-2:NGz+z_dim_vec(i)-1,NGx+1:x_dim_vec(i)+NGx,2),&
													var(i)%p_var(NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,2),&
													x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(1,1),z_grid(i-1)%p_grid(1,2),&
													z_grid(i)%p_grid(1,z_dim_vec(i)-2),z_grid(i)%p_grid(1,z_dim_vec(i)-1),x_grid(i)%p_grid(2,:),&
													z_grid(i)%p_grid(1,z_dim_vec(i)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-1)))

	
	call PeriodicBoundary(var(i)%p_var(:,:,1),var(i)%p_var(:,:,2),var(i)%p_var(:,:,3),var(i)%p_var(:,:,4),&
						var(i)%p_var(2:z_dim_vec(i)+NGz+1,NGx:NGx+x_dim_vec(i)+1,5),z_dim_vec(i),x_dim_vec(i))
	
	deallocate(Bc_p,ptbx2)

elseif(step_vec(i)<stop_vec(i) .and. (rat*step_vec(i) .eq. step_vec(i+1)) .and. i < N_grids)then

	Bc(i)%p_var(:,:,1)=var(i)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i),1)
	Bc(i)%p_var(:,:,2)=var(i)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i),2)
	Bc(i)%p_var(:,:,3)=var(i)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i),3)
	Bc(i)%p_var(:,:,4)=var(i)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i),4)
	Bc(i)%p_var(:,:,5)=var(i)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i),5)
	BC_PT(i)%p_var(:,:,1)=var(i)%p_var(NGz+1:NGz+4,NGx+1:NGx+x_dim_vec(i),5)

	if(i==1)then
		call UpperBoundary(var(i)%p_var(:,:,1),var(i)%p_var(:,:,2),var(i)%p_var(:,:,3),&
						var(i)%p_var(:,:,4),var(i)%p_var(NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),1),&
						var(i)%p_var(NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),3),&
						var(i)%p_var(NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),4),&
						hx/(rat**(i-1)),hz/(rat**(i-1)),z_dim_vec(i),x_dim_vec(i))

	endif

	allocate(U_T(z_dim_vec_ext(i)+2*NGz,x_dim_vec(i)+2*NGx),W_T(z_dim_vec_ext(i)+2*NGz-1,x_dim_vec(i)+2*NGx),&
			Bu_T(z_dim_vec_ext(i)+2*NGz,x_dim_vec(i)+2*NGx),Bs_T(z_dim_vec_ext(i)+2*NGz,x_dim_vec(i)+2*NGx),&
			P_T(z_dim_vec_ext(i)+2,x_dim_vec(i)+2))
	
	allocate(BC_u(2,x_dim_vec(i-1)),BC_w(2,x_dim_vec(i-1)),BC_bu(2,x_dim_vec(i-1)),BC_bs(2,x_dim_vec(i-1)),Bc_p(x_dim_vec(i)),&
			ptbx2(x_dim_vec(i)))


call Grid_average(U_T(NGz+1:NGz+z_dim_vec_ext(i),NGx+1:NGx+x_dim_vec(i)),&
					W_T(NGz+1:NGz+z_dim_vec_ext(i),NGx+1:NGx+x_dim_vec(i)),&
					Bu_T(NGz+1:NGz+z_dim_vec_ext(i),NGx+1:NGx+x_dim_vec(i)),&
					Bs_T(NGz+1:NGz+z_dim_vec_ext(i),NGx+1:NGx+x_dim_vec(i)),&
					P_T(1+1:1+z_dim_vec_ext(i),1+1:1+x_dim_vec(i)),z_dim_vec_ext(i),x_dim_vec(i),N_grids-i+1,i)


	
	if(i==1)then
		P_T(z_dim_vec_ext(i)+2,2:x_dim_vec(i)+1)=P_T(z_dim_vec_ext(i)+1,2:x_dim_vec(i)+1)
		P_T(1,2:x_dim_vec(i)+1)=P_T(2,2:x_dim_vec(i)+1)
		U_T(NGz+z_dim_vec_ext(i)+1,NGx+1:x_dim_vec(i)+NGx)=U_T(NGz+z_dim_vec_ext(i),NGx+1:x_dim_vec(i)+NGx)
		U_T(NGz+z_dim_vec_ext(i)+2,NGx+1:x_dim_vec(i)+NGx)=U_T(NGz+z_dim_vec_ext(i)-1,NGx+1:x_dim_vec(i)+NGx)
		U_T(1,NGx+1:NGx+x_dim_vec(i))=U_T(4,NGx+1:NGx+x_dim_vec(i))
		U_T(2,NGx+1:NGx+x_dim_vec(i))=U_T(3,NGx+1:NGx+x_dim_vec(i))
		Bu_T(NGz+z_dim_vec_ext(i)+1,NGx+1:x_dim_vec(i)+NGx)=Bu_T(NGz+z_dim_vec_ext(i),NGx+1:x_dim_vec(i)+NGx)
		Bu_T(NGz+z_dim_vec_ext(i)+2,NGx+1:x_dim_vec(i)+NGx)=Bu_T(NGz+z_dim_vec_ext(i)-1,NGx+1:x_dim_vec(i)+NGx)
		Bs_T(NGz+z_dim_vec_ext(i)+1,NGx+1:x_dim_vec(i)+NGx)=Bs_T(NGz+z_dim_vec_ext(i),NGx+1:x_dim_vec(i)+NGx)
		Bs_T(NGz+z_dim_vec_ext(i)+2,NGx+1:x_dim_vec(i)+NGx)=Bs_T(NGz+z_dim_vec_ext(i)-1,NGx+1:x_dim_vec(i)+NGx)
		Bu_T(2,NGx+1:NGx+x_dim_vec(i))=Bu_T(3,NGx+1:NGx+x_dim_vec(i))
		Bs_T(2,NGx+1:NGx+x_dim_vec(i))=Bs_T(3,NGx+1:NGx+x_dim_vec(i))
		Bu_T(1,NGx+1:NGx+x_dim_vec(i))=Bu_T(4,NGx+1:NGx+x_dim_vec(i))
		Bs_T(1,NGx+1:NGx+x_dim_vec(i))=Bs_T(4,NGx+1:NGx+x_dim_vec(i))	
		W_T(2,NGx+1:NGx+x_dim_vec(i))=0
		W_T(1,NGx+1:NGx+x_dim_vec(i))=-W_T(3,NGx+1:NGx+x_dim_vec(i))
		W_T(NGz+z_dim_vec_ext(i),NGx+1:NGx+x_dim_vec(i))=0
		W_T(NGz+z_dim_vec_ext(i)+1,NGx+1:NGx+x_dim_vec(i))=0
	endif

	if(i>1)then
		if(2.0*(step_vec(i-1)-1)==step_vec(i))then
			BC_u=Bc(i-1)%p_var(:,:,1)
			BC_w=Bc(i-1)%p_var(:,:,2)
			BC_bu=Bc(i-1)%p_var(:,:,3)
			BC_bs=Bc(i-1)%p_var(:,:,4)

			call NeumannBoundary(Bc_p,0.5*(var(i-1)%p_var(3:6,NGx+1:x_dim_vec(i-1)+NGx,5)+BC_PT(i-1)%p_var(:,:,1)),&
						x_dim_vec(i-1),x_grid(i-1)%p_grid(2,:),&
						x_grid(i)%p_grid(2,:),hx/rat**(i-2),hz/(rat**(i-2)))
			Bc_p=Bc_p/(hz/(rat**(i-1)))

		else
			BC_u=0.5*(Bc(i-1)%p_var(:,:,1)+var(i-1)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i-1),1))
			BC_w=0.5*(Bc(i-1)%p_var(:,:,2)+var(i-1)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i-1),2))
			BC_bu=0.5*(Bc(i-1)%p_var(:,:,3)+var(i-1)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i-1),3))
			BC_bs=0.5*(Bc(i-1)%p_var(:,:,4)+var(i-1)%p_var(NGz+1:NGz+2,NGx+1:NGx+x_dim_vec(i-1),4))

			call NeumannBoundary(Bc_p,var(i-1)%p_var(3:6,NGx+1:x_dim_vec(i-1)+NGx,5),&
						x_dim_vec(i-1),x_grid(i-1)%p_grid(2,:),&
						x_grid(i)%p_grid(2,:),hx/rat**(i-2),hz/(rat**(i-2)))
			Bc_p=Bc_p/(hz/(rat**(i-1)))

	
		endif

		call CubicInterpT_coarse_fine(Bc_bu,&
											Bu_T(NGz+z_dim_vec_ext(i)-1:NGz+z_dim_vec_ext(i),NGx+1:x_dim_vec(i)+NGx),&
											Bu_T(NGz+z_dim_vec_ext(i)+1,NGx+1:x_dim_vec(i)+NGx),&
											x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
											z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
											z_grid(i)%p_grid(2,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz)

		call CubicInterpT_coarse_fine(Bc_bu,&
											Bu_T(NGz+z_dim_vec_ext(i)-1:NGz+z_dim_vec_ext(i),NGx+1:x_dim_vec(i)+NGx),&
											Bu_T(NGz+z_dim_vec_ext(i)+2,NGx+1:x_dim_vec(i)+NGx),&
											x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
											z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
											z_grid(i)%p_grid(2,z_dim_vec(i))+2.0*hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz)

		call CubicInterpT_coarse_fine(Bc_bs,&
											Bs_T(NGz+z_dim_vec_ext(i)-1:NGz+z_dim_vec_ext(i),NGx+1:x_dim_vec(i)+NGx),&
											Bs_T(NGz+z_dim_vec_ext(i)+1,NGx+1:x_dim_vec(i)+NGx),&
											x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
											z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
											z_grid(i)%p_grid(2,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz)

		call CubicInterpT_coarse_fine(Bc_bs,&
											Bs_T(NGz+z_dim_vec_ext(i)-1:NGz+z_dim_vec_ext(i),NGx+1:x_dim_vec(i)+NGx),&
											Bs_T(NGz+z_dim_vec_ext(i)+2,NGx+1:x_dim_vec(i)+NGx),&
											x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
											z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(2,:),&
											z_grid(i)%p_grid(2,z_dim_vec(i))+2.0*hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz)

		call CubicInterpT_coarse_fine(Bc_u,&
											U_T(NGz+z_dim_vec_ext(i)-1:NGz+z_dim_vec_ext(i),NGx+1:x_dim_vec(i)+NGx),&
											U_T(NGz+z_dim_vec_ext(i)+1,NGx+1:x_dim_vec(i)+NGx),&
											x_grid(i-1)%p_grid(1,:),x_grid(i)%p_grid(1,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
											z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(1,:),&
											z_grid(i)%p_grid(2,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz)



		call CubicInterpT_coarse_fine(Bc_u,&
											U_T(NGz+z_dim_vec_ext(i)-1:NGz+z_dim_vec_ext(i),NGx+1:x_dim_vec(i)+NGx),&
											U_T(NGz+z_dim_vec_ext(i)+2,NGx+1:x_dim_vec(i)+NGx),&
											x_grid(i-1)%p_grid(1,:),x_grid(i)%p_grid(1,:),z_grid(i-1)%p_grid(2,1),z_grid(i-1)%p_grid(2,2),&
											z_grid(i)%p_grid(2,z_dim_vec(i)-1),z_grid(i)%p_grid(2,z_dim_vec(i)),x_grid(i)%p_grid(1,:),&
											z_grid(i)%p_grid(2,z_dim_vec(i))+2.0*hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz)
	


		call CubicInterpT_coarse_fine(Bc_w,&
											W_T(NGz+z_dim_vec_ext(i)-2:NGz+z_dim_vec_ext(i)-1,NGx+1:x_dim_vec(i)+NGx),&
											W_T(NGz+z_dim_vec_ext(i)+1,NGx+1:x_dim_vec(i)+NGx),&
											x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(1,1),z_grid(i-1)%p_grid(1,2),&
											z_grid(i)%p_grid(1,z_dim_vec(i)-2),z_grid(i)%p_grid(1,z_dim_vec(i)-1),x_grid(i)%p_grid(2,:),&
											z_grid(i)%p_grid(1,z_dim_vec(i))+hz/(rat**(i-1)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz)
	


	endif
		P_T(1,2:x_dim_vec(i)+1)=P_T(2,2:x_dim_vec(i)+1)
		U_T(1,NGx+1:NGx+x_dim_vec(i))=U_T(4,NGx+1:NGx+x_dim_vec(i))
		U_T(2,NGx+1:NGx+x_dim_vec(i))=U_T(3,NGx+1:NGx+x_dim_vec(i))
		Bu_T(2,NGx+1:NGx+x_dim_vec(i))=Bu_T(3,NGx+1:NGx+x_dim_vec(i))
		Bs_T(2,NGx+1:NGx+x_dim_vec(i))=Bs_T(3,NGx+1:NGx+x_dim_vec(i))
		Bu_T(1,NGx+1:NGx+x_dim_vec(i))=Bu_T(4,NGx+1:NGx+x_dim_vec(i))
		Bs_T(1,NGx+1:NGx+x_dim_vec(i))=Bs_T(4,NGx+1:NGx+x_dim_vec(i))
		W_T(2,NGx+1:NGx+x_dim_vec(i))=0
		W_T(1,NGx+1:NGx+x_dim_vec(i))=-W_T(3,NGx+1:NGx+x_dim_vec(i))

	call PeriodicBoundary(U_T,W_T,Bu_T,Bs_T,P_T,z_dim_vec_ext(i),x_dim_vec(i))


	call rk_stepper(U_T,W_T,Bu_T,Bs_T,P_T,z_dim_vec_ext(i),x_dim_vec(i),&
					dt/(rat**(i-1)),hx/(rat**(i-1)),hz/(rat**(i-1)),i,Bc_p,Bc_p)

	var(i)%p_var(NGz+1:NGz+z_dim_vec(i),:,1)=U_T(NGz+z_dim_vec_ext(i)-z_dim_vec(i)+1:NGz+z_dim_vec_ext(i),:)
	var(i)%p_var(NGz+1:NGz+z_dim_vec(i),:,2)=W_T(NGz+z_dim_vec_ext(i)-z_dim_vec(i)+1:NGz+z_dim_vec_ext(i),:)
	var(i)%p_var(NGz+1:NGz+z_dim_vec(i),:,3)=Bu_T(NGz+z_dim_vec_ext(i)-z_dim_vec(i)+1:NGz+z_dim_vec_ext(i),:)
	var(i)%p_var(NGz+1:NGz+z_dim_vec(i),:,4)=Bs_T(NGz+z_dim_vec_ext(i)-z_dim_vec(i)+1:NGz+z_dim_vec_ext(i),:)
	var(i)%p_var(NGz+1:NGz+z_dim_vec(i),2:NGx+x_dim_vec(i)+1,5)=P_T(1+z_dim_vec_ext(i)-z_dim_vec(i)+1:1+z_dim_vec_ext(i),:)
	
	
	
	step_vec(i)=step_vec(i)+1
	
	if(i .ne. 1)then
		call CubicInterpT_coarse_fine(Bc_w,&
									var(i)%p_var(NGz+z_dim_vec(i)-2:NGz+z_dim_vec(i)-1,NGx+1:x_dim_vec(i)+NGx,2),&
									var(i)%p_var(NGz+z_dim_vec(i),NGx+1:x_dim_vec(i)+NGx,2),&
									x_grid(i-1)%p_grid(2,:),x_grid(i)%p_grid(2,:),z_grid(i-1)%p_grid(1,1),z_grid(i-1)%p_grid(1,2),&
									z_grid(i)%p_grid(1,z_dim_vec(i)-2),z_grid(i)%p_grid(1,z_dim_vec(i)-1),x_grid(i)%p_grid(2,:),&
									z_grid(i)%p_grid(1,z_dim_vec(i)),z_dim_vec(i-1),x_dim_vec(i-1),hx/(rat**(i-1)),hz/(rat**(i-1)))

	endif
	deallocate(U_T,W_T,Bu_T,Bs_T,P_T,Bc_u,Bc_w,Bc_bu,Bc_bs,Bc_p,ptbx2)
	i=i+1
	call recursive_step(var,step_vec,stop_vec)
else
	i=i+1
	call recursive_step(var,step_vec,stop_vec)
endif

i=1 !reset i to 1 for use in the next recursive call

end subroutine


end module step

