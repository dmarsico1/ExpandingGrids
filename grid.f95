module grid

use pointer_type

implicit none

integer :: N_grids=4 !number of grids
type(var_pointer), dimension(:), allocatable :: var
type(var_pointer), dimension(:), allocatable :: ql
type(var_pointer), dimension(:), allocatable :: spge,sf
type(grid_pointer), dimension(:), allocatable :: x_grid,z_grid
type(var_pointer), dimension(:), allocatable :: Bc
type(var_pointer), dimension(:), allocatable :: BC_PT
type(single_var_pointer), dimension(:), allocatable :: qvs
type(single_var_pointer), dimension(:), allocatable :: qt_vec,dqt_vec,dtheta_vec,dtheta_e_vec
type(single_var_pointer), dimension(:), allocatable :: z_grid_ext
real*8, dimension(:), allocatable :: write_freq_vec
integer, dimension(:), allocatable :: x_dim_vec, z_dim_vec,z_dim_vec_ext !dimensions of the grids
integer :: N_vert=50 !number of coarsest vertical momentum points.
integer, parameter :: nls=100 !max number of vertical levels
integer, dimension(nls) :: ref_levels !vertical levels at which the grid is refined
integer :: rat=2 !grid refinement ratio (always 2)
integer :: N_x=64 !number of coarsest grid horizontal cells
integer :: N_scalars=5 !number of variables, including pressure
real*8 :: hz=200 !vertical grid spacing
real*8 :: hx=200 !horizontal grid spacing
integer :: NGx=2 !number of ghost cells in horizontal direction
integer :: NGz=2 !number of ghost cells in vertical direction
character(len=305) :: file_name='test'
real*8 :: Tend=200 !stopping time
real*8 :: dt=0.5 !coarsest grid time step
real*8 :: eb=0.0!1 !buoyancy diffusion coefficient
real*8 :: ev=0.0!1 !velocity diffusion coefficient
real*8 :: mu1=20.0
real*8 :: mu2=1.0 !horizontal viscosity coefficient
!real*8 :: mu1=5.0
!real*8 :: mu2=0.25 !horizontal viscosity coefficient
real*8 :: dtheta=4.75/1000.0 !Note that for grid14, this was 2/1000
real*8 :: dthetae
real*8 :: dqt=(-20.0/15000.0)*10.0**(-3)  !background total water gradient
real*8 :: qt0=22.0*10.0**(-3)
real*8 :: VT=0.0!2.0
integer :: qt_tilde_type=1
real*8 :: th0=300.0 !constant background potential temperature
real*8 :: max_moist=10.0*10.0**(-3) !for grid13 this was 25
real*8 :: max_moist_z=0.05 !was originally 0.09
real*8 :: max_t=-10.0!-12.0 !originally -2 !I think for grid13 this was -20
real*8 :: max_t_z=0.05!was originally 0.1
real*8 :: spge_h=0.3
real*8 :: BL_top=0.0 !Height of top of lower boundary layer
real*8 :: cp=1000
real*8 :: Rvd=0.61
real*8 :: Rv=461.5
real*8 :: Rd=287.1
real*8 :: Lv=2500.0*10.0**3
real*8 :: g=9.8
real*8 :: relax_u=86400.0 !relaxation parameter for horizontal velocity
real*8 :: therm_relax=86400.0 
real*8 :: write_freq=1.0
real*8 :: write_freq_dt=0.000001

contains

subroutine define_grid

implicit none

integer :: i

allocate(var(N_grids),x_grid(N_grids),z_grid(N_grids),x_dim_vec(N_grids),z_dim_vec(N_grids),z_dim_vec_ext(N_grids),qvs(N_grids),&
		BC_PT(N_grids),z_grid_ext(N_grids),ql(N_grids),dqt_vec(N_grids),dtheta_vec(N_grids),dtheta_e_vec(N_grids),qt_vec(N_grids),&
		spge(N_grids),sf(N_grids))

call grid_dims

do i=1,N_grids
	allocate(var(i)%p_var(z_dim_vec(i)+2*NGz,x_dim_vec(i)+2*NGx,N_scalars))  !It be 2*NGz and 2*NGx
	allocate(spge(i)%p_var(z_dim_vec_ext(i),x_dim_vec(i),1))
	allocate(sf(i)%p_var(z_dim_vec_ext(i),x_dim_vec(i),1))
	allocate(qvs(i)%p_var(z_dim_vec_ext(i)))
	allocate(dqt_vec(i)%p_var(z_dim_vec_ext(i)),dtheta_vec(i)%p_var(z_dim_vec_ext(i)),dtheta_e_vec(i)%p_var(z_dim_vec_ext(i)),&
			qt_vec(i)%p_var(z_dim_vec_ext(i)))
	allocate(x_grid(i)%p_grid(2,x_dim_vec(i)))
	allocate(z_grid(i)%p_grid(2,z_dim_vec(i)))
	allocate(BC_PT(i)%p_var(4,x_dim_vec(i),1))
	allocate(z_grid_ext(i)%p_var(z_dim_vec_ext(i)))
	allocate(ql(i)%p_var(z_dim_vec(i),x_dim_vec(i),1))
enddo

end subroutine define_grid

subroutine define_Bc

implicit none

integer :: i

allocate(Bc(N_grids))

do i=1,N_grids
	allocate(Bc(i)%p_var(NGz,x_dim_vec(i),N_scalars))
enddo

end subroutine define_Bc

subroutine grid_dims

implicit none

integer ::i

z_dim_vec(1)=ref_levels(1)-ref_levels(2)

do i=2,N_grids-1
	z_dim_vec(i)=2.0*ref_levels(i)-ref_levels(i+1)
enddo

if(N_grids>1)then
z_dim_vec(N_grids)=2.0*ref_levels(N_grids)
endif

do i=1,N_grids
	x_dim_vec(i)=N_x*(rat**(i-1))
enddo

z_dim_vec_ext(1)=ref_levels(1)

do i=2,N_grids-1
	z_dim_vec_ext(i)=2.0*ref_levels(i)
enddo

z_dim_vec_ext(N_grids)=z_dim_vec(N_grids)

end subroutine grid_dims

subroutine init_grid

implicit none

integer :: i

do i=1,N_grids
	call init_grid_k(i)
enddo

end subroutine init_grid

subroutine init_grid_k(k)

implicit none

integer, intent(in) :: k
integer :: i
real*8 :: dz_k,dx_k
real*8 :: vert_shift

dx_k=hx/(rat**(k-1))
dz_k=hz/(rat**(k-1))

if(k<N_grids)then
	vert_shift=dz_k*ref_levels(k+1)
else
	vert_shift=0
endif

do i=1,z_dim_vec(k)
	z_grid(k)%p_grid(1,i)=vert_shift+dz_k*dble(i)!dz_k*(i+N_vert-z_dim_vec(i)-1)
	z_grid(k)%p_grid(2,i)=vert_shift+dz_k*(dble(i)-0.5)!dz_k*(i+N_vert-z_dim_vec(i)-0.5)
enddo

do i=1,x_dim_vec(k)
	x_grid(k)%p_grid(1,i)=(real(i)-1.0)*dx_k
	x_grid(k)%p_grid(2,i)=(real(i)-0.5)*dx_k
enddo

do i=1,z_dim_vec_ext(k)
	z_grid_ext(k)%p_var(i)=dz_k*(dble(i)-0.5)
enddo

end subroutine init_grid_k

subroutine calc_write_dim_freq_vec

implicit none

integer :: dim_freq,i

dim_freq=dint(Tend/write_freq)

write(*,*) dim_freq

allocate(write_freq_vec(dim_freq))

do i=1,dim_freq
	write_freq_vec(i)=(i-1.0)*write_freq
enddo

end subroutine calc_write_dim_freq_vec

subroutine define_name_vars

implicit none

namelist /model/ &
	N_grids,N_vert,ref_levels,N_scalars,N_x,hx,hz,rat,NGx,NGz,Tend,dt,write_freq,eb,ev,file_name
open(unit=1,file='NAMELIST')
read(1,nml=model)
 close(1)

end subroutine define_name_vars


end module grid

