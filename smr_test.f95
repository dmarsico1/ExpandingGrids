program smr_test

use pointer_type
use grid
use netcdf_data
use init
use step
use therm

implicit none

real*8 :: start,finish

call cpu_time(start)

call define_name_vars

call define_grid

call define_Bc

call init_grid

call calc_write_dim_freq_vec

call init_qvs

call calc_dtheta_e_vec

call init_sponge

call init_surf_flux

call init_var(.true.)

call init_netcdf_vars


do i=1,N_grids
	call write_netcdf_grids(x_grid(i)%p_grid(1,:),z_grid(i)%p_grid(1,:),x_grid(i)%p_grid(2,:),z_grid(i)%p_grid(2,:),i)
enddo

call time_stepper(.true.)

call init_var(.false.)

call time_stepper(.false.)

nc_rec = nf90_close(ncid)

call cpu_time(finish)

write(*,*) '("Time = ",f8.3," seconds.")',finish-start

end program smr_test

