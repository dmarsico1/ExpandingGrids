module netcdf_data

use netcdf
use grid

implicit none

integer, dimension(3) :: u_dimids,w_dimids,t_dimids
integer :: rec_dimid,ncid,nc_rec
character(len=100), dimension(:,:), allocatable :: grid_str
character(len=100), dimension(:,:), allocatable :: var_str
integer, dimension(:,:), allocatable :: varid_array
integer, dimension(:,:), allocatable :: dimid_array
integer, dimension(:,:), allocatable :: vardimid_array

contains

subroutine init_netcdf_vars

implicit none

character(len=100) :: str_form,temp_str
integer :: i,j

allocate(grid_str(4,N_grids),var_str(6,N_grids),varid_array(6,N_grids),dimid_array(4,N_grids),vardimid_array(4,N_grids))

if(N_grids<10)then
	str_form='(I1)'
elseif(N_grids>=10)then
	str_form='(I2)'
endif

do i=1,N_grids
	write(temp_str,str_form) i
	grid_str(1,i)="xm"//trim(temp_str)
	grid_str(2,i)="zm"//trim(temp_str)
	grid_str(3,i)="xt"//trim(temp_str)
	grid_str(4,i)="zt"//trim(temp_str)
enddo

do i=1,N_grids
	write(temp_str,str_form) i
	var_str(1,i)="u"//trim(temp_str)
	var_str(2,i)="w"//trim(temp_str)
	var_str(3,i)="bu"//trim(temp_str)
	var_str(4,i)="bs"//trim(temp_str)
	var_str(5,i)="p"//trim(temp_str)
	var_str(6,i)="ql"//trim(temp_str)
enddo

nc_rec = nf90_create(file_name,NF90_CLOBBER,ncid)
nc_rec = nf90_def_dim(ncid,"time",nf90_unlimited,rec_dimid)

do i=1,N_grids
	write(temp_str,str_form) i
	nc_rec = nf90_def_dim(ncid,"xm"//trim(temp_str),x_dim_vec(i),dimid_array(1,i))
	nc_rec = nf90_def_dim(ncid,"zm"//trim(temp_str),z_dim_vec(i),dimid_array(2,i))
	nc_rec = nf90_def_dim(ncid,"xt"//trim(temp_str),x_dim_vec(i),dimid_array(3,i))
	nc_rec = nf90_def_dim(ncid,"zt"//trim(temp_str),z_dim_vec(i),dimid_array(4,i))
enddo

do i=1,N_grids
	u_dimids=(/dimid_array(4,i),dimid_array(1,i),rec_dimid/)
	w_dimids=(/dimid_array(2,i),dimid_array(3,i),rec_dimid/)
	t_dimids=(/dimid_array(4,i),dimid_array(3,i),rec_dimid/)
	nc_rec = nf90_def_var(ncid,var_str(1,i),nf90_double,u_dimids,varid_array(1,i))
	nc_rec = nf90_def_var(ncid,var_str(2,i),nf90_double,w_dimids,varid_array(2,i))
	nc_rec = nf90_def_var(ncid,var_str(3,i),nf90_double,t_dimids,varid_array(3,i))
	nc_rec = nf90_def_var(ncid,var_str(4,i),nf90_double,t_dimids,varid_array(4,i))
	nc_rec = nf90_def_var(ncid,var_str(5,i),nf90_double,t_dimids,varid_array(5,i))
	nc_rec = nf90_def_var(ncid,var_str(6,i),nf90_double,t_dimids,varid_array(6,i))
	nc_rec = nf90_def_var(ncid,grid_str(1,i),nf90_double,dimid_array(1,i),vardimid_array(1,i))
	nc_rec = nf90_def_var(ncid,grid_str(2,i),nf90_double,dimid_array(2,i),vardimid_array(2,i))
	nc_rec = nf90_def_var(ncid,grid_str(3,i),nf90_double,dimid_array(3,i),vardimid_array(3,i))
	nc_rec = nf90_def_var(ncid,grid_str(4,i),nf90_double,dimid_array(4,i),vardimid_array(4,i))
enddo

nc_rec = nf90_enddef(ncid)

end subroutine init_netcdf_vars

subroutine write_netcdf_vars(U,W,Bu,Bs,P,ql,k,i)

implicit none

integer :: start(3), ucounts(3), wcounts(3),tcounts(3),nc_rec,k,varid,i
real*8 :: time
real*8, dimension(:,:) :: U,W,Bu,Bs,P,ql
ucounts = (/z_dim_vec(i),x_dim_vec(i),1/)
wcounts = (/z_dim_vec(i),x_dim_vec(i),1/)
tcounts = (/z_dim_vec(i),x_dim_vec(i),1/)
start = (/1,1,k/)

nc_rec= nf90_put_var(ncid,varid_array(1,i),U,start,ucounts)
nc_rec= nf90_put_var(ncid,varid_array(2,i),W,start,wcounts)
nc_rec = nf90_put_var(ncid,varid_array(3,i),Bu,start,tcounts)
nc_rec = nf90_put_var(ncid,varid_array(4,i),Bs,start,tcounts)
nc_rec = nf90_put_var(ncid,varid_array(5,i),P,start,tcounts)
nc_rec = nf90_put_var(ncid,varid_array(6,i),ql,start,tcounts)


end subroutine write_netcdf_vars

subroutine write_full_grid_vars(k)

implicit none

integer :: k,i

do i=1,N_grids
		call write_netcdf_vars(var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),1),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),2),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),3),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),4),&
						var(i)%p_var(NGz+1:NGz+z_dim_vec(i),NGx+1:NGx+x_dim_vec(i),5),&
						ql(i)%p_var(:,:,1),k,i)
enddo

end subroutine write_full_grid_vars

subroutine write_netcdf_grids(xm,zm,xt,zt,i)

implicit none
real*8, dimension(:) :: xm,xt,zm,zt
integer :: i

nc_rec = nf90_put_var(ncid,vardimid_array(1,i),xm,(/1/),(/x_dim_vec(i)/))
nc_rec = nf90_put_var(ncid,vardimid_array(2,i),zm,(/1/),(/z_dim_vec(i)/))
nc_rec = nf90_put_var(ncid,vardimid_array(3,i),xt,(/1/),(/x_dim_vec(i)/))
nc_rec = nf90_put_var(ncid,vardimid_array(4,i),zt,(/1/),(/z_dim_vec(i)/))

end subroutine

end module netcdf_data


