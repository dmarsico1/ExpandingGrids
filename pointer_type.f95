module pointer_type

type var_pointer
real*8, dimension(:,:,:), pointer :: p_var
end type var_pointer

type grid_pointer
real*8, dimension(:,:), pointer :: p_grid
end type grid_pointer

type single_var_pointer
real*8, dimension(:), pointer :: p_var
end type single_var_pointer

end module pointer_type

