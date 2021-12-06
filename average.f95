module average

use grid
use interp
use pointer_type

contains



subroutine TGrid_average(Ng1,Ng2,Tc,Tf,M,N)

implicit none

real*8, dimension(M,N) :: Tc
real*8, dimension(M*rat**(Ng2-Ng1),N*rat**(Ng2-Ng1)) :: Tf
integer :: Ng1,Ng2
integer :: M,N
integer :: i,j
integer :: diff

diff=Ng2-Ng1

do i=1,M
	do j=1,N
		Tc(i,j)=(1.0/((rat**2)**diff))*(sum(Tf((i-1)*(rat**diff)+1:i*(rat**diff),(j-1)*(rat**diff)+1:j*(rat**diff))))
	enddo
enddo

end subroutine TGrid_average

subroutine UGrid_average(Ng1,Ng2,Uc,Uf,M,N)

implicit none

real*8, dimension(M,N) :: Uc
real*8, dimension(M*rat**(Ng2-Ng1),N*rat**(Ng2-Ng1)) :: Uf
real*8, dimension(M*rat**(Ng2-Ng1),rat**(Ng2-Ng1-1)) :: UfB
integer :: Ng1,Ng2
integer :: M,N
integer :: i,j,r
integer :: diff

diff=Ng2-Ng1

UfB=Uf(:,N*rat**diff-(rat**(diff-1))+1:N*rat**diff)

j=1

do i=1,M
	r=rat**(diff-1)
	Uc(i,j)=(1.0/((rat**diff)*(1+rat**diff)))*sum(Uf(((i-1)*rat**diff)+1:(i*rat**diff),1:r+1))+&
			(1.0/((rat**diff)*(1+rat**diff)))*sum(UfB(((i-1)*rat**diff)+1:(i*rat**diff),1:r))
enddo

do i=1,M
	do j=2,N
		r=(rat**(diff-1))
		Uc(i,j)=(1.0/((rat**diff)*(1+rat**diff)))*sum(Uf(((i-1)*rat**diff)+1:(i*rat**diff),((j-1)*rat**diff)+1-r:((j-1)*rat**diff)+1+r))
	enddo
enddo


end subroutine UGrid_average


subroutine WGrid_average(Ng1,Ng2,Wc,Wf,M,N)

implicit none

real*8, dimension(M,N) :: Wc
real*8, dimension(M*rat**(Ng2-Ng1),N*rat**(Ng2-Ng1)) :: Wf
integer :: Ng1,Ng2
integer :: M,N
integer :: i,j,r
integer :: diff

diff=Ng2-Ng1

do i=1,M-1
	do j=1,N
		r=(rat**(diff-1))
		Wc(i,j)=(1.0/((rat**diff)*(1+rat**diff)))*sum(Wf((i*rat**diff)+1-r:(i*rat**diff)+1+r,((j-1)*rat**diff)+1:(j*rat**diff)))
	enddo
enddo

i=M

do j=1,N
	Wc(i,j)=(1.0/(rat**diff))*(sum(Wf(M*rat**diff,(j-1)*rat**diff+1:j*rat**diff)))
enddo

end subroutine WGrid_average


subroutine UGrid_cons_average(Ng1,Ng2,Uc,Uf,M,N)

implicit none

real*8, dimension(M,N) :: Uc
real*8, dimension(M*rat**(Ng2-Ng1),N*rat**(Ng2-Ng1)) :: Uf
real*8, dimension(M*rat**(Ng2-Ng1),rat**(Ng2-Ng1-1)) :: UfB
integer :: Ng1,Ng2
integer :: M,N
integer :: i,j,r
integer :: diff

diff=Ng2-Ng1

do i=1,M
	do j=1,N
		r=(rat**(diff-1))
		Uc(i,j)=(1.0/(rat**diff))*sum(Uf(((i-1)*rat**diff)+1:(i*rat**diff),((j-1)*rat**diff)+1))
	enddo
enddo


end subroutine UGrid_cons_average



subroutine WGrid_cons_average(Ng1,Ng2,Wc,Wf,M,N)

implicit none

real*8, dimension(M,N) :: Wc
real*8, dimension(M*rat**(Ng2-Ng1),N*rat**(Ng2-Ng1)) :: Wf
integer :: Ng1,Ng2
integer :: M,N
integer :: i,j,r
integer :: diff

diff=Ng2-Ng1

do i=1,M
	do j=1,N
		Wc(i,j)=(1.0/rat**diff)*sum(Wf(((i-1)*rat**diff)+rat*diff,((j-1)*rat**diff)+1:j*rat**diff))
	enddo
enddo


end subroutine



subroutine TwoGrid_T_Average(Tc,Tf,M1,N1)

implicit none

real*8, dimension(M1,N1) :: Tc
real*8, dimension(2*M1,2*N1) :: Tf
integer :: M1,M2,N1,N2
integer :: i,j

do i=1,M1
	do j=1,N1
		Tc(i,j)=0.25*(sum(Tf((i-1)*2+1:i*2,(j-1)*2+1:j*2)))
	enddo
enddo

end subroutine

subroutine Grid_average(U_T,W_T,Bu_T,Bs_T,P_T,M,N,grids,grid_num)

implicit none
!type(var_pointer), dimension(N_grids) :: var 
!type(var_pointer), dimension(:), allocatable :: var_T
real*8, dimension(M,N) :: U_T,W_T,Bu_T,Bs_T,P_T
real*8, dimension(:,:), allocatable :: Wc,Uc,Buc,Bsc,Pc,Wf,Uf,Buf,Bsf,Pf
integer :: i,j,M,N,grids,grid_num,Mz
integer :: r,k




allocate(Uf(z_dim_vec_ext(N_grids),x_dim_vec(N_grids)),Wf(z_dim_vec_ext(N_grids),x_dim_vec(N_grids)),&
		Buf(z_dim_vec_ext(N_grids),x_dim_vec(N_grids)),Bsf(z_dim_vec_ext(N_grids),x_dim_vec(N_grids)),&
		Pf(z_dim_vec_ext(N_grids),x_dim_vec(N_grids)))

Uf=var(N_grids)%p_var(NGz+1:NGz+z_dim_vec(N_grids),NGx+1:NGx+x_dim_vec(N_grids),1)
Wf=var(N_grids)%p_var(NGz+1:NGz+z_dim_vec(N_grids),NGx+1:NGx+x_dim_vec(N_grids),2)
Buf=var(N_grids)%p_var(NGz+1:NGz+z_dim_vec(N_grids),NGx+1:NGx+x_dim_vec(N_grids),3)
Bsf=var(N_grids)%p_var(NGz+1:NGz+z_dim_vec(N_grids),NGx+1:NGx+x_dim_vec(N_grids),4)
Pf=var(N_grids)%p_var(NGz+1:NGz+z_dim_vec(N_grids),NGx+1:NGx+x_dim_vec(N_grids),5)

do i=N_grids,grid_num+1,-1

	allocate(Uc(z_dim_vec_ext(i-1),x_dim_vec(i-1)),Wc(z_dim_vec_ext(i-1),x_dim_vec(i-1)),Buc(z_dim_vec_ext(i-1),x_dim_vec(i-1)),&
			Bsc(z_dim_vec_ext(i-1),x_dim_vec(i-1)),Pc(z_dim_vec_ext(i-1),x_dim_vec(i-1)))

	Uc(z_dim_vec_ext(i-1)-z_dim_vec(i-1)+1:z_dim_vec_ext(i-1),:)=&
			var(i-1)%p_var(Ngz+1:NGz+z_dim_vec(i-1),NGx+1:NGx+x_dim_vec(i-1),1)
	Wc(z_dim_vec_ext(i-1)-z_dim_vec(i-1)+1:z_dim_vec_ext(i-1),:)=&
			var(i-1)%p_var(Ngz+1:NGz+z_dim_vec(i-1),NGx+1:NGx+x_dim_vec(i-1),2)
	Buc(z_dim_vec_ext(i-1)-z_dim_vec(i-1)+1:z_dim_vec_ext(i-1),:)=&
			var(i-1)%p_var(Ngz+1:NGz+z_dim_vec(i-1),NGx+1:NGx+x_dim_vec(i-1),3)
	Bsc(z_dim_vec_ext(i-1)-z_dim_vec(i-1)+1:z_dim_vec_ext(i-1),:)=&
			var(i-1)%p_var(Ngz+1:NGz+z_dim_vec(i-1),NGx+1:NGx+x_dim_vec(i-1),4)
	Pc(z_dim_vec_ext(i-1)-z_dim_vec(i-1)+1:z_dim_vec_ext(i-1),:)=&
			var(i-1)%p_var(Ngz+1:NGz+z_dim_vec(i-1),NGx+1:NGx+x_dim_vec(i-1),5)
	
	Mz=z_dim_vec_ext(i-1)-z_dim_vec(i-1)

	call TGrid_average(1,2,Buc(1:Mz,:),Buf,Mz,x_dim_vec(i-1))
	call TGrid_average(1,2,Bsc(1:Mz,:),Bsf,Mz,x_dim_vec(i-1))
	call TGrid_average(1,2,Pc(1:Mz,:),Pf,Mz,x_dim_vec(i-1))
	call UGrid_cons_average(1,2,Uc(1:Mz,:),Uf,Mz,x_dim_vec(i-1))
	call WGrid_cons_average(1,2,Wc(1:Mz,:),Wf,Mz,x_dim_vec(i-1))

	deallocate(Uf,Wf,Buf,Bsf,Pf)
	allocate(Uf(z_dim_vec_ext(i-1),x_dim_vec(i-1)),Wf(z_dim_vec_ext(i-1),x_dim_vec(i-1)),Buf(z_dim_vec_ext(i-1),x_dim_vec(i-1)),&
			Bsf(z_dim_vec_ext(i-1),x_dim_vec(i-1)),Pf(z_dim_vec_ext(i-1),x_dim_vec(i-1)))
	Uf=Uc
	Wf=Wc
	Buf=Buc
	Bsf=Bsc
	Pf=Pc
	deallocate(Uc,Wc,Buc,Bsc,Pc)
enddo

U_T=Uf
W_T=Wf
Bu_T=Buf
Bs_T=Bsf
P_T=Pf

deallocate(Uf,Wf,Buf,Bsf,Pf)

end subroutine

end module average


