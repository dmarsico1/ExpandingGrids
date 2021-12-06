module matrix

use grid
use Boundary

contains

function SumMat(M)

implicit none

real*8 :: M(:,:)
integer :: numrows,numcols,i,j
real*8 :: sum
real*8 :: SumMat
	
numrows=size(M,1)
numcols=size(M,2)
sum=0

do i=1,numrows
	do j=1,numcols
		sum=sum+M(i,j)
	enddo
enddo

SumMat=sum

end function SumMat

subroutine MatMult(p,q,M,N,hx,hz)

implicit none
real*8, dimension(M+2,N+2), intent(in) :: q !this is what gets INPUT.  Layer of ghost cells around q.
real*8, dimension(M,N), intent(out) :: p !this is what gets OUTPUT.  No ghost cells around p.
real*8, intent(in) :: hx,hz
integer, intent(in) :: M,N
integer :: i,j

do i=2,M+1
	do j=2,N+1
		p(i-1,j-1)=(1.0/hx**2)*(q(i,j+1)-2.0*q(i,j)+q(i,j-1))+(1.0/hz**2)*(q(i+1,j)-2.0*q(i,j)+q(i-1,j))
	enddo
enddo


end subroutine MatMult

subroutine ConjGrad(x,f,BcT,BcB,tol,N,M,hx,hz,dt,grid_num,fu,fl)

implicit none
real*8, dimension(M+2,N+2), intent(out) :: x !This is what gets output
real*8, dimension(M,N), intent(in) :: f !This is the right hand side
real*8, intent(in) :: tol
integer :: i,j,k
real*8, dimension(M,N) :: r,q,Ax
real*8, intent(in), dimension(N) :: fu,fl
real*8, dimension(N) :: fu_temp
real*8, dimension(M+2,N+2) :: p
real*8 :: rold,rnew,alpha
integer :: N,M,BcT,BcB,grid_num
real*8 :: hx,hz,dt


r=f

x=0

if(grid_num==1)then
	fu_temp=0
endif

r(M,:)=r(M,:)-fu
r(1,:)=r(1,:)+fl

r=r-(1.0/(M*N))*summat(r)


if (SumMat(r*r)<10.0**(-10))then
	do i=1,M
		do j=1,N
			x(i,j)=0.0
		enddo
	enddo
else

p(2:M+1,2:N+1)=r

rold=SumMat(r*r)

do k=1,N**2
	call HomogeneousBoundary(p,BcT,BcB,M,N)
	call MatMult(q,p,M,N,hx,hz)
	alpha=rold/(SumMat(p(2:M+1,2:N+1)*q))
	x(2:M+1,2:N+1) = x(2:M+1,2:N+1) + alpha*p(2:M+1,2:N+1)
	r = r - alpha*q
	rnew = SumMat(r*r)
	if (sqrt(rnew) < tol) then
		exit
	end if
	p(2:M+1,2:N+1) = r + (rnew/rold)*p(2:M+1,2:N+1)
	rold=rnew
enddo

endif

	
end subroutine ConjGrad

end module matrix
