module interp

use grid

contains

subroutine InterpU_coarse_fine(Uc,Uf,zc,zf,N)

implicit none

real*8, dimension(:,:) :: Uc
real*8, dimension(:) :: Uf
real*8, dimension(:) :: zc
real*8 :: zf
integer :: N
real*8 :: p1,p2,m,y,z
integer :: i

do i=1,N-1
	call LineInterp(Uc(1,i),Uc(2,i),y,zc(1),zc(2),zf)
	Uf(2*i-1)=y
	call LineInterp(0.5*(Uc(1,i)+Uc(1,i+1)),0.5*(Uc(2,i)+Uc(2,i+1)),y,zc(1),zc(2),zf)
	Uf(2*i)=y
enddo

call LineInterp(Uc(1,N),Uc(2,N),y,zc(1),zc(2),zf)
Uf(2*N-1)=y
call LineInterp(0.5*(Uc(1,N)+Uc(1,1)),0.5*(Uc(2,N)+Uc(2,1)),y,zc(1),zc(2),zf)
Uf(2*N)=y

end subroutine


subroutine InterpU_fine_coarse(Uc,Uf,N)

implicit none

real*8, dimension(:,:) :: Uf
real*8, dimension(:) :: Uc
integer :: N
integer :: i

do i=1,N
	Uc(i)=0.5*(Uf(1,2*i-1)+Uf(2,2*i-1))
enddo

end subroutine

subroutine InterpT_fine_coarse(Tc,Tf,N)

implicit none

real*8, dimension(2,rat*N) :: Tf
real*8, dimension(N) :: Tc
integer :: N,i

do i=1,N
	Tc(i)=0.25*(Tf(1,2*i-1)+Tf(1,2*i-1)+Tf(2,2*i)+Tf(2,2*i))
enddo

end subroutine


subroutine InterpT_coarse_fine(Uc,Uf,xc,xf,z1,z2,zf,N)

implicit none

real*8, dimension(2,N) :: Uc
real*8, dimension(rat*N) :: Uf,xf
real*8, dimension(N) :: xc
real*8 :: z1,z2,zf
integer :: N,i
real*8 :: y1,hx

hx=xc(2)-xc(1)

call BilinearInterp_coarse_fine(Uc(1,N),Uc(1,1),Uc(2,N),Uc(2,1),y1,xc(1)-hx,xc(1),z1,z2,xf(1),zf) 
Uf(1)=y1

do i=1,N-1
	call BilinearInterp_coarse_fine(Uc(1,i),Uc(1,i+1),Uc(2,i),Uc(2,i+1),y1,xc(i),xc(i+1),z1,z2,xf(2*i),zf)
	Uf(2*i)=y1
	call BilinearInterp_coarse_fine(Uc(1,i),Uc(1,i+1),Uc(2,i),Uc(2,i+1),y1,xc(i),xc(i+1),z1,z2,xf(2*i+1),zf)
	Uf(2*i+1)=y1
enddo

call BilinearInterp_coarse_fine(Uc(1,N),Uc(1,1),Uc(2,N),Uc(2,1),y1,xc(N),xc(N)+hx,z1,z2,xf(2*N),zf)
Uf(2*N)=y1

!Uf(1)=Uf(2)
!Uf(2*N)=Uf(2*N-1)

end subroutine

subroutine InterpW_coarse_fine(Wc,Wf,zc,zf,xc,xf,N)

implicit none

real*8, dimension(:,:) :: Wc
real*8, dimension(:) :: Wf,zc,xc,xf
real*8 :: zf
real*8 :: y1,hx
integer :: N,i

hx=xc(2)-xc(1)

!call BilinearInterp_coarse_fine(Wc(1,N),Wc(1,1),Wc(2,N),Wc(2,1),y1,xc(1)-hx,xc(1),zc(1),zc(2),xf(1),zf)
!Wf(1)=y1

do i=1,N-1
	call BilinearInterp_coarse_fine(Wc(1,i),Wc(1,i+1),Wc(2,i),Wc(2,i+1),y1,xc(i),xc(i+1),zc(1),zc(2),xf(2*i),zf)
	Wf(2*i)=y1
	call BilinearInterp_coarse_fine(Wc(1,i),Wc(1,i+1),Wc(2,i),Wc(2,i+1),y1,xc(i),xc(i+1),zc(1),zc(2),xf(2*i+1),zf)
	Wf(2*i+1)=y1
enddo

!call BilinearInterp_coarse_fine(Wc(1,N),Wc(1,1),Wc(2,N),Wc(2,1),y1,xc(N),xc(N)+hx,zc(1),zc(2),xf(2*N),zf)
!Wf(2*N)=y1
wf(1)=wf(2)
wf(2*N)=wf(2*N-1)

end subroutine


subroutine InterpW_fine_coarse(Uc,Uf,N)

implicit none

real*8, dimension(:) :: Uc,Uf
integer :: N,i

do i=1,N
	Uc(i)=0.5*(Uf(2*i-1)+Uf(2*i))
enddo


end subroutine

subroutine LineInterp(p1,p2,y1,x1,x2,x3)

implicit none

real*8 :: p1,p2
real*8 :: y1
real*8 :: x1,x2,h,x3 !x1/x2 are the coords of p1/p2 and x3 is the coord of y1
real*8 :: m

m=(p2-p1)/(x2-x1)

y1=p1+m*(x3-x1)

end subroutine

subroutine BilinearInterp_coarse_fine(pLL,pLR,pTl,pTR,y1,x1,x2,z1,z2,x3,z3)

implicit none

real*8 :: pLL,pLR,pTL,pTR
real*8 :: x1,x2,z1,z2,x3,z3
real*8 :: y1,y2
real*8 ::p1,p2,m1,m2,m3
integer :: i

m1=(pLR-pLL)/(x2-x1)
m2=(pTR-PTL)/(x2-x1)

p1=pLL+m1*(x3-x1)
p2=pTL+m2*(x3-x1)

m3=(p2-p1)/(z2-z1)

y1=p1+m3*(z3-z1)

end subroutine


subroutine BilinearInterp_fine_coarse(x1,x2,x3,x4,y1)

implicit none

real*8 :: x1,x2,x3,x4
real*8 :: y1
integer :: i

y1=(0.25)*(x1+x2+x3+x4)


end subroutine

subroutine CubicInterpT_coarse_fine(Tc,Tf_interp,Tf,xc,xf,zc1,zc2,zf1,zf2,x0,z0,M,N,hx,hz)

implicit none

real*8, dimension(2,N) :: Tc
real*8, dimension(2,rat*N) :: Tf_interp
real*8, dimension(rat*N) :: Tf
real*8, dimension(N) :: xc
real*8, dimension(rat*N) :: xf,x0
real*8 :: zc1,zc2,zf1,zf2,z0,hx,hz
integer :: M,N
integer :: i
			

call CubicInterp_coarse_fine((/Tc(1,1),Tc(1,2),Tc(1,3),Tc(1,4)/),(/Tc(2,1),Tc(2,2),Tc(2,3),Tc(2,4)/),&
											(/Tf_interp(1,1),Tf_interp(1,2),Tf_interp(1,3),Tf_interp(1,4)/),&
											(/Tf_interp(2,1),Tf_interp(2,2),Tf_interp(2,3),Tf_interp(2,4)/),&
											(/xc(1),xc(2),xc(3),xc(4)/),(/xc(1),xc(2),xc(3),xc(4)/),&
											(/xf(1),xf(2),xf(3),xf(4)/),(/xf(1),xf(2),xf(3),xf(4)/),&
											zc1,zc2,zf1,zf2,x0(1),z0,Tf(1))

call CubicInterp_coarse_fine((/Tc(1,1),Tc(1,2),Tc(1,3),Tc(1,4)/),(/Tc(2,1),Tc(2,2),Tc(2,3),Tc(2,4)/),&
											(/Tf_interp(1,1),Tf_interp(1,2),Tf_interp(1,3),Tf_interp(1,4)/),&
											(/Tf_interp(2,1),Tf_interp(2,2),Tf_interp(2,3),Tf_interp(2,4)/),&
											(/xc(1),xc(2),xc(3),xc(4)/),(/xc(1),xc(2),xc(3),xc(4)/),&
											(/xf(1),xf(2),xf(3),xf(4)/),(/xf(1),xf(2),xf(3),xf(4)/),&
											zc1,zc2,zf1,zf2,x0(2),z0,Tf(2))				

do i=2,N-2
	call CubicInterp_coarse_fine(Tc(1,i-1:i+2),Tc(2,i-1:i+2),Tf_interp(1,2*i-2:2*i+1),Tf_interp(2,2*i-2:2*i+1),&
												xc(i-1:i+2),xc(i-1:i+2),xf(2*i-2:2*i+1),xf(2*i-2:2*i+1),zc1,zc2,zf1,zf2,x0(2*i-1),z0,Tf(2*i-1))
	call CubicInterp_coarse_fine(Tc(1,i-1:i+2),Tc(2,i-1:i+2),Tf_interp(1,2*i-2:2*i+1),Tf_interp(2,2*i-2:2*i+1),&
												xc(i-1:i+2),xc(i-1:i+2),xf(2*i-2:2*i+1),xf(2*i-2:2*i+1),zc1,zc2,zf1,zf2,x0(2*i),z0,Tf(2*i))
enddo

call CubicInterp_coarse_fine((/Tc(1,N-3),Tc(1,N-2),Tc(1,N-1),Tc(1,N)/),(/Tc(2,N-3),Tc(2,N-2),Tc(2,N-1),Tc(2,N)/),&
											(/Tf_interp(1,rat*N-3),Tf_interp(1,rat*N-2),Tf_interp(1,rat*N-1),Tf_interp(1,rat*N)/),&
											(/Tf_interp(2,rat*N-3),Tf_interp(2,rat*N-2),Tf_interp(2,rat*N-1),Tf_interp(2,rat*N)/),&
											(/xc(N-3),xc(N-2),xc(N-1),xc(N)/),(/xc(N-3),xc(N-2),xc(N-1),xc(N)/),&
											(/xf(rat*N-3),xf(rat*N-2),xf(rat*N-1),xf(rat*N)/),(/xf(rat*N-3),xf(rat*N-2),xf(rat*N-1),xf(rat*N)/),&
											zc1,zc2,zf1,zf2,x0(rat*N-3),z0,Tf(rat*N-3))

call CubicInterp_coarse_fine((/Tc(1,N-3),Tc(1,N-2),Tc(1,N-1),Tc(1,N)/),(/Tc(2,N-3),Tc(2,N-2),Tc(2,N-1),Tc(2,N)/),&
											(/Tf_interp(1,rat*N-3),Tf_interp(1,rat*N-2),Tf_interp(1,rat*N-1),Tf_interp(1,rat*N)/),&
											(/Tf_interp(2,rat*N-3),Tf_interp(2,rat*N-2),Tf_interp(2,rat*N-1),Tf_interp(2,rat*N)/),&
											(/xc(N-3),xc(N-2),xc(N-1),xc(N)/),(/xc(N-3),xc(N-2),xc(N-1),xc(N)/),&
											(/xf(rat*N-3),xf(rat*N-2),xf(rat*N-1),xf(rat*N)/),(/xf(rat*N-3),xf(rat*N-2),xf(rat*N-1),xf(rat*N)/),&
											zc1,zc2,zf1,zf2,x0(rat*N-2),z0,Tf(rat*N-2))

call CubicInterp_coarse_fine((/Tc(1,N-3),Tc(1,N-2),Tc(1,N-1),Tc(1,N)/),(/Tc(2,N-3),Tc(2,N-2),Tc(2,N-1),Tc(2,N)/),&
											(/Tf_interp(1,rat*N-3),Tf_interp(1,rat*N-2),Tf_interp(1,rat*N-1),Tf_interp(1,rat*N)/),&
											(/Tf_interp(2,rat*N-3),Tf_interp(2,rat*N-2),Tf_interp(2,rat*N-1),Tf_interp(2,rat*N)/),&
											(/xc(N-3),xc(N-2),xc(N-1),xc(N)/),(/xc(N-3),xc(N-2),xc(N-1),xc(N)/),&
											(/xf(rat*N-3),xf(rat*N-2),xf(rat*N-1),xf(rat*N)/),(/xf(rat*N-3),xf(rat*N-2),xf(rat*N-1),xf(rat*N)/),&
											zc1,zc2,zf1,zf2,x0(rat*N-1),z0,Tf(rat*N-1))


call CubicInterp_coarse_fine((/Tc(1,N-3),Tc(1,N-2),Tc(1,N-1),Tc(1,N)/),(/Tc(2,N-3),Tc(2,N-2),Tc(2,N-1),Tc(2,N)/),&
											(/Tf_interp(1,rat*N-3),Tf_interp(1,rat*N-2),Tf_interp(1,rat*N-1),Tf_interp(1,rat*N)/),&
											(/Tf_interp(2,rat*N-3),Tf_interp(2,rat*N-2),Tf_interp(2,rat*N-1),Tf_interp(2,rat*N)/),&
											(/xc(N-3),xc(N-2),xc(N-1),xc(N)/),(/xc(N-3),xc(N-2),xc(N-1),xc(N)/),&
											(/xf(rat*N-3),xf(rat*N-2),xf(rat*N-1),xf(rat*N)/),(/xf(rat*N-3),xf(rat*N-2),xf(rat*N-1),xf(rat*N)/),&
											zc1,zc2,zf1,zf2,x0(rat*N),z0,Tf(rat*N))


end subroutine CubicInterpT_coarse_fine


subroutine CubicInterpT_fine_coarse(Tc_interp,Tf,Tc,xc,xf,zc1,zc2,zf1,zf2,x0,z0,M,N,hx,hz)

implicit none

real*8, dimension(2,N) :: Tc_interp
real*8, dimension(2,rat*N) :: Tf
real*8, dimension(N) :: Tc
real*8, dimension(rat*N) :: xf
real*8, dimension(N) :: xc,x0
real*8 :: zc1,zc2,zf1,zf2,z0,hx,hz
integer :: M,N
integer :: i


call CubicInterp_coarse_fine((/Tc_interp(1,1),Tc_interp(1,2),Tc_interp(1,3),Tc_interp(1,4)/),&
											(/Tc_interp(2,1),Tc_interp(2,2),Tc_interp(2,3),Tc_interp(2,4)/),&
											(/Tf(1,1),Tf(1,2),Tf(1,3),Tf(1,4)/),&
											(/Tf(2,1),Tf(2,2),Tf(2,3),Tf(2,4)/),&
											(/xc(1),xc(2),xc(3),xc(4)/),(/xc(1),xc(2),xc(3),xc(4)/),&
											(/xf(1),xf(2),xf(3),xf(4)/),(/xf(1),xf(2),xf(3),xf(4)/),&
											zc1,zc2,zf1,zf2,x0(1),z0,Tc(1))

			

do i=2,N-1
	call CubicInterp_coarse_fine(Tc_interp(1,i-1:i+2),Tc_interp(2,i-1:i+2),Tf(1,2*i-2:2*i+1),Tf(2,2*i-2:2*i+1),&
												xc(i-1:i+2),xc(i-1:i+2),xf(2*i-2:2*i+1),xf(2*i-2:2*i+1),zc1,zc2,zf1,zf2,x0(i),z0,Tc(i))
enddo

call CubicInterp_coarse_fine((/Tc_interp(1,N-3),Tc_interp(1,N-2),Tc_interp(1,N-1),Tc_interp(1,N)/),&
											(/Tc_interp(2,N-3),Tc_interp(2,N-2),Tc_interp(2,N-1),Tc_interp(2,N)/),&
											(/Tf(1,rat*N-3),Tf(1,rat*N-2),Tf(1,rat*N-1),Tf(1,rat*N)/),&
											(/Tf(2,rat*N-3),Tf(2,rat*N-2),Tf(2,rat*N-1),Tf(2,rat*N)/),&
											(/xc(N-3),xc(N-2),xc(N-1),xc(N)/),(/xc(N-3),xc(N-2),xc(N-1),xc(N)/),&
											(/xf(rat*N-3),xf(rat*N-2),xf(rat*N-1),xf(rat*N)/),(/xf(rat*N-3),xf(rat*N-2),xf(rat*N-1),xf(rat*N)/),&
											zc1,zc2,zf1,zf2,x0(N),z0,Tc(N))




end subroutine CubicInterpT_fine_coarse

subroutine CubicInterp_fine_coarse_avg(Tf,Tc,xc,xf,zf1,zf2,zf3,zf4,z0,M,N)

implicit none

real*8, dimension(4,rat*N+2) :: Tf !extra layer of ghost cells on the lateral sides
real*8, dimension(N) :: Tc
real*8, dimension(rat*N) :: xf
real*8, dimension(N) :: xc
real*8 :: zf1,zf2,zf3,zf4,z0
integer :: M,N
integer :: i

do i=1,N
	call CubicInterp_coarse_fine(Tf(2*i-2,2*i-1:2*i+2),Tf(2*i-1,2*i-1:2*i+2),Tf(2*i,2*i-1:2*i+2),Tf(2*i+1,2*i-1:2*i+2),&
								xf(2*i-1),xf(2*i),xf(2*i+1),xf(2*i+2),zf1,zf2,zf3,zf4,xc(i),z0,Tc(i))
									
enddo

end subroutine CubicInterp_fine_coarse_avg

subroutine CubicInterp_coarse_fine(p1,p2,p3,p4,x1,x2,x3,x4,z1,z2,z3,z4,x0,z0,p0)

implicit none

real*8, dimension(4) :: p1,p2,p3,p4,x1,x2,x3,x4
real*8 :: z1,z2,z3,z4,x0,z0,p0
real*8 :: y1,y2,y3,y4
real*8, dimension(4) :: y_vec,z_vec

y1=Cubic_Lagrange(p1,x1,x0)
y2=Cubic_Lagrange(p2,x2,x0)
y3=Cubic_Lagrange(p3,x3,x0)
y4=Cubic_Lagrange(p4,x4,x0)

y_vec(1)=y1
y_vec(2)=y2
y_vec(3)=y3
y_vec(4)=y4

z_vec(1)=z1
z_vec(2)=z2
z_vec(3)=z3
z_vec(4)=z4

p0=Cubic_Lagrange(y_vec,z_vec,z0)

end subroutine CubicInterp_coarse_fine


function Cubic_Lagrange(p,x,x0)

!Computes cubic lagrange polynomial and evaluates at a specified point

implicit none
real*8, dimension(4) :: p,x
real*8 :: x0,Cubic_Lagrange

	Cubic_Lagrange=p(1)*(x0-x(2))*(x0-x(3))*(x0-x(4))/((x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))) &
		+p(2)*(x0-x(1))*(x0-x(3))*(x0-x(4))/((x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))) &
		+p(3)*(x0-x(1))*(x0-x(2))*(x0-x(4))/((x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4))) &
		+p(4)*(x0-x(1))*(x0-x(2))*(x0-x(3))/((x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3)))	


end function Cubic_Lagrange

function quadratic_interp(p1,p2,p3,x1,x2,x3,x0)

implicit none

real*8, intent(in) :: p1,p2,p3,x1,x2,x3,x0
real*8 :: quadratic_interp

quadratic_interp=quadratic_lagrange(p1,p2,p3,x1,x2,x3,x0)

end function quadratic_interp


function Quadratic_Lagrange(p1,p2,p3,x1,x2,x3,x0)

implicit none
real*8, intent(in) :: p1,p2,p3,x1,x2,x3,x0
real*8 :: Quadratic_Lagrange

Quadratic_lagrange=p1*(x0-x2)*(x0-x3)/((x1-x2)*(x1-x3)) + p2*(x0-x1)*(x0-x3)/((x2-x1)*(x2-x3)) &
					+ p3*(x0-x1)*(x0-x2)/((x3-x1)*(x3-x2))

end function Quadratic_Lagrange

end module interp

