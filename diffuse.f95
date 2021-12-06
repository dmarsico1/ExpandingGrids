module diffuse

use grid

contains

subroutine diffVel(Ud,U,Wd,W,M,N,hx,hz)

!This subroutine assumes that U,W already have the appropriate layers of ghost cells

implicit none

real*8, intent(out) :: Ud(M,N),Wd(M-1,N)
real*8, intent(inout) :: W(:,:),U(:,:)
real*8, intent(in) :: hx,hz
integer, intent(in) :: M,N
integer :: i,j

do i=1,M
	do j=1,N
		Ud(i,j)=ev*(1.0/hx**2)*(U(i+1,j+2)-2.0*U(i+1,j+1)+U(i+1,j)) &
 				+ev*(1.0/hz**2)*(U(i+2,j+1)-2.0*U(i+1,j+1)+U(i,j+1))
	enddo
enddo

do i=1,M-1
	do j=1,N
		Wd(i,j)=ev*(1.0/hx**2)*(W(i+1,j+2)-2.0*W(i+1,j+1)+W(i+1,j)) &
 				+ev*(1.0/hz**2)*(W(i+2,j+1)-2.0*W(i+1,j+1)+W(i,j+1))
	enddo
enddo

end subroutine diffVel


subroutine diffB(Bud,Bu,Bsd,Bs,M,N,hx,hz)

!This subroutine already assumes that Bu,Bs already have the approprate layer of ghost cells

implicit none

real*8, intent(out) :: Bud(M,N),Bsd(M,N)
real*8, intent(inout) :: Bu(M+2,N+2),Bs(M+2,N+2)
real*8 :: hx,hz
integer, intent(in) :: M,N
integer :: i,j

do i=1,M
	do j=1,N
		Bud(i,j)=eb*(1.0/hx**2)*(Bu(i+1,j+2)-2.0*Bu(i+1,j+1)+Bu(i+1,j)) &
 				+eb*(1.0/hz**2)*(Bu(i+2,j+1)-2.0*Bu(i+1,j+1)+Bu(i,j+1))
		Bsd(i,j)=eb*(1.0/hx**2)*(Bs(i+1,j+2)-2.0*Bs(i+1,j+1)+Bs(i+1,j)) &
 				+eb*(1.0/hz**2)*(Bs(i+2,j+1)-2.0*Bs(i+1,j+1)+Bs(i,j+1))
	enddo
enddo

end subroutine diffB

subroutine Horiz_HyperDiff(U,W,Uhd,Whd,M,N,hx)

implicit none

real*8, dimension(M,N) :: Uhd
real*8, dimension(M,N+4) :: U
real*8, dimension(M-1,N+4) :: W
real*8, dimension(M-1,N) :: Whd
integer, intent(in) :: M,N
real*8 :: hx
integer :: i,j

do i=NGz+1,NGz+M
	do j=NGx+1,NGx+N
		Uhd(i-NGx,j-NGx)=-mu1*(1.0/hx**4)*(U(i,j-2)-4.0*U(i,j-1)+6.0*U(i,j)-4.0*U(i,j+1)+U(i,j+2))
	enddo
enddo

do i=NGz+1,NGz+M-1
	do j=NGx+1,NGx+N
		Whd(i-NGx,j-NGx)=-mu1*(1.0/hx**4)*(W(i,j-2)-4.0*W(i,j-1)+6.0*W(i,j)-4.0*W(i,j+1)+W(i,j+2))
	enddo
enddo

end subroutine Horiz_HyperDiff


subroutine Horiz_HyperDiffB(Bu,Bs,Buhd,Bshd,M,N,hx)

real*8, dimension(M,N) :: Buhd,Bshd
real*8, dimension(M,N+4) :: Bu,Bs
integer, intent(in) :: M,N
real*8 :: hx
integer :: i,j


do i=NGz+1,NGz+M
	do j=NGx+1,NGx+N
		Buhd(i-NGx,j-NGx)=-mu1*(1.0/hx**4)*(Bu(i,j-2)-4.0*Bu(i,j-1)+6.0*Bu(i,j)-4.0*Bu(i,j+1)+Bu(i,j+2))
		Bshd(i-NGx,j-NGx)=-mu1*(1.0/hx**4)*(Bs(i,j-2)-4.0*Bs(i,j-1)+6.0*Bs(i,j)-4.0*Bs(i,j+1)+Bs(i,j+2))
	enddo
enddo


end subroutine Horiz_HyperDiffB


subroutine Normal_viscVel(Ud,U,Wd,W,M,N,hz)

!This subroutine assumes that U,W already have the appropriate layers of ghost cells

implicit none

real*8, intent(out) :: Ud(M,N),Wd(M-1,N)
real*8, intent(inout) :: W(M+1,N),U(M+2,N)
real*8, intent(in) :: hz
integer, intent(in) :: M,N
integer :: i,j

do i=1,M
	do j=1,N
		Ud(i,j)=mu2*(1.0/hz**2)*(U(i+2,j+1)-2.0*U(i+1,j+1)+U(i,j+1))
	enddo
enddo

do i=1,M-1
	do j=1,N
		Wd(i,j)=mu2*(1.0/hz**2)*(W(i+2,j+1)-2.0*W(i+1,j+1)+W(i,j+1))
	enddo
enddo


end subroutine Normal_viscVel

subroutine Normal_viscB(Bud,Bu,Bsd,Bs,M,N,hz)

!This subroutine already assumes that Bu,Bs already have the approprate layer of ghost cells

implicit none

real*8, intent(out) :: Bud(M,N),Bsd(M,N)
real*8, intent(inout) :: Bu(M+2,N+2),Bs(M+2,N+2)
real*8 :: hx,hz
integer, intent(in) :: M,N
integer :: i,j

do i=1,M
	do j=1,N
		Bud(i,j)=mu2*(1.0/hz**2)*(Bu(i+2,j+1)-2.0*Bu(i+1,j+1)+Bu(i,j+1))
		Bsd(i,j)=mu2*(1.0/hz**2)*(Bs(i+2,j+1)-2.0*Bs(i+1,j+1)+Bs(i,j+1))
	enddo
enddo

end subroutine Normal_viscB

end module diffuse
