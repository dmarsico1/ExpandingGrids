module advect

use grid

contains

subroutine advectB(Bua,Bsa,Bu,Bs,U,W,M,N,hx,hz)

implicit none

real*8, dimension(:,:) :: Bu,Bs,U
real*8, dimension(:,:) :: W
real*8, intent(out) :: Bua(M,N),Bsa(M,N)
real*8 , dimension(M,N) :: U_Bux,W_Buz,U_Bsx,W_Bsz
integer, intent(in) :: M,N
real*8, intent(in) :: hx,hz
real*8 :: Uavg,Wavg
integer :: i,j

do i=3,M+2
	do j=3,N+2
		Uavg=0.5*(U(i,j+1)+U(i,j))
		U_Bux(i-2,j-2)=Uavg*(Bu(i,j+1)-Bu(i,j-1))/(2.0*hx)
		U_Bux(i-2,j-2)=Uavg*(1.0/(12.0*hx))*(-Bu(i,j+2)+8.0*(Bu(i,j+1)-Bu(i,j-1))+Bu(i,j-2)) &
					+(1.0/(12.0*hx))*abs(Uavg)*(Bu(i,j+2)-4.0*Bu(i,j+1)+6.0*Bu(i,j)-4.0*Bu(i,j-1)+Bu(i,j-2))

	enddo
enddo

do i=3,M+2
	do j=3,N+2
		Wavg=0.5*(W(i-1,j)+W(i,j))
		W_Buz(i-2,j-2)=Wavg*(Bu(i+1,j)-Bu(i-1,j))/(2.0*hz)
		W_Buz(i-2,j-2)=Wavg*(1.0/(12.0*hz))*(-Bu(i+2,j)+8.0*(Bu(i+1,j)-Bu(i-1,j))+Bu(i-2,j)) &
					+(1.0/(12.0*hz))*abs(Wavg)*(Bu(i+2,j)-4.0*Bu(i+1,j)+6.0*Bu(i,j)-4.0*Bu(i-1,j)+Bu(i-2,j))

	enddo
enddo

do i=3,M+2
	do j=3,N+2
		Uavg=0.5*(U(i,j+1)+U(i,j))
		U_Bsx(i-2,j-2)=Uavg*(Bs(i,j+1)-Bs(i,j-1))/(2.0*hx)
		U_Bsx(i-2,j-2)=Uavg*(1.0/(12.0*hx))*(-Bs(i,j+2)+8.0*(Bs(i,j+1)-Bs(i,j-1))+Bs(i,j-2)) &
					+(1.0/(12.0*hx))*abs(Uavg)*(Bs(i,j+2)-4.0*Bs(i,j+1)+6.0*Bs(i,j)-4.0*Bs(i,j-1)+Bs(i,j-2))

	enddo
enddo

do i=3,M+2
	do j=3,N+2
		Wavg=0.5*(W(i-1,j)+W(i,j))
		W_Bsz(i-2,j-2)=Wavg*(Bs(i+1,j)-Bs(i-1,j))/(2.0*hz)
		W_Bsz(i-2,j-2)=Wavg*(1.0/(12.0*hz))*(-Bs(i+2,j)+8.0*(Bs(i+1,j)-Bs(i-1,j))+Bs(i-2,j)) &
					+(1.0/(12.0*hz))*abs(Wavg)*(Bs(i+2,j)-4.0*Bs(i+1,j)+6.0*Bs(i,j)-4.0*Bs(i-1,j)+Bs(i-2,j))

	enddo
enddo

Bua=U_Bux+W_Buz
Bsa=U_Bsx+W_Bsz

end subroutine advectB

subroutine advectVel(Ua,Wa,U,W,M,N,hx,hz)

!U,W have the appropriate layers of ghost cells surruonding them

implicit none

real*8, dimension(M+4,N+4) :: U
real*8, dimension(M+3,N+4) :: W
real*8, intent(out) :: Ua(M,N),Wa(M-1,N)
real*8 :: U_Ux(M,N),W_Uz(M,N),U_Wx(M-1,N),W_Wz(M-1,N)
real*8 :: hx,hz
integer :: M,N
integer :: i,j
real*8 :: Wbar,Ubar


do i=3,M+2
	do j=3,N+2
		U_Ux(i-2,j-2)=U(i,j)*(U(i,j+1)-U(i,j-1))/(2.0*hx)
		U_Ux(i-2,j-2)=U(i,j)*(1.0/(12.0*hx))*(-U(i,j+2)+8.0*(U(i,j+1)-U(i,j-1))+U(i,j-2)) &
					+(1.0/(12.0*hx))*abs(U(i,j))*(U(i,j+2)-4.0*U(i,j+1)+6.0*U(i,j)-4.0*U(i,j-1)+U(i,j-2))

	enddo
enddo

do i=3,M+2
	do j=3,N+2
		Wbar=(1.0/4.0)*(W(i-1,j)+W(i-1,j-1)+W(i,j)+W(i,j-1))
		W_Uz(i-2,j-2)=Wbar*(U(i+1,j)-U(i-1,j))/(2.0*hz)
		W_Uz(i-2,j-2)=Wbar*(1.0/(12.0*hz))*(-U(i+2,j)+8.0*(U(i+1,j)-U(i-1,j))+U(i-2,j)) &
					+(1.0/(12.0*hz))*abs(Wbar)*(U(i+2,j)-4.0*U(i+1,j)+6.0*U(i,j)-4.0*U(i-1,j)+U(i-2,j))

	enddo
enddo

do i=3,M+1
	do j=3,N+2
		Ubar=(1.0/4.0)*(U(i,j)+U(i,j+1)+U(i+1,j)+U(i+1,j+1))
		U_Wx(i-2,j-2)=Ubar*(W(i,j+1)-W(i,j-1))/(2.0*hx)
		U_Wx(i-2,j-2)=Ubar*(1.0/(12.0*hx))*(-W(i,j+2)+8.0*(W(i,j+1)-W(i,j-1))+W(i,j-2)) &
					+(1.0/(12.0*hx))*abs(Ubar)*(W(i,j+2)-4.0*W(i,j+1)+6.0*W(i,j)-4.0*W(i,j-1)+W(i,j-2))
	enddo
enddo

do i=3,M+1
	do j=3,N+2
		W_Wz(i-2,j-2)=0.5*(W(i,j+1)+W(i,j-1))*(W(i+1,j)-W(i-1,j))/(2.0*hz)
		W_Wz(i-2,j-2)=W(i,j)*(1.0/(12.0*hz))*(-W(i+2,j)+8.0*(W(i+1,j)-W(i-1,j))+W(i-2,j)) &
					+(1.0/(12.0*hz))*abs(W(i,j))*(W(i+2,j)-4.0*W(i+1,j)+6.0*W(i,j)-4.0*W(i-1,j)+W(i-2,j))
	enddo
enddo


Ua=U_Ux+W_Uz
Wa=U_Wx+W_Wz

end subroutine advectVel


end module advect
