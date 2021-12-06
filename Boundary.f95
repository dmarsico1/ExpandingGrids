module Boundary

use grid
use interp

contains

subroutine Fine_Coarse_VertBoundary(Uc,Wc,Buc,Bsc,Uf,Wf,Buf,Bsf,M,N,bound)

!The coarse arguments have the appropriate layers of ghost cells surrounding them

implicit none
real*8, dimension(:,:) :: Uc,Wc,Buc,Bsc,Uf,Wf,Buf,Bsf
real*8 :: hx,hz
integer :: M,N
integer :: bound
integer :: i

do i=1,NGz
	call InterpU_fine_coarse(Uc(i,NGx+1:N+NGx),Uf(2*i-1:2*i,:),N) !Need NGz to be even here?
	call InterpW_fine_coarse(Wc(i,NGx+1:N+NGx),Wf(2*i-1,:),N)
	call InterpT_fine_coarse(Buc(i,NGx+1:N+NGx),Buf(2*i-1:2*i,:),N)
	call InterpT_fine_coarse(Bsc(i,NGx+1:N+NGx),Bsf(2*i-1:2*i,:),N)
enddo

end subroutine Fine_Coarse_VertBoundary

subroutine Coarse_Fine_VertBoundary(Uc,Wc,Buc,Bsc,Uf,Wf,Buf,Bsf,xtc,xtf,xmc,xmf,ztc,ztf,zmc,zmf,Mf,N,hx,hz)

!The fine arguments have the appropriate layer of ghost cells surrounding them

implicit none
real*8, dimension(2,N) :: Uc,Wc,Buc,Bsc
real*8, dimension(:,:) :: Uf,Wf,Buf,Bsf
real*8, dimension(N) :: xtc,xmc
real*8, dimension(rat*N) :: xtf,xmf
real*8, dimension(:) :: ztc,zmc
real*8, dimension(:) :: ztf,zmf
integer :: Mf,N
real*8 :: hx,hz
integer :: i

do i=1,NGz
	call InterpU_coarse_fine(Uc,Uf(Mf+NGz+i,NGx+1:NGx+rat*N),ztc(1:2),ztf(Mf)+0.5*real(i)*hz,N)
	call InterpT_coarse_fine(Buc,Buf(Mf+NGz+i,NGx+1:NGx+rat*N),xtc,xtf,ztc(1),ztc(2),ztf(Mf)+0.5*real(i)*hz,N)
	call InterpT_coarse_fine(Bsc,Bsf(Mf+NGz+i,NGx+1:NGx+rat*N),xtc,xtf,ztc(1),ztc(2),ztf(Mf)+0.5*real(i)*hz,N)
enddo

call InterpW_coarse_fine(Wc,Wf(Mf+NGz+1,NGx+1:NGx+rat*N),zmc(1:2),zmf(Mf)+0.5*hz,xtc,xtf,N)

end subroutine Coarse_Fine_VertBoundary


subroutine CubicVertBoundary_coarse_fine

implicit none



end subroutine CubicVertBoundary_coarse_fine

subroutine UpperBoundary(U,W,Bu,Bs,UT,BuT,BsT,hx,hz,M,N)

implicit none

real*8, dimension(M+4,N+4) :: U,Bu,Bs
real*8, dimension(M+4,N+4) :: W
real*8, dimension(N) :: UT,BuT,BsT
real*8 :: hx,hz
integer :: M,N
integer :: i

do i=1,N
	U(M+2*NGz-1,NGx+1:N+NGx)=UT
	W(M+2*NGz-2,NGx+1:NGx+N)=0
	Bu(M+2*NGz-1,NGx+1:N+NGx)=BuT
	Bs(M+2*Ngz-1,NGx+1:N+NGx)=BsT
enddo

do i=1,N
	U(M+2*NGz,NGx+1:N+NGx)=UT
	W(M+2*NGz-1,NGx+1:NGx+N)=0
	Bu(M+2*NGz,NGx+1:N+NGx)=BuT
	Bs(M+2*Ngz,NGx+1:N+NGx)=BsT
enddo

end subroutine UpperBoundary


subroutine Lowerboundary(U,W,Bu,Bs,UB,BuB,BsB,hx,hz,M,N)

implicit none

real*8, dimension(:,:) :: U,W,Bu,Bs
real*8, dimension(N) :: UB,BuB,BsB
real*8 :: hx,hz
integer :: M,N
integer :: i

do i=1,N
	U(2,NGx+1:N+NGx)=UB
	W(2,NGx+1:NGx+N)=0
	Bu(2,NGx+1:N+NGx)=Bu(3,NGx+1:NGx+N)
	Bs(2,NGx+1:N+NGx)=Bs(3,NGx+1:NGx+N)
enddo

do i=1,N
	U(1,NGx+1:N+NGx)=UB
	W(1,NGx+1:NGx+N)=-W(NGz+1,NGx+1:NGx+N)
	Bu(1,NGx+1:N+NGx)=Bu(4,NGx+1:NGx+N)
	Bs(1,NGx+1:N+NGx)=Bs(4,NGx+1:NGx+N)
enddo

end subroutine LowerBoundary


subroutine PeriodicBoundary(U,W,Bu,Bs,P,M,N)

!arguments have the appropriate layer of ghost cells surrounding it

implicit none

real*8, dimension(M+4,N+4) :: U,Bu,Bs
real*8, dimension(:,:) :: W
real*8, dimension(:,:) :: P
integer :: M,N
integer :: i

do i=1,NGx
	U(:,i)=U(:,N+i)
	U(:,N+NGx+i)=U(:,NGx+i)
	W(:,i)=W(:,N+i)
	W(:,N+NGx+i)=W(:,NGx+i)
	Bu(:,i)=Bu(:,N+i)
	Bu(:,N+NGx+i)=Bu(:,NGx+i)
	Bs(:,i)=Bs(:,N+i)
	Bs(:,N+NGx+i)=Bs(:,NGx+i)
enddo

P(:,1)=P(:,N+1)
P(:,N+2)=P(:,2)

end subroutine PeriodicBoundary

subroutine PredictBoundaryVel(PUB,PWB,U,W,B,P,N,hx,hz)

implicit none

real*8, dimension(N) :: PUB,PWB
real*8, dimension(3,N) :: U,W
real*8, dimension(2,N) :: B,P
real*8, intent(in) :: hx,hz
integer, intent(in) :: N
real*8, dimension(N) :: Uxx,Wxx,Uzz,Wzz,UB,WB,BB,PxB,PzB,UxB,WxB,UzB,WzB
integer :: i

UxB(1)=(0.5/hx)*(U(2,2)+U(1,2)-U(2,N)-U(1,N))
WxB(1)=(0.5/hx)*(W(2,2)+W(1,2)-W(2,N)-W(1,N))
PxB(1)=(0.5/hx)*(P(2,2)+P(1,2)-P(2,N)-P(1,N))
Uxx(1)=(0.5/hx**2)*(U(2,2)+U(2,2)-2.0*(U(2,1)+U(1,1))+U(2,N)+U(1,N))
Wxx(1)=(0.5/hx**2)*(W(2,2)+W(2,2)-2.0*(W(2,1)+W(1,1))+W(2,N)+W(1,N))

do i=2,N-1
	UxB(i)=(0.5/hx)*(U(2,i+1)+U(1,i+1)-U(2,i-1)-U(1,i-1))
	WxB(i)=(0.5/hx)*(W(2,i+1)+W(1,i+1)-W(2,i-1)-W(1,i-1))
	PxB(i)=(0.5/hx)*(P(2,i+1)+P(1,i)-P(2,i-1)-P(1,i-1))
	Uxx(i)=(0.5/hx**2)*(U(2,i+1)+U(2,i+1)-2.0*(U(2,i)+U(1,i))+U(2,i-1)+U(1,i-1))
	Wxx(i)=(0.5/hx**2)*(W(2,i+1)+W(2,i+1)-2.0*(W(2,i)+W(1,i))+W(2,i-1)+W(1,i-1))

enddo

UxB(N)=(0.5/hx)*(U(2,1)+U(1,1)-U(2,N-1)-U(1,N-1))
WxB(N)=(0.5/hx)*(W(2,1)+W(1,1)-W(2,N-1)-W(1,N-1))
PxB(N)=(0.5/hx)*(P(2,1)+P(1,1)-P(2,N-1)-P(1,N-1))
Uxx(N)=(0.5/hx**2)*(U(2,2)+U(2,2)-2.0*(U(2,N)+U(1,N))+U(2,N-1)+U(1,N-1))
Wxx(N)=(0.5/hx**2)*(W(2,2)+W(2,2)-2.0*(W(2,N)+W(1,N))+W(2,N-1)+W(1,N-1))


Uzz=(1.0/hz**2)*(U(3,:)-2.0*U(2,:)+U(1,:))
Wzz=(1.0/hz**2)*(W(3,:)-2.0*W(2,:)+W(1,:))

UzB=(1.0/hz)*(U(2,:)-U(1,:))
WzB=(1.0/hz)*(W(2,:)-W(1,:))
PzB=(1.0/hz)*(P(2,:)-P(1,:))

UB=0.5*(U(2,:)+U(1,:))
WB=0.5*(W(2,:)+W(1,:))
BB=0.5*(B(1,:)+B(2,:))


do i=1,N
	PUB(i)=UB(i)-dt*UB(i)*UxB(i)-dt*WB(i)*UzB(i)-dt*ev*Uzz(i)-dt*ev*Uxx(i)
	PWB(i)=WB(i)-dt*UB(i)*WxB(i)-dt*WB(i)*WzB(i)-dt*ev*Wzz(i)-dt*ev*Wxx(i)
enddo


end subroutine PredictBoundaryVel


subroutine PredictBoundaryB

implicit none



end subroutine PredictBoundaryB

subroutine HomogeneousBoundary(P,BcT,BcB,M,N)

!Sets up appropriate ghost values of pressure for use in the poisson equation
!P has one layer of ghost cells surrounding all sides; these are the onle things that are changed.

implicit none
real*8, dimension(M+2,N+2) :: P
integer :: M,N
integer :: BcT,BcB

if(BcT==1)then
	P(M+2,2:N+1)=0
elseif(BcT==2)then
	P(M+2,2:N+1)=P(M+1,2:N+1)
endif

if(BcB==1)then
	P(1,2:N+1)=0
elseif(BcB==2)then
	P(1,2:N+1)=P(2,2:N+1)
endif

!periodic lateral boundaries
P(:,1)=P(:,N+1)
P(:,N+2)=P(:,2)

end subroutine HomogeneousBoundary

subroutine NeumannBoundary(Pz,P,N,xc,xf,hx,hz)

implicit none

real*8, dimension(rat*N) :: pz
real*8, dimension(4,N) :: P
integer, intent(in) :: N
real*8, dimension(N) :: xc
real*8, dimension(rat*N) :: xf
real*8, dimension(3,N) :: pz_c
real*8, dimension(N) :: pz_cb
real*8 :: hx,hz
integer :: i

pz_c(1,:)=(1.0/hz)*(P(2,:)-P(1,:))

pz_c(2,:)=(1.0/hz)*(P(3,:)-P(2,:))

pz_c(3,:)=(1.0/hz)*(P(4,:)-P(3,:))

pz_cb=3.0*pz_c(1,:)-3.0*pz_c(2,:)+pz_c(3,:)

pz(1)=quadratic_interp(pz_cb(N),pz_cb(1),pz_cb(2),xc(1)-hx/rat,xc(1),xc(2),xf(1))
pz(2)=quadratic_interp(pz_cb(N),pz_cb(1),pz_cb(2),xc(1)-hx/rat,xc(1),xc(2),xf(2))

do i=2,N-1
	pz(2*i-1)=quadratic_interp(pz_cb(i-1),pz_cb(i),pz_cb(i+1),xc(i-1),xc(i),xc(i+1),xc(i)-hx/rat)
	pz(2*i)=quadratic_interp(pz_cb(i-1),pz_cb(i),pz_cb(i+1),xc(i-1),xc(i),xc(i+1),xc(i)+hx/rat)
enddo

pz(2*N-1)=quadratic_interp(pz_cb(N-1),pz_cb(N),pz_cb(1),xc(N-1),xc(N),xc(N)+hx/rat,xf(2*N-1))
pz(2*N)=quadratic_interp(pz_cb(N-1),pz_cb(N),pz_cb(1),xc(N-1),xc(N),xc(N)+hx/rat,xf(2*N))


pz(1)=0.75*pz_cb(1)+0.25*pz_cb(N)

do i=1,N-1
	pz(2*i)=0.75*pz_cb(i)+0.25*pz_cb(i+1)
	pz(2*i+1)=0.25*pz_cb(i)+0.75*pz_cb(i+1)
enddo

pz(2*N)=0.25*pz_cb(1)+0.75*pz_cb(N)

end subroutine NeumannBoundary

end module Boundary

