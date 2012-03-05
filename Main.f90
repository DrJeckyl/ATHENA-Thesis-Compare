!Program for taking in simulation data that has already been broken down into a column format text file

!CALL with 3 arguments: NX,NY,NZ

!The goal is to calculate:
!j . b = \bar{J . B} - \bar{J} . \bar{B}
!v x b = \bar{v x B} - \bar{v} x \bar{B}
!a . b = \bar{a . b} - \bar{a} . \bar{b}
!A x ( v x B + grad(phi)) - \bar{A} x (\bar{v x B} + grad{\bar{phi}})
!grad^2 phi = div (v x b)
!grad^2 a = -curl(B)
Program Main
use Funcs
IMPLICIT NONE
complex, parameter :: i = (0.,1.)
real,allocatable, dimension(:) :: bx_raw,by_raw,bz_raw,ux_raw,uy_raw,uz_raw
real,allocatable :: x(:),y(:),z(:)
real,allocatable, dimension(:,:,:) :: bx,by,bz,ux,uy,uz,r,k
complex,allocatable, dimension(:,:,:) :: bx_,by_,bz_,ux_,uy_,uz_
complex, allocatable, dimension(:,:,:) :: bx2_,bxr_,bx2r_
real, dimension(32,128,128) :: ibxr_,ibx2r_
complex, dimension(32,128,128) :: ibx_,ibx2_
integer :: NX,NY,NZ
character *100 buffer
integer :: ii,jj,kk

call getarg(1,buffer)
read(buffer,*) NX
call getarg(2,buffer)
read(buffer,*) NY
call getarg(3,buffer)
read(buffer,*) NZ

!NX=32
!NY=128
!NZ = NY

allocate(x(NX),y(NY),z(NZ))
allocate(bx_raw(NX*NY*NZ),by_raw(NX*NY*NZ),bz_raw(NX*NY*NZ))
allocate(ux_raw(NX*NY*NZ),uy_raw(NX*NY*NZ),uz_raw(NX*NY*NZ))
allocate(bx(NX,NY,NZ),by(NX,NY,NZ),bz(NX,NY,NZ))
allocate(ux(NX,NY,NZ),uy(NX,NY,NZ),uz(NX,NY,NZ))
allocate(bx_(NX,NY,NZ),by_(NX,NY,NZ),bz_(NX,NY,NZ))
allocate(ux_(NX,NY,NZ),uy_(NX,NY,NZ),uz_(NX,NY,NZ))
allocate(r(NX,NY,NZ),k(NX,NY,NZ))
allocate(bx2_(NX,NY,NZ),bxr_(NX,NY,NZ),bx2r_(NX,NY,NZ))

!Need to read in the outputs from simulation
!open the files
open(unit=10,file='../AthenaDumpsFromShane/0500.small.formatted/bx.dat')
open(unit=11,file='../AthenaDumpsFromShane/0500.small.formatted/by.dat')
open(unit=12,file='../AthenaDumpsFromShane/0500.small.formatted/bz.dat')
open(unit=13,file='../AthenaDumpsFromShane/0500.small.formatted/vx.dat')
open(unit=14,file='../AthenaDumpsFromShane/0500.small.formatted/vy.dat')
open(unit=15,file='../AthenaDumpsFromShane/0500.small.formatted/vz.dat')
open(unit=16,file='../AthenaDumpsFromShane/0500.small.formatted/x.dat')
open(unit=17,file='../AthenaDumpsFromShane/0500.small.formatted/y.dat')
open(unit=18,file='../AthenaDumpsFromShane/0500.small.formatted/z.dat')

!Read the data into a temporary array
read(10,*) bx_raw
read(11,*) by_raw
read(12,*) bz_raw
read(13,*) ux_raw
read(14,*) uy_raw
read(15,*) uz_raw
read(16,*) x
read(17,*) y
read(18,*) z

!Reshape the temp array into a 3D 'cube' array
bx = reshape(bx_raw,(/NX,Ny,NZ/))
by = reshape(by_raw,(/NX,Ny,NZ/))
bz = reshape(bz_raw,(/NX,Ny,NZ/))
ux = reshape(ux_raw,(/NX,Ny,NZ/))
uy = reshape(uy_raw,(/NX,Ny,NZ/))
uz = reshape(uz_raw,(/NX,Ny,NZ/))

!Ok now everything is put into a big matrix
!Want to Fourier Transform these arrays
bx_=cmplx(bx)
call fft_real_3d(bx,bxr_)
call fft_3d(cmplx(bx),bx_)
!!$call fft_3d(cmplx(by),by_)
!!$call fft_3d(cmplx(bz),bz_)
!!$call fft_3d(cmplx(ux),ux_)
!!$call fft_3d(cmplx(uy),uy_)
!!$call fft_3d(cmplx(uz),uz_)

!!$!Testing the multi-d fft
!!$!Using alex's method
!!$bx2_ = cmplx(bx)
!!$do kk=1,NZ
!!$   do jj=1,NY
!!$      xt = cmplx(bx(:,jj,kk))
!!$      call fft(xt,xt_)
!!$      bx2_(:,jj,kk) = xt_
!!$   end do
!!$   do iii=1,NX
!!$      yt = cmplx(bx(iii,:,kk))
!!$      call fft(yt,yt_)
!!$      bx2_(iii,:,kk) = yt_
!!$   end do
!!$end do
!!$do jj=1,NY
!!$   do iii=1,NX
!!$      zt = cmplx(bx(iii,jj,:))
!!$      call fft(zt,zt_)
!!$      bx2_(iii,jj,:) = zt_
!!$   end do
!!$end do
!!$
bx2_=cmplx(bx)
do kk=1,NZ
   do jj=1,NY
      call fft(cmplx(bx(:,jj,kk)),bx2_(:,jj,kk))
   end do
end do
do kk=1,NZ
   do ii=1,NX
      call fft(cmplx(bx(ii,:,kk)),bx2_(ii,:,kk))
   end do
end do
do jj=1,NY
   do ii=1,NX
      call fft(cmplx(bx(ii,jj,:)),bx2_(ii,jj,:))
   end do
end do

!Testing the new real to complex fftw's
do kk=1,NZ
   do jj=1,NY
      call fft_real_1d(bx(:,jj,kk),bx2r_(:,jj,kk))
   end do
end do
do kk=1,NZ
   do ii=1,NX
      call fft_real_1d(bx(ii,:,kk),bx2r_(ii,:,kk))
   end do
end do
do jj=1,NY
   do ii=1,NX
      call fft_real_1d(bx(ii,jj,:),bx2r_(ii,jj,:))
   end do
end do
!!!!

!Do the inverse transforms
call ifft_3d(bx_,ibx_)
call ifft_real_3d(bxr_,ibxr_)
do kk=1,NZ
   do jj=1,NY
      call ifft(bx_(:,jj,kk),ibx2_(:,jj,kk))
   end do
end do
do kk=1,NZ
   do ii=1,NX
      call ifft(bx_(ii,:,kk),ibx2_(ii,:,kk))
   end do
end do
do jj=1,NY
   do ii=1,NX
      call ifft(bx_(ii,jj,:),ibx2_(ii,jj,:))
   end do
end do

do kk=1,NZ
   do jj=1,NY
      call ifft_real_1d(ibx_(:,jj,kk),ibx2r_(:,jj,kk))
   end do
end do
do kk=1,NZ
   do ii=1,NX
      call ifft_real_1d(ibx_(ii,:,kk),ibx2r_(ii,:,kk))
   end do
end do
do jj=1,NY
   do ii=1,NX
      call ifft_real_1d(ibx_(ii,jj,:),ibx2r_(ii,jj,:))
   end do
end do

print *,'bx:'
print *,bx(:,1,1)
print *,'bx_:'
print *,bx_(:,1,1)
print *,'bx2_:'
print *,bx2_(:,1,1)
print *,'bxr_:'
print *,bxr_(:,1,1)
print *,'bx2r_:'
print *,bx2r_(:,1,1)
print *, ''
print *, ''
print *,'ibx_:'
print *,ibx_(:,1,1)
print *,'ibx2_:'
print *,ibx2_(:,1,1)
print *,'ibxr_:'
print *,ibxr_(:,1,1)
print *,'ibx2r_:'
print *,ibx2r_(:,1,1)
print *,'***'
print *,'***'
print *,'bx/ibx_:'
print *,bx(:,1,1)/abs(ibx_(:,1,1))
print *,'bx/ibxr_:'
print *,bx(:,1,1)/abs(ibxr_(:,1,1))
print *,'bx/ibx2_:'
print *,bx(:,1,1)/abs(ibx2_(:,1,1))
print *,'bx/ibx2r_:'
print *,bx(:,1,1)/abs(ibx2r_(:,1,1))




!Calculate in k-space
!j = curl(B)
!j = k x B


!Inverse Fourier Transform
End Program Main
