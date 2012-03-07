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
real,allocatable, dimension(:) :: bx_raw,by_raw,bz_raw,ux_raw,uy_raw,uz_raw !input from the text files
real,allocatable :: x(:),y(:),z(:) !x,y,z vectors
real,allocatable, dimension(:,:,:) :: bx,by,bz,ux,uy,uz,r,kx,ky,kz !Real, input arrays
complex,allocatable, dimension(:,:,:) :: bx_,by_,bz_,ux_,uy_,uz_ !Fourier Transformed Variables
complex, allocatable, dimension(:,:,:,:) :: k,u,b,j,E,a,jh
complex, allocatable, dimension(:,:,:) :: j_dot_b, a_dot_b, phi
complex, allocatable :: divj(:,:,:)
complex, allocatable :: grad2(:,:,:)
integer :: NX,NY,NZ
character *100 buffer
integer :: ii,jj,kk,counter
integer :: PGOPEN
real(kind=4), allocatable :: varx(:), vary(:), varz(:), varF(:,:,:)

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
allocate(r(NX,NY,NZ))
allocate(kx(NX,NY,NZ),ky(NX,NY,NZ),kz(NX,NY,NZ))
allocate(k(3,NX,NY,NZ),u(3,NX,NY,NZ),b(3,NX,NY,NZ),j(3,NX,NY,NZ),E(3,NX,NY,NZ),a(3,NX,NY,NZ))
allocate(jh(3,NX,NY,NZ))
allocate(divj(NX,NY,NZ))
allocate(varx(NX),vary(NY),varz(NZ),varF(NX,NY,NZ))
allocate(j_dot_b(NX,NY,NZ),a_dot_b(NX,NY,NZ),phi(NX,NY,NZ))
allocate(grad2(NX,NY,NZ))

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
bx_=cmplx(bx,0.)
by_=cmplx(by,0.)
bz_=cmplx(bz,0.)
ux_=cmplx(ux,0.)
uy_=cmplx(uy,0.)
uz_=cmplx(uz,0.)
call fft3d(bx_,bx_)
call fft3d(by_,by_)
call fft3d(bz_,bz_)
call fft3d(ux_,ux_)
call fft3d(uy_,uy_)
call fft3d(uz_,uz_)
b(1,:,:,:) = bx_
b(2,:,:,:) = by_
b(3,:,:,:) = bz_
u(1,:,:,:) = ux_
u(2,:,:,:) = uy_
u(3,:,:,:) = uz_

!Define the k-vector
do ii=1,NX
   kx(ii,:,:) = x(ii)-x(ii)/2.
end do
do jj=1,NY
   ky(:,jj,:) = y(jj)-y(jj)/2.
end do
do kk=1,NZ
   kz(:,:,kk) = z(kk)-z(kk)/2.
end do
k(1,:,:,:) = kx
k(2,:,:,:) = ky
k(3,:,:,:) = kz
!Calculate in k-space
!j = curl(B)
!j = k x B

!Do a cross product
!Writing it in the main program to test with the intention to write a general function for it later
!Maybe have it work for 1,2,3 dimensions?
!Defined as:
!C_x = A_y*B_z - A_z*B_y
!C_y = A_z*B_x - A_x*B_z
!C_z = A_x*B_y - A_y*B_x
!jx= i*ky*bz_ - i*kz*by_
!jy= i*kz*bx_ - i*kx*bz_
!jz= i*kx*by_ - i*ky*bx_
j = cross(i*k,b)

!As per Ethan's suggestion, I should take the divergence of this quantitiy to see that I'm doing the curl correctly

!Divergence is:
! div(j) = ik dot j
!divj = i*kx*jx + i*ky*jy + i*kz*jz
divj = dot(i*k,j)

!Ok, presuming that is working, lets do j dot b
!j_dot_b = jx*bx_ + jy*by_ + jz*bz_
j_dot_b = dot(j,b)

!And v x b
E = cross(u,b)

!Calculate the vector potential
!grad^2 a = -j
!a = -j/grad^2
grad2 = dot(k,k)
forall(ii=1:3) a(ii,:,:,:) = -j(ii,:,:,:)/grad2

!a dot b
a_dot_b = dot(a,b)

!scalar potential, phi
phi = dot(i*k,E)/grad2

!A x ( v x B + grad(phi)) - \bar{A} x (\bar{v x B} + grad{\bar{phi}})
jh = cross(a,(E+ grad(i*k,phi)))


!Now plot a slice and check that it's zero
!Book keeping for the PGPLOT routines(required)
!!$call PGSLCT(PGOPEN('/ps'))
!!$call PGASK(.false.)
!!$
!!$varx = real(kx(:,1,1),kind=4)
!!$vary = real(ky(1,:,1),kind=4)
!!$varz = real(kz(1,1,:),kind=4)
varF = abs(divj)
!!$

open(unit=1,file='divj.dat',form='formatted')
do jj=1,NY
   do ii=1,NX
      write(1,*) real(kx(ii,1,1)), real(ky(1,jj,1)), varF(ii,jj,NZ/2)
   end do
end do
close(1)


!Inverse Fourier Transform
End Program Main
