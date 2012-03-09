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
complex(kind=8), parameter :: i = (0._8,1._8)
real(kind=8),allocatable, dimension(:) :: bx_raw,by_raw,bz_raw,ux_raw,uy_raw,uz_raw !input from the text files
real(kind=8),allocatable :: x(:),y(:),z(:) !x,y,z vectors
real(kind=8),allocatable, dimension(:,:,:) :: bx,by,bz,ux,uy,uz,r,kx,ky,kz !Real, input arrays
complex(kind=8),allocatable, dimension(:,:,:) :: bx_,by_,bz_,ux_,uy_,uz_ !Fourier Transformed Variables
complex(kind=8), allocatable, dimension(:,:,:,:) :: k,u,b,j,E,a,jh
complex(kind=8), allocatable, dimension(:,:,:) :: j_dot_b, a_dot_b, phi
complex(kind=8), allocatable :: divj(:,:,:)
complex(kind=8), allocatable :: grad2(:,:,:)
integer :: NX,NY,NZ
character *100 buffer
integer :: ii,jj,kk,counter
integer :: PGOPEN
integer :: slice
real(kind=4), allocatable :: varx(:), vary(:), varz(:), varF(:,:,:)
real(kind=8), parameter :: dx = 0.01_8
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
bx_=cmplx(bx,0._8)
by_=cmplx(by,0._8)
bz_=cmplx(bz,0._8)
ux_=cmplx(ux,0._8)
uy_=cmplx(uy,0._8)
uz_=cmplx(uz,0._8)
call fft3d(bx_,bx_,.true.)
call fft3d(by_,by_,.true.)
call fft3d(bz_,bz_,.true.)
call fft3d(ux_,ux_,.true.)
call fft3d(uy_,uy_,.true.)
call fft3d(uz_,uz_,.true.)
b(1,:,:,:) = bx_
b(2,:,:,:) = by_
b(3,:,:,:) = bz_
u(1,:,:,:) = ux_
u(2,:,:,:) = uy_
u(3,:,:,:) = uz_

!Define the k-vector
do ii=1,NX
   kx(ii,:,:) = x(ii)-x(ii)/2._8
end do
do jj=1,NY
   ky(:,jj,:) = y(jj)-y(jj)/2._8
end do
do kk=1,NZ
   kz(:,:,kk) = z(kk)-z(kk)/2._8
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
!j_dot_b = dot(j,b)
j_dot_b = average3(dot(j,b)) - dot( average4(j), average4(b) )

!And v x b
E = cross(u,b)
E = average4(E) - cross( average4(u), average4(b) )

!Calculate the vector potential
!grad^2 a = -j
!a = -j/grad^2
grad2 = dot(k,k)
!Add the tiniest number to grad2 where it is zero so that we can avoid the divide by zero later
do ii=1,nx
   do jj=1,ny
      do kk=1,nz
         if( grad2(ii,jj,kk) == 0._8 )then
            grad2(ii,jj,kk) = grad2(ii,jj,kk) + tiny(dx)
         end if
      end do
   end do
end do
forall(ii=1:3) a(ii,:,:,:) = -j(ii,:,:,:)/grad2

!a dot b
a_dot_b = dot(a,b)
a_dot_b = average3(a_dot_b) - dot( average4(a), average4(b) )

!scalar potential, phi
phi = dot(i*k,E)/grad2

!A x ( v x B + grad(phi)) - \bar{A} x (\bar{v x B} + grad{\bar{phi}})
jh = cross(a,(E+ grad(i*k,phi)))
jh = average4(jh) - cross( average4(a), (average4(E) + grad(k,average3(phi))))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!Inverse Fourier Transform

do ii=1,3
   call ifft3d(j(i,:,:,:),j(i,:,:,:),.true.)
   call ifft3d(a(i,:,:,:),a(i,:,:,:),.true.)
   call ifft3d(E(i,:,:,:),E(i,:,:,:),.true.)
   call ifft3d(jh(i,:,:,:),jh(i,:,:,:),.true.)
end do
call ifft3d(j_dot_b,j_dot_b,.true.)
call ifft3d(a_dot_b,a_dot_b,.true.)
call ifft3d(phi,phi,.true.)

!Now plot a slice and check that it's zero
!Book keeping for the PGPLOT routines(required)
!!$call PGSLCT(PGOPEN('/ps'))
!!$call PGASK(.false.)
!!$
!!$varx = real(kx(:,1,1),kind=4)
!!$vary = real(ky(1,:,1),kind=4)
!!$varz = real(kz(1,1,:),kind=4)
varF = abs(j(1,:,:,:))
slice = NZ/2
open(unit=1,file='jx_slice.dat',form='formatted')
do jj=1,NY
   do ii=1,NX
      write(1,*) real(kx(ii,1,slice)), real(ky(1,jj,slice)), varF(ii,jj,slice)
   end do
end do
close(1)

varF = abs(a(1,:,:,:))
open(unit=1,file='ax_slice.dat',form='formatted')
do jj=1,NY
   do ii=1,NX
      write(1,*) real(kx(ii,1,slice)), real(ky(1,jj,slice)), varF(ii,jj,slice)
   end do
end do
close(1)

varF = abs(E(1,:,:,:))
open(unit=1,file='Ex_slice.dat',form='formatted')
do jj=1,NY
   do ii=1,NX
      write(1,*) real(kx(ii,1,slice)), real(ky(1,jj,slice)), varF(ii,jj,slice)
   end do
end do
close(1)

varF = abs(jh(1,:,:,:))
open(unit=1,file='jhx_slice.dat',form='formatted')
do jj=1,NY
   do ii=1,NX
      write(1,*) real(kx(ii,1,slice)), real(ky(1,jj,slice)), varF(ii,jj,slice)
   end do
end do
close(1)

varF = abs(j_dot_b)
open(unit=1,file='j_dot_b_slice.dat',form='formatted')
do jj=1,NY
   do ii=1,NX
      write(1,*) real(kx(ii,1,slice)), real(ky(1,jj,slice)), varF(ii,jj,slice)
   end do
end do
close(1)

varF = abs(a_dot_b)
open(unit=1,file='a_dot_b_slice.dat',form='formatted')
do jj=1,NY
   do ii=1,NX
      write(1,*) real(kx(ii,1,slice)), real(ky(1,jj,slice)), varF(ii,jj,slice)
   end do
end do
close(1)

varF = abs(phi)
open(unit=1,file='phi_slice.dat',form='formatted')
do jj=1,NY
   do ii=1,NX
      write(1,*) real(kx(ii,1,slice)), real(ky(1,jj,slice)), varF(ii,jj,slice)
   end do
end do
close(1)
End Program Main
