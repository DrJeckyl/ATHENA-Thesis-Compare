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
  real, parameter :: pi = atan(1.)*4.
  complex, parameter :: i = (0.,1.)
  real,allocatable, dimension(:) :: bx_raw,by_raw,bz_raw,ux_raw,uy_raw,uz_raw !input from the text files
  real,allocatable :: x_raw(:),y_raw(:),z_raw(:) !x,y,z vectors
  real,allocatable, dimension(:,:,:) :: x, y, z
  real,allocatable, dimension(:,:,:) ::r,kx,ky,kz !Real, input arrays
  complex,allocatable,dimension(:,:,:) :: ibx
  complex, allocatable, dimension(:,:,:) :: kxc,kyc,kzc !complex versions of the k-vector to shift
  complex, allocatable, dimension(:,:,:,:) :: k,u,b,j,E,a,jh
  complex, allocatable, dimension(:,:,:,:) :: u_,b_,j_,E_,a_,jh_
  complex, allocatable, dimension(:,:,:) :: j_dot_b, a_dot_b, phi
  complex, allocatable, dimension(:,:,:) :: j_dot_b_,a_dot_b_,phi_
  complex, allocatable, dimension(:,:) :: u_av, b_av, E_av, j_av, a_av,jh_av
  complex, allocatable, dimension(:) :: phi_av,j_dot_b_av,a_dot_b_av
  complex, allocatable :: divj(:,:,:),divb(:,:,:)
  complex, allocatable :: divj_(:,:,:), divb_(:,:,:)
  complex, allocatable :: grad2(:,:,:)
  complex, allocatable,dimension(:,:,:) :: Gauss, Gauss_, iGauss

  real :: minx, maxx,miny,maxy,minz,maxz
  integer :: NX,NY,NZ
  character *100 buffer
  integer :: ii,jj,kk,counter
  integer :: PGOPEN
  integer :: slice
  real(kind=4), allocatable :: varx(:), vary(:), varz(:), varF(:,:,:)
  real, parameter :: dx = 0.01
  logical, parameter :: shift_both = .false.
  logical, parameter :: shift_k = .true.

  call getarg(1,buffer)
  read(buffer,*) NX
  call getarg(2,buffer)
  read(buffer,*) NY
  call getarg(3,buffer)
  read(buffer,*) NZ

  !NX=32
  !NY=128
  !NZ = NY

  allocate(iGauss(NX,NY,NZ),Gauss_(NX,NY,NZ),Gauss(NX,NY,NZ))
  allocate(ibx(NX,NY,NZ))
  allocate(x_raw(NX),y_raw(NY),z_raw(NZ))
  allocate(x(NX,NY,NZ),y(NX,NY,NZ),z(NX,NY,NZ))
  allocate(bx_raw(NX*NY*NZ),by_raw(NX*NY*NZ),bz_raw(NX*NY*NZ))
  allocate(ux_raw(NX*NY*NZ),uy_raw(NX*NY*NZ),uz_raw(NX*NY*NZ))
  allocate(r(NX,NY,NZ))
  allocate(kx(NX,NY,NZ),ky(NX,NY,NZ),kz(NX,NY,NZ))
  allocate(kxc(NX,NY,NZ),kyc(NX,NY,NZ),kzc(NX,NY,NZ))
  allocate(k(3,NX,NY,NZ),u(3,NX,NY,NZ),b(3,NX,NY,NZ),j(3,NX,NY,NZ),E(3,NX,NY,NZ),a(3,NX,NY,NZ))
  allocate(jh(3,NX,NY,NZ))
  allocate(divj(NX,NY,NZ),divb(NX,NY,NZ))
  allocate(varx(NX),vary(NY),varz(NZ),varF(NX,NY,NZ))
  allocate(j_dot_b(NX,NY,NZ),a_dot_b(NX,NY,NZ),phi(NX,NY,NZ))
  allocate(grad2(NX,NY,NZ))
  allocate(u_(3,NX,NY,NZ),b_(3,NX,NY,NZ),j_(3,NX,NY,NZ),E_(3,NX,NY,NZ),a_(3,NX,NY,NZ))
  allocate(jh_(3,NX,NY,NZ))
  allocate(divj_(NX,NY,NZ),divb_(NX,NY,NZ))
  allocate(j_dot_b_(NX,NY,NZ),a_dot_b_(NX,NY,NZ),phi_(NX,NY,NZ))
  allocate(a_av(3,NZ),b_av(3,NZ),u_av(3,NZ),j_av(3,NZ),E_av(3,NZ))
  allocate(j_dot_b_av(NZ),a_dot_b_av(NZ),phi_av(NZ))

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
  read(16,*) x_raw
  read(17,*) y_raw
  read(18,*) z_raw

  !Reshape the temp array into a 3D 'cube' array
  b(1,:,:,:) = cmplx(reshape(bx_raw,(/NX,Ny,NZ/)),0.)
  b(2,:,:,:) = cmplx(reshape(by_raw,(/NX,Ny,NZ/)),0.)
  b(3,:,:,:) = cmplx(reshape(bz_raw,(/NX,Ny,NZ/)),0.)
  u(1,:,:,:) = cmplx(reshape(ux_raw,(/NX,Ny,NZ/)),0.)
  u(2,:,:,:) = cmplx(reshape(uy_raw,(/NX,Ny,NZ/)),0.)
  u(3,:,:,:) = cmplx(reshape(uz_raw,(/NX,Ny,NZ/)),0.)

  !Ok now everything is put into a big matrix
  !Want to Fourier Transform these arrays
  do ii=1,3
     call fft_3d(b(ii,:,:,:),b_(ii,:,:,:),shift_both)
     call fft_3d(u(ii,:,:,:),u_(ii,:,:,:),shift_both)
  end do
 

  !Define the x and k-vectors

  maxx = maxval(x_raw)
  minx = minval(x_raw)
  maxy = maxval(y_raw)
  miny = minval(y_raw)
  maxz = maxval(z_raw)
  minz = minval(z_raw)

  do ii=1,NX
     x(ii,:,:) = x_raw(ii) + abs(minx)! - (maxx-minx)/2.
     k(1,ii,:,:) = real(ii-NX/2-1)*2.*pi/(maxx+abs(minx))
  end do
  do jj=1,NY
     y(:,jj,:) = y_raw(jj)+ abs(miny)!- (maxy-miny)/2.
     k(2,:,jj,:) = real(jj-NY/2-1)*2.*pi/(maxy+abs(miny))
  end do
  do kk=1,NZ
     z(:,:,kk) = z_raw(kk)+ abs(minz)!- (maxz-minz)/2.
     k(3,:,:,kk) = real(kk-NZ/2-1)*2.*pi/(maxz+abs(minz))
  enddo

!forall(ii=1:NX) x(ii,:,:) = x_raw(ii)
!forall(jj=1:NY) x(:,jj,:) = x_raw(jj)
!forall(kk=1:NZ) x(:,:,kk) = x_raw(kk)

  if(shift_k)then
     do ii=1,3
        call ifftshift3d(k(ii,:,:,:))
     end do
  end if

!!$  print *, 'x(1:5), x(n-5:n)'
!!$  print *, x(1:5,1,1), x(nx-5:nx,1,1)
!!$  print *, 'kx(1:5), kx(n-5:n)'
!!$  print *, k(1,1:5,1,1), k(1,nx-5:nx,1,1)

  !Calculate in k-space
  !j = curl(B)
  !j = k x B
  !j_ = (i*k .x. b_)
  j_ = cross(i*k,b_)
  do ii=1,3
     j_av(ii,:) = average(j(ii,:,:,:),3)
  end do
  !As per Ethan's suggestion, I should take the divergence of this quantitiy to see that I'm doing the curl correctly

  !Divergence is:
  ! div(j) = ik dot j
  !divj = i*kx*jx + i*ky*jy + i*kz*jz
  !divj_ = (i*k .dot. j_)
  divj_ = dot(i*k,j_)
  
  !Calculate div(b) too
  !divb_ = i*k .dot. b_
  divb_ = dot(i*k,b_)

  !Take the average of b and u
  do ii=1,3
     b_av(ii,:) = average(b_(ii,:,:,:),3)
     u_av(ii,:) = average(u_(ii,:,:,:),3)
  end do

  !Ok, presuming that is working, lets do j dot b
  !j_dot_b = jx*bx_ + jy*by_ + jz*bz_
  !j_dot_b = dot(j,b)

  !j_dot_b_ = average3((j_.dot.b_)) - (average4(j_).dot.average4(b_))
  !j_dot_b_ = average3(dot(j_,b_)) - dot(average4(j_),average4(b_))
  j_dot_b_ = j_ .dot. b_
  j_dot_b_av = average(j_dot_b_,3) - (j_av .dot. b_av ) 
  
  !And v x b
  !E_ = average4(u_ .x. b_)
  !E_ = average4(cross(u_,b_)) - cross(average4(u_),average4(b_))
  E_ = u_ .x. b_
  do ii=1,3
     E_av(ii,:) = average( E_(ii,:,:,:), 3 )
  end do
  E_av = E_av - (u_av .x. b_av) 


!!$ print *, '\bar{E}'
!!$  print *, average(E_(1,:,:,:),3), average(E_(2,:,:,:),3), average(E(3,:,:,:),3)
!!$  print *, '\bar{u}'
!!$  print *, average(u_(1,:,:,:),3), average(u_(2,:,:,:),3), average(u_(3,:,:,:),3)
!!$  print *, '\bar{b}'
!!$  print *, average(b_(1,:,:,:),3), average(b_(2,:,:,:),3), average(b_(3,:,:,:),3)
!!$    

  !Calculate the vector potential
  !grad^2 a = -j
  !a = -j/grad^2
  grad2 = -dot(k, k)
  !Add the tiniest number to grad2 where it is zero so that we can avoid the divide by zero later
  do ii=1,nx
     do jj=1,ny
        do kk=1,nz
           if( grad2(ii,jj,kk) == 0. )then
              grad2(ii,jj,kk) = grad2(ii,jj,kk) + tiny(dx)
           end if
        end do
     end do
  end do
  forall(ii=1:3) a_(ii,:,:,:) = -j_(ii,:,:,:)/grad2
  do ii=1,3 
     a_av(ii,:) = average(a_(ii,:,:,:),3)
  end do

  !a dot b
  !a_dot_b_ = average3((a_.dot.b_)) - (average4(a_).dot.average4(b_))
  !a_dot_b_ = average3(dot(a_,b_)) - dot(average4(a_),average4(b_))
  a_dot_b_ = a_ .dot. b_
  a_dot_b_av = average(a_dot_b_,3) - ( a_av .dot. b_av )

  !scalar potential, phi
  !phi_ = -(i*k .dot. (u_ .x. b_))/grad2
  phi_ = -(i*k .dot. (u_.x.b_))/grad2
  phi_av = average(phi_,3)

  !A x ( v x B + grad(phi)) - \bar{A} x (\bar{v x B} + grad{\bar{phi}})
  !jh_ = average4(((u_ .x. b_) + grad(i*k,phi_))) -  (average4(a_) .x. (average4(u_ .x. b_) + grad(i*k,average3(phi_))))

  !jh_ = average4( cross(a_,(cross(u_,b_) + grad(i*k,phi_))) )
  !jh_ = jh_ - cross(average4(a_), (average4(cross(u_,b_) + grad(i*k,average3(phi_)))))
  jh_ = a_ .x. ( E_ + (i*k .V. phi_))
  do ii=1,3
     jh_av(ii,:) = average( jh(ii,:,:,:),3 )
  end do
  jh_av = jh_av - (a_av .x. ( E_av + (i*k(:,1,1,:) .V. phi_av) ) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !Inverse Fourier Transform

  do ii=1,3
     call ifft_3d(j_(ii,:,:,:),j(ii,:,:,:),shift_both)
     call ifft_3d(a_(ii,:,:,:),a(ii,:,:,:),shift_both)
     call ifft_3d(E_(ii,:,:,:),E(ii,:,:,:),shift_both)
     call ifft_3d(jh_(ii,:,:,:),jh(ii,:,:,:),shift_both)
  end do
  call ifft_3d(j_dot_b_,j_dot_b,shift_both)
  call ifft_3d(a_dot_b_,a_dot_b,shift_both)
  call ifft_3d(phi_,phi,shift_both)
  call ifft_3d(b_(1,:,:,:),ibx,shift_both)
  call ifft_3d(divb_,divb,shift_both)
  call ifft_3d(divj_,divj,shift_both)

  !Now plot a slice and check that it's zero
  !Book keeping for the PGPLOT routines(required)
  call PGSLCT(PGOPEN('/xwin'))
!!$call PGASK(.false.)
!!$
  !varx = real(x(:,1,1),kind=4)
  !vary = real(y(1,:,1),kind=4)
  !varz = real(z(1,1,:),kind=4)

  varx = real(x_raw,kind=4)
  vary = real(y_raw,kind=4)
  varz = real(z_raw,kind=4)

  slice = NZ/2
  varF = real(b(1,:,:,:),kind=4)
  open(unit=1,file='bx_input.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(varx(ii)), real(vary(jj)), varF(ii,jj,slice)
        !write(1,*) real(x(ii,1,1)),real(y(1,jj,1)),varF(ii,jj,slice)
     end do
  end do
  close(1) 

  varF = real(b_(1,:,:,:),kind=4)
  open(unit=1,file='bx_fftd.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(k(1,ii,1,1)), real(k(2,1,jj,1)), varF(ii,jj,slice)
        !write(1,*) real(x(ii,1,1)),real(y(1,jj,1)),varF(ii,jj,slice)
     end do
  end do
  close(1)

  varF = real(ibx,kind=4)
  open(unit=1,file='bx_output.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(varx(ii)), real(vary(jj)), varF(ii,jj,slice)
        !write(1,*) real(x(ii,1,1)),real(y(1,jj,1)),varF(ii,jj,slice)
     end do
  end do
  close(1)

  varF = abs(j(1,:,:,:))
  open(unit=1,file='jx_slice.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(varx(ii)), real(vary(jj)), varF(ii,jj,slice)
        !write(1,*) real(x(ii,1,1)),real(y(1,jj,1)),varF(ii,jj,slice)
     end do
  end do
  close(1)

  varF = abs(a(1,:,:,:))
  open(unit=1,file='ax_slice.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(varx(ii)), real(vary(jj)), varF(ii,jj,slice)
     end do
  end do
  close(1)

  varF = abs(E(1,:,:,:))
  open(unit=1,file='Ex_slice.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(varx(ii)), real(vary(jj)), varF(ii,jj,slice)
     end do
  end do
  close(1)

  call plot(varz,real(average(divb,3),kind=4),'z','\bar{Ex}','\bar{Ex} vs. z')

  varF = abs(jh(1,:,:,:))
  open(unit=1,file='jhx_slice.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(varx(ii)), real(vary(jj)), varF(ii,jj,slice)
     end do
  end do
  close(1)

  varF = abs(j_dot_b)
  open(unit=1,file='j_dot_b_slice.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(varx(ii)), real(vary(jj)), varF(ii,jj,slice)
     end do
  end do
  close(1)

  varF = abs(a_dot_b)
  open(unit=1,file='a_dot_b_slice.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(varx(ii)), real(vary(jj)), varF(ii,jj,slice)
     end do
  end do
  close(1)

  varF = abs(phi)
  open(unit=1,file='phi_slice.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(varx(ii)), real(vary(jj)), varF(ii,jj,slice)
     end do
  end do
  close(1)

  varF = abs(divj)
  open(unit=1,file='divj_slice.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(varx(ii)), real(vary(jj)), varF(ii,jj,slice)
     end do
  end do
  close(1)

  varF = abs(divb)
  open(unit=1,file='divb_slice.dat',form='formatted')
  do jj=1,NY
     do ii=1,NX
        write(1,*) real(varx(ii)), real(vary(jj)), varF(ii,jj,slice)
     end do
  end do
  close(1)

call PGCLOS
End Program Main
