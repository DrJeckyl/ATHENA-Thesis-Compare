Program Gauss
  use Funcs
  implicit none
  
  integer :: PGOPEN
  real, parameter :: pi=atan(1.)*4.
  complex, parameter :: i=(0.,1.)
  integer, parameter :: nx = 100
  integer, parameter :: ny = 100
  integer, parameter :: nz = 100
  complex, dimension(nx,ny,nz) :: F_, iF_
  complex, dimension(nx,ny,nz) :: F 
  complex, dimension(3,nx,ny,nz) :: b,b_,ib,b_anal,k,a,a_,ib_anal
  real, dimension(nx,ny,nz) :: x,y,z,R
  real, dimension(nx,ny,nz) :: kx,ky,kz
  complex, dimension(nx,ny,nz) :: kxt,kyt,kzt
  integer :: ii,jj,kk
  integer :: mid
  real,parameter :: dx=0.1
  integer, parameter :: slice=nx/2 + 3
  mid = 0!nx/2
 
  do ii=1,nx
     x(ii,:,:) = (real(ii)-1.)*dx - 5.
     y(:,ii,:) = (real(ii)-1.)*dx - 5.
     z(:,:,ii) = (real(ii)-1.)*dx - 5.
  end do
  R = sqrt((x-real(mid)*dx)**2 + (y-real(mid)*dx)**2 + (z-real(mid)*dx)**2)
  F = exp(-R**2)

  !kx = x - real(mid)*dx
  !ky = y - real(mid)*dx
  !kz = z - real(mid)*dx

  forall(ii=1:nx) kx(ii,:,:) = real(ii-nx/2-1)!*2.*pi/10.
  forall(jj=1:ny) ky(:,jj,:) = real(jj-ny/2-1)!*2.*pi/10.
  forall(kk=1:nz) kz(:,:,kk) = real(kk-nz/2-1)!*2.*pi/10.

  kxt = cmplx(kx,0.)
  kyt = cmplx(ky,0.)
  kzt = cmplx(kz,0.)
  call fftshift3d(kxt)
  call fftshift3d(kyt)
  call fftshift3d(kzt)
  kx = real(kxt)
  ky = real(kyt)
  kz = real(kzt)

  print *, 'x(1:5), x(n-5:n)'
  print *, x(1:5,1,1), x(nx-5:nx,1,1)
  print *, 'kx(1:5), kx(n-5:n)'
  print *, kx(1:5,1,1), kx(nx-5:nx,1,1)

  open(unit=1,file='gauss.dat',form='formatted')
  do kk=1,nz
     do jj=1,ny
        write(1,*) x(jj,1,1), y(1,kk,1), real(F(jj,kk,slice))
     end do
  end do
  close(1)
  
  !Now Fourier Transform that Gaussian!
  call fft_3d(F,F_,.false.)
  
  open(unit=1,file='gauss1_.dat',form='formatted')
  do kk=1,nz
     do jj=1,ny
        write(1,*) kx(jj,1,1), ky(1,kk,1), abs(F_(jj,kk,slice))
     end do
  end do
  close(1)

  !For completeness, inverse back
  call ifft_3d(F_,iF_,.false.)
  
  open(unit=1,file='igauss1_.dat',form='formatted')
  do kk=1,nz
     do jj=1,ny
         write(1,*) x(jj,1,1), y(1,kk,1), abs(iF_(jj,kk,slice))
     end do
  end do
  close(1)

  !Set Ax = F and find the magnetic field derived from that
  a(1,:,:,:) = F
  a(2,:,:,:) = (0.,0.)
  a(3,:,:,:) = (0.,0.)

  b(1,:,:,:) = (0.,0.)
  b(2,:,:,:) = -2.*z*a(1,:,:,:)
  b(3,:,:,:) = 2.*y*a(1,:,:,:)

  !Take the fourier transform of the answer to compare to the k-space method
  do ii=1,3
     call fft_3d(b(ii,:,:,:),b_anal(ii,:,:,:),.false.)
  end do

  !k-space method
  !Take the fft of the vector potential
  do ii=1,3
     call fft_3d(a(ii,:,:,:),a_(ii,:,:,:),.false.)
  end do
  !Take the "curl" of a
  k(1,:,:,:) = kx
  k(2,:,:,:) = ky
  k(3,:,:,:) = kz

  b_ = i*k .x. a_

  !Now do an ifft back to real space
  do ii=1,3
     call ifft_3d(b_(ii,:,:,:),ib(ii,:,:,:),.false.)
     call ifft_3d(b_anal(ii,:,:,:),ib_anal(ii,:,:,:),.false.)
  end do
  
  open(unit=1,file='test_by.dat',form='formatted')
  do kk=1,nx
     do jj=1,ny
        write(1,*) x(jj,1,1), y(1,kk,1), abs(b(3,jj,kk,slice))
     end do
  end do
  close(1)

  open(unit=1,file='test_by_fft_analytic.dat',form='formatted')
  do kk=1,nx
     do jj=1,ny
        write(1,*) kx(jj,1,1), ky(1,kk,1), abs(b_anal(3,jj,kk,slice))
     end do
  end do
  close(1)

  open(unit=1,file='test_by_k-space.dat',form='formatted')
  do kk=1,nx
     do jj=1,ny
        write(1,*) kx(jj,1,1), ky(1,kk,1), abs(b_(3,jj,kk,slice))
     end do
  end do
  close(1)

 open(unit=1,file='test_by_ifft.dat',form='formatted')
  do kk=1,nx
     do jj=1,ny
        write(1,*) x(jj,1,1), y(1,kk,1), abs(ib(3,jj,kk,slice))
     end do
  end do
  close(1)

 open(unit=1,file='test_iby_anal.dat',form='formatted')
  do kk=1,nx
     do jj=1,ny
        write(1,*) x(jj,1,1), y(1,kk,1), abs(ib_anal(3,jj,kk,slice))
     end do
  end do
  close(1)

  print *, abs(ib_anal(2,1:10,slice,slice))/abs(ib(2,1:10,slice,slice))
end Program Gauss
