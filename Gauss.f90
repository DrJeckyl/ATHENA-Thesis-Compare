Program Gauss
  use Funcs
  implicit none
  
  integer :: PGOPEN
  complex(kind=8), parameter :: i=(0._8,1._8)
  complex(kind=8), allocatable, dimension(:,:,:) :: F_, iF_
  complex(kind=8), allocatable, dimension(:,:,:) :: F
  real(kind=8), allocatable, dimension(:,:,:) :: x,y,z,R
  real(kind=8), allocatable, dimension(:,:,:) :: kx,ky,kz
  complex(kind=8), allocatable, dimension(:,:,:) :: kxt,kyt,kzt
  integer :: ii,jj,kk
  integer :: nx,ny,nz,mid
  real(kind=8),parameter :: dx=0.1_8
  real(kind=4), allocatable :: varx(:), vary(:), varz(:), varF(:,:,:)

  nx=100
  ny=nx
  nz=nx
  mid = nx/2
  
  allocate(F(nx,ny,nz),F_(nx,ny,nz),iF_(nx,ny,nz))
  allocate(x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz),R(nx,ny,nz))
  allocate(kx(nx,ny,nz),ky(nx,ny,nz),kz(nx,ny,nz))
  allocate(kxt(nx,ny,nz),kyt(nx,ny,nz),kzt(nx,ny,nz))

  do ii=1,nx
     x(ii,:,:) = (real(ii,kind=8)-1._8)*dx
     y(:,ii,:) = (real(ii,kind=8)-1._8)*dx
     z(:,:,ii) = (real(ii,kind=8)-1._8)*dx
  end do
  R = sqrt((x-real(mid,kind=8)*dx)**2 + (y-real(mid,kind=8)*dx)**2 + (z-real(mid,kind=8)*dx)**2)
  F = exp(-R**2)
  
!!$  do ii=1,nx/2
!!$     kx(ii,:,:) = real((ii-1.)/(nx*1.)/dx)
!!$     kx(nx/2+ii,:,:) = real((-nx+ii-1.+nx/2.)/(nx*1.)/dx)
!!$  end do
!!$  do ii=1,ny/2
!!$     ky(:,ii,:) = real((ii-1.)/(ny*1.)/dx)
!!$     ky(:,ny/2+ii,:) = real((-ny+ii-1.+ny/2.)/(ny*1.)/dx)
!!$  end do
!!$  do ii=1,nz/2
!!$     kz(:,:,ii) = real((ii-1.)/(nz*1.)/dx)
!!$     kz(:,:,nz/2+ii) = real((-nz+ii-1.+nz/2.)/(nz*1.)/dx)
!!$  end do

  kx = x - real(mid,kind=8)*dx
  ky = y - real(mid,kind=8)*dx
  kz = z - real(mid,kind=8)*dx

  kxt = cmplx(kx,0._8,kind=8)
  kyt = cmplx(ky,0._8,kind=8)
  kzt = cmplx(kz,0._8,kind=8)
  call fftshift3d(kxt)
  call fftshift3d(kyt)
  call fftshift3d(kzt)
  kx = real(kxt,kind=8)
  ky = real(kyt,kind=8)
  kz = real(kzt,kind=8)

  print *, 'x(1:5), x(n-5:n)'
  print *, x(1:5,1,1), x(nx-5:nx,1,1)
  print *, 'kx(1:5), kx(n-5:n)'
  print *, kx(1:5,1,1), kx(nx-5:nx,1,1)

  open(unit=1,file='gauss.dat',form='formatted')
  do kk=1,nz
     do jj=1,ny
        write(1,*) x(jj,1,1), y(1,kk,1), real(F(jj,kk,1))
     end do
  end do
  close(1)
  
  !Now Fourier Transform that Gaussian!
  call fft_3d(F,F_,.false.)
  
  open(unit=1,file='gauss1_.dat',form='formatted')
  do kk=1,nz
     do jj=1,ny
        write(1,*) kx(jj,1,1), ky(1,kk,1), abs(F_(jj,kk,1))
     end do
  end do
  close(1)

  !For completeness, inverse back
  call ifft_3d(F_,iF_,.false.)
  
  open(unit=1,file='igauss1_.dat',form='formatted')
  do kk=1,nz
     do jj=1,ny
         write(1,*) x(jj,1,1), y(1,kk,1), abs(iF_(jj,kk,1))
     end do
  end do
  close(1)
  
end Program Gauss
