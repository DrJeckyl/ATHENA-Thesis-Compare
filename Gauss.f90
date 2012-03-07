Program Gauss
  use Funcs
  implicit none
  
  integer :: PGOPEN
  complex, parameter :: i=(0.,1.)
  complex, allocatable, dimension(:,:,:) :: F_, iF_
  complex, allocatable, dimension(:,:,:) :: F
  real, allocatable, dimension(:,:,:) :: x,y,z,R
  integer :: ii,jj,kk
  integer :: nx,ny,nz,mid
  real,parameter :: dx=0.1
  real(kind=4), allocatable :: varx(:), vary(:), varz(:), varF(:,:,:)

  nx=100
  ny=nx
  nz=nx
  mid = nx/2
  
  allocate(F(nx,ny,nz),F_(nx,ny,nz),iF_(nx,ny,nz))
  allocate(x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz),R(nx,ny,nz))

  do ii=1,nx
     x(ii,:,:) = real(ii-1.)*dx
     y(:,ii,:) = real(ii-1.)*dx
     z(:,:,ii) = real(ii-1.)*dx
  end do
  R = sqrt((x-mid*dx)**2 + (y-mid*dx)**2 + (z-mid*dx)**2)
  F = exp(-R**2)
  
  open(unit=1,file='gauss.dat',form='formatted')
  do kk=1,nz
     do jj=1,ny
        write(1,*) y(1,jj,1), z(1,1,kk), real(F(1,jj,kk))
     end do
  end do
  close(1)
  
  !Now Fourier Transform that Gaussian!
  call fft3d(F,F_)
  
  open(unit=1,file='gauss1_.dat',form='formatted')
  do kk=1,nz
     do jj=1,ny
        write(1,*) y(1,jj,1), z(1,1,kk), abs(F_(1,jj,kk))
     end do
  end do
  close(1)

  !For completeness, inverse back
  call ifft3d(F_,iF_)
  open(unit=1,file='igauss1_.dat',form='formatted')
  do kk=1,nz
     do jj=1,ny
        write(1,*) y(1,jj,1), z(1,1,kk), abs(iF_(1,jj,kk))
     end do
  end do
  close(1)
  
end Program Gauss
