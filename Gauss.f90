Program Gauss
  use Funcs
  implicit none
  
  integer :: PGOPEN
  complex, parameter :: i=(0.,1.)
  real, allocatable, dimension(:,:,:) :: F
  complex, allocatable, dimension(:,:,:) :: F1_,F2_,F3_,F4_
  real, allocatable, dimension(:,:,:) :: x,y,z,R
  integer :: ii,jj,kk
  integer :: nx,ny,nz,mid
  real,parameter :: dx=0.1
  real(kind=4), allocatable :: varx(:), vary(:), varz(:), varF(:,:,:)
  complex(kind=8) :: test(2,2,2)

  test = reshape((/1,2,3,4,5,6,7,8,9/),(/2,2,2/))
  print *, 'test(:,:,1):'
  print *, real(test(1,1,1)),real(test(2,1,1))
  print *, real(test(1,2,1)),real(test(2,2,1))
  print *, 'test(:,:,2):'
  print *, real(test(1,1,2)),real(test(2,1,2))
  print *, real(test(1,2,2)),real(test(2,2,2))
  call fftshift3d(test)
    print *, 'test(:,:,1):'
  print *, real(test(1,1,1)),real(test(2,1,1))
  print *, real(test(1,2,1)),real(test(2,2,1))
  print *, 'test(:,:,2):'
  print *, real(test(1,1,2)),real(test(2,1,2))
  print *, real(test(1,2,2)),real(test(2,2,2))
  nx=100
  ny=nx
  nz=nx
  mid = nx/2

  allocate(F1_(nx,ny,nz),F2_(nx,ny,nz),F3_(nx,ny,nz),F4_(nx,ny,nz))
  allocate(F(nx,ny,nz),x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz),R(nx,ny,nz))
  allocate(varx(nx),vary(ny),varz(nz),varF(nx,ny,nz))

  do ii=1,nx
     x(ii,:,:) = real(ii-1)*dx
  end do
  do jj=1,ny
     y(:,jj,:) = real(jj-1)*dx
  end do
  do kk=1,nz
     z(:,:,kk) = real(kk-1)*dx
  end do
  R = sqrt((x-mid*dx)**2 + (y-mid*dx)**2 + (z-mid*dx)**2)
  F = exp(-R**2)

  !Book keeping for the PGPLOT routines(required)
  call PGSLCT(PGOPEN('/ps'))
  !call PGASK(.false.)
  varx = real(x(:,1,1),kind=4)
  vary = real(y(1,:,1),kind=4)
  varz = real(z(1,1,:),kind=4)
  varF = real(F,kind=4)

  call plot(varx,varF(mid,mid,:),'x','F(x,1,1)','Gaussian')
  
  !print *, varF(1,1,:)
  
  open(unit=1,file='gauss.dat',form='formatted')
  do kk=1,nz
     do jj=1,ny
        write(1,*) y(1,jj,1), z(1,1,kk), F(mid,jj,kk)
     end do
  end do
  close(1)

  !Fourier Transform...
  !call ifftshift3d(cmplx(F,kind=8))
  call fft_3d(cmplx(F),F1_)
  call fft_real_3d(F,F2_)

  do kk=1,nz
     do jj=1,ny
        call fft(cmplx(F(:,jj,kk)),F3_(:,jj,kk))
     end do
  end do
  do kk=1,nz
     do ii=1,nx
        call fft(cmplx(F(ii,:,kk)),F3_(ii,:,kk))
     end do
  end do
  do jj=1,ny
     do ii=1,nx
        call fft(cmplx(F(ii,jj,:)),F3_(ii,jj,:))
     end do
  end do

  do kk=1,nz
     do jj=1,ny
        call fft_real_1d((F(:,jj,kk)),F4_(:,jj,kk))
     end do
  end do
  do kk=1,nz
     do ii=1,nx
        call fft_real_1d((F(ii,:,kk)),F4_(ii,:,kk))
     end do
  end do
  do jj=1,ny
     do ii=1,nx
        call fft_real_1d((F(ii,jj,:)),F4_(ii,jj,:))
     end do
  end do

  call fftshift3d(cmplx(F1_,kind=8))
  open(unit=1,file='gauss_.dat',form='formatted')
  do kk=1,nz
     do jj=1,ny
        write(1,*) y(1,jj,1), z(1,1,kk), real(F1_(1,jj,kk))
     end do
  end do
  close(1)

end Program Gauss
