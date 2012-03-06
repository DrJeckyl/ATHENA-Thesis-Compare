module Funcs
  implicit none

contains
  subroutine fftshift3d(stuffin)
    implicit none
    complex(kind=8) :: stuffin(:,:,:)
    integer :: numx,numy,numz,ix,iy,iz
    numx = size(stuffin,1)
    numy = size(stuffin,2)
    numz = size(stuffin,3)
    do iz=1,numz
       do iy=1,numy
          call fft_shift(stuffin(:,iy,iz))
       end do
    end do
    do iz=1,numz
       do ix=1,numx
          call fft_shift(stuffin(ix,:,iz))
       end do
    end do
    do iy=1,numy
       do ix=1,numx
          call fft_shift(stuffin(ix,iy,:))
       end do
    end do
  end subroutine fftshift3d

  subroutine ifftshift3d(stuffin)
    implicit none
    complex(kind=8) :: stuffin(:,:,:)
    integer :: numx,numy,numz,ix,iy,iz

    numx = size(stuffin,1)
    numy = size(stuffin,2)
    numz = size(stuffin,3)
    do iz=1,numz
       do iy=1,numy
          call ifft_shift(stuffin(:,iy,iz))
       end do
    end do
    do iz=1,numz
       do ix=1,numx
          call ifft_shift(stuffin(ix,:,iz))
       end do
    end do
    do iy=1,numy
       do ix=1,numx
          call ifft_shift(stuffin(ix,iy,:))
       end do
    end do
  end subroutine ifftshift3d


  function t_diff(time_start, time_end)
    implicit none
    integer, dimension(8), intent(in) :: time_start, time_end
    integer, dimension(8) :: t_diff
    integer, dimension(8) :: dummy_time
    !time_in is a array of 8 numbers
    !order: year, month, day, UTC offset, hours, minutes, seconds, miliseconds
    !Calculate the number of minutes difference between time_start and time_end and convert to hours/minutes/seconds

    !year,month can be considered monotone, so just subtract
    dummy_time(1) = time_end(1) - time_start(1)
    dummy_time(2) = time_end(2) - time_start(2)

    !If the year and month haven't changed, subtract the days
    if(dummy_time(1) == 0 .and. dummy_time(2) == 0)then
       dummy_time(3) = time_end(3) - time_start(3)
    else
       !give up
       write(*,*) 'I give up!'
       stop
    end if

    !If day hasn't changed, subtract the hours
    if(dummy_time(3) == 0)then
       dummy_time(5) = time_end(5) - time_start(5)
    else
       !Day changed, add 24 to time_end and subtract 1 from days
       dummy_time(5) = time_end(5) + 24 - time_start(5)
       dummy_time(4) = dummy_time(4) - 1
    end if

    !If hour hasn't changed, subtract the minutes
    if(dummy_time(5) == 0)then
       dummy_time(6) = time_end(6) - time_start(6)
    else
       !hour changed, add 60 and subtract
       dummy_time(6) = time_end(6) - time_start(6) + 60
       dummy_time(5) = dummy_time(5) -1
    end if

    !If minutes haven't changed, subtract the seconds
    if(dummy_time(6) == 0)then
       dummy_time(7) = time_end(7) - time_start(7)
    else
       !Minutes changed, subtract after adding 60*minute_change
       dummy_time(7) = time_end(7) - time_start(7) + 60
       dummy_time(6) = dummy_time(6) - 1
    end if

    !If seconds haven't changed, subtract the miliseconds
    if(dummy_time(7) == 0)then
       dummy_time(8) = time_end(8) - time_start(8)
    else
       !You get the idea
       dummy_time(8) = time_end(8) - time_start(8) + 1000
       dummy_time(7) = dummy_time(7) - 1
    end if

    !Check all the "bounds"
    if(dummy_time(8) > 1000)then
       dummy_time(8) = dummy_time(8) - 1000
       dummy_time(7) = dummy_time(7) + 1
    end if

    if(dummy_time(7) > 60)then
       dummy_time(7) = dummy_time(7) - 60
       dummy_time(6) = dummy_time(6) + 1
    end if

    if(dummy_time(6) > 60)then
       dummy_time(6) = dummy_time(6) - 60
       dummy_time(5) = dummy_time(5) + 1
    end if


    t_diff = dummy_time
  end function t_diff

  subroutine ifft_3d(data_in,data_out)
    implicit none
    include "fftw3.f"
    !3d fftw
    complex, intent(in) :: data_in(:,:,:)
    complex :: data_out(:,:,:)
    complex(kind=8),allocatable :: in(:,:,:),out(:,:,:)
    integer*8 :: plan
    integer :: nx,ny,nz

    nx = size(data_in,1)
    ny = size(data_in,2)
    nz = size(data_in,3)

    allocate(in(nx,ny,nz))
    allocate(out(nx,ny,nz))

    in = cmplx(0.,kind=8)
    out = cmplx(0.,kind=8)

    in = cmplx(data_in,kind=8)

    call dfftw_plan_dft_3d(plan,nx,ny,nz,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)

    call dfftw_execute_dft(plan,in,out)

    call dfftw_destroy_plan(plan)

    out = cmplx(out/sqrt(nx*ny*nz*1.),kind=8)
    data_out = out
  end subroutine ifft_3d

  subroutine ifft(data_in, data_out)
    implicit none
    include "fftw3.f"
    !This is a rewrite of the original ifft subroutine to incorporate Mischa's idea and include the shift function in the routine.
    integer :: ndata,n_start,n_end,middle
    complex,intent(in) :: data_in(:) 
    complex :: data_out(:)
    complex(kind=8) :: in(size(data_in)), out(size(data_out))
    integer*8 :: plan

    ndata = size(data_in)
    middle = size(in)/2
    n_start = middle -  floor(ndata/2.)
    n_end = middle +  floor(ndata/2.)

    !Make a dummy array that is pad x the size of the input array to pad zeros
    in = cmplx(0.,kind=8)
    out = cmplx(0.,kind=8)

    !First change the precision to standard double and put it in the middle
    in = cmplx(data_in,kind=8)

    !shift the data
    !call fft_shift(in)

    !Make the plan
    call dfftw_plan_dft_1d(plan,size(in),in,out,FFTW_BACKWARD,FFTW_ESTIMATE)

    !Do the inverse transform
    call dfftw_execute_dft(plan,in,out)

    !Destroy the plan
    call dfftw_destroy_plan(plan)

    !shift again
    !call ifft_shift(out)

    !Normalize
    out = out/sqrt(ndata*1.)

    !data_out = out(ndata*n_start+1:ndata*n_end)*1.
    data_out = out
  end subroutine ifft


!!$  subroutine fft(data_in,data_out)
!!$    
!!$    use TypesAndDefs
!!$    implicit none
!!$    include "fftw3.f"
!!$    integer :: ndata
!!$    complex, dimension(:) :: data_in
!!$    complex, dimension(:) :: data_out
!!$    complex(kind=8), dimension(size(data_in)) :: in, out
!!$    ndata=size(data_in)
!!$    in = cmplx(data_in,kind=8)
!!$    out = cmplx(data_out,kind=8)
!!$    !in is the data you want to transform
!!$    !Do a forward transform
!!$    call dfftw_execute_dft(plan,in,out)
!!$    !Normalize
!!$    out=out/sqrt(ndata*1.)
!!$    data_out = out
!!$  end subroutine fft

  subroutine fft(data_in, data_out)
    implicit none
    include "fftw3.f"
    !This is a rewrite of the original ifft subroutine to incorporate Mischa's idea and include the shift function in the routine.
    complex,intent(in) :: data_in(:)
    complex :: data_out(:)
    complex(kind=8) :: in(size(data_in)), out(size(data_out))
    integer*8 :: plan
    integer :: ndata

    ndata=size(data_in)

    !Make a dummy array that is 11x the size of the input array to pad zeros
    in = cmplx(0.,kind=8)
    out = cmplx(0.,kind=8)

    !First change the precision to standard double and put it in the middle
    !in(ndata*n_start+1:ndata*n_end) = cmplx(data_in,kind=8)
    !out(ndata*n_start+1:ndata*n_end) = cmplx(data_out,kind=8)
    in = cmplx(data_in,kind=8)

    !shift the data
    !call fft_shift(in) 

    !Make Plan
    call dfftw_plan_dft_1d(plan,size(in),in,out,FFTW_FORWARD,FFTW_ESTIMATE)

    !Do the inverse transform
    call dfftw_execute_dft(plan,in,out)

    !Destroy the plan
    call dfftw_destroy_plan(plan)

    !shift again
    !call ifft_shift(out)

    !Normalize
    out = cmplx(out/sqrt(ndata*1.),kind=8)

    !data_out = out(ndata*n_start+1:ndata*n_end)*1.
    data_out = out!(n_start:n_end)
  end subroutine fft

  subroutine fft_3d(data_in,data_out)
    implicit none
    include "fftw3.f"
    !3d fftw
    complex, intent(in) :: data_in(:,:,:)
    complex :: data_out(:,:,:)
    complex(kind=8),allocatable :: in(:,:,:),out(:,:,:)
    integer*8 :: plan
    integer :: nx,ny,nz

    nx = size(data_in,1)
    ny = size(data_in,2)
    nz = size(data_in,3)

    allocate(in(nx,ny,nz))
    allocate(out(nx,ny,nz))

    in = cmplx(data_in,kind=8)

    call dfftw_plan_dft_3d(plan,nx,ny,nz,in,out,FFTW_FORWARD,FFTW_ESTIMATE)

    call dfftw_execute_dft(plan,in,out)

    call dfftw_destroy_plan(plan)

    out = out/sqrt(nx*ny*nz*1.)
    data_out = cmplx(out)
  end subroutine fft_3d

  subroutine fft_shift(stuff_in)

    implicit none
    complex(kind=8) :: stuff_in(:)
    complex(kind=8) :: dummy(size(stuff_in))
    integer :: fin, middle

    fin = size(stuff_in)
    middle = ceiling(fin/2.)

    if(mod(fin,2) == 0)then
       dummy(1:middle) = stuff_in(middle+1:fin)
       dummy(middle+1:fin) = stuff_in(1:middle)
    else
       dummy(1:middle) = stuff_in(middle:fin)
       dummy(middle+1:fin) = stuff_in(1:middle-1)
    end if
    stuff_in = dummy
  end subroutine fft_shift

  subroutine ifft_shift(stuff_in)

    implicit none
    complex(kind=8) :: stuff_in(:)
    complex(kind=8) :: dummy(size(stuff_in))
    integer :: fin
    integer :: middle

    fin = size(stuff_in)
    middle = floor(fin/2.)

!!$    dummy(1:middle) = stuff_in(middle+1:fin)
!!$    dummy(middle+1:fin) = stuff_in(1:middle)
!!$    stuff_in = dummy
!!$

    if(mod(fin,2) == 0)then
       !Its even, the fft_shifts are symmetric
       call fft_shift(stuff_in)
    else
       !Its odd, so we do the inverse of the shift we did before
       dummy(middle+1:fin) = stuff_in(1:middle+1)
       dummy(1:middle) = stuff_in(middle+2:fin)
       stuff_in = dummy
    end if
  end subroutine ifft_shift

  function tridiag(a,b,c,d) result(x)
    implicit none
    real, intent(in), dimension(:) :: a,b,c
    complex,intent(in), dimension(:) :: d
    complex, dimension(size(a)) :: x
    real,dimension(size(a)) :: at,bt,ct,dtR,dtI,xR,xI
    real :: id
    integer(kind=4) :: m,n
    complex :: i
    i = (0.,1.)

    at = a
    bt = b
    ct = c
    dtR = real(d)
    dtI = aimag(d)

    n=size(a)
    ct(1) = ct(1)/bt(1)
    dtR(1) = dtR(1)/bt(1)
    do m=2,n
       id = 1.0/( bt(m) - ct(m-1)*at(m) )
       ct(m) = ct(m)*id
       dtR(m) = ( dtR(m) - dtR(m-1)*at(m) )*id
    end do    !Back substitute
    xR(n) = dtR(n)
    do m=n-1,1,-1
       xR(m) = dtR(m) - ct(m)*xR(m+1)
    end do

    at = a
    bt = b
    ct = c

    ct(1) = ct(1)/bt(1)
    dtI(1) = dtI(1)/bt(1)
    do m=2,n
       id = 1.0/( bt(m) - ct(m-1)*at(m) )
       ct(m) = ct(m)*id
       dtI(m) = ( dtI(m) - dtI(m-1)*at(m) )*id
    end do
    !Back substitute
    xI(n) = dtI(n)
    do m=n-1,1,-1
       xI(m) = dtI(m) - ct(m)*xI(m+1)
    end do

    x = xR + i*xI
  end function tridiag

  function Romberg(F,a,b)
    use TypesAndDefs
    implicit none
    complex Romberg
    complex :: F(:)
    real, intent(in) :: a, b
    integer :: j,k1
    complex :: R1(N+1,N+1)
    real :: tol

    tol = 2.0**-30

    R1(1,1) = trap(F,a,b,2**0)
    do k1=2,N+1
       R1(k1,1) = traprefine(R1(k1-1,1),F,a,b,2**(k1-1))
       !R(k,1) = trap(F,a,b,2**(k-1))
       !Calculate the kth row
       do j = 2,k1
          R1(k1,j) = (R1(k1,j-1)*4.**(j-1)-R1(k1-1,j-1))/(4.**(j-1) - 1.)
       end do
       !check for early convergence
       !if(abs(R(k,k)-R(k,k-1))<tol)then
       !   Romberg = R(k,k)
       !write(*,*) R(k,k), R(k,k-1)
       !   return
       !end if
    end do

    Romberg = R1(N+1,N+1)
    return
  end function Romberg


  function trap(F,a,b,p)
    implicit none
    complex trap
    real :: a,b,h
    complex :: F(:)
    complex :: S
    integer :: p
    integer :: j
    real :: rend
    rend = size(F)
    S = (0.0,0.0)
    h = (b-a)/p
    do j=1,rend,int((rend-1)/p)
       if(j==1)then
          S = S + F(1)
       end if
       if(j==rend)then
          S = S + F(rend)
       end if
       S = S + 2.0*F(j)
    end do
    trap = h*S/2.0
    return
  end function trap
  function traprefine(S,F,a,b,p) result(T)
    implicit none
    real :: a,b,h
    complex :: F(:), S
    integer :: temp, j, p, rend
    complex :: dummy
    complex :: T
    dummy = S
    h=(b-a)/p
    dummy = dummy/2.0
    rend=size(F)
    do j=1,p,2
       temp = (size(F)-1)/p*j + 1
       if(temp > rend)then
          write(*,*) 'Out of bounds in traprefine'
          exit
       end if
       dummy = dummy + h*F(temp)
    end do
    T = dummy
  end function traprefine

  subroutine plot(x,y,xaxis,yaxis,title)
    implicit none
    real(kind=4), intent(in) :: x(:), y(:)
    character(LEN=*) :: xaxis, yaxis, title
    call PGENV(minval(x),maxval(x),minval(y),maxval(y),0,0)
    call PGLAB(xaxis,yaxis,title)
    call PGLINE(size(x),x,y)
    return
  end subroutine plot

  function mag(vec)
    implicit none
    complex, intent(in) :: vec(:,:)
    complex :: mag(size(vec,1))
    mag = sqrt(sum(vec*conjg(vec)))
  end function mag

  subroutine fft_real_1d(data_in,data_out)
    implicit none
    include "fftw3.f"
    real :: data_in(:)
    complex :: data_out(:)
    real(kind=8) :: in(size(data_in))
    complex(kind=8) :: out(size(data_out))
    integer*8 :: plan
    integer :: ndata

    ndata=size(data_in)

    !Make a dummy array that is 11x the size of the input array to pad zeros
    in = real(0.,kind=8)
    out = cmplx(0.,kind=8)

    in = real(data_in,kind=8)

    !shift the data
    call ifft_shift(cmplx(in,kind=8))

    !Make Plan
    call dfftw_plan_dft_r2c_1d(plan,size(in),in,out,FFTW_ESTIMATE)

    !Do the inverse transform
    call dfftw_execute_dft_r2c(plan,in,out)

    !Destroy the plan
    call dfftw_destroy_plan(plan)

    !shift the data
    call fft_shift(cmplx(out,kind=8))

    !Normalize
    out = cmplx(out/sqrt(ndata*1.),kind=8)

    data_out = out
  end subroutine fft_real_1d

  subroutine ifft_real_1d(data_in,data_out)
    implicit none
    include "fftw3.f"
    complex :: data_in(:)
    real :: data_out(:)
    complex(kind=8) :: in(size(data_in))
    real(kind=8) :: out(size(data_out))
    integer*8 :: plan
    integer :: ndata
    ndata=size(data_in)

    in = cmplx(0.,kind=8)
    out = real(0.,kind=8)

    in = cmplx(data_in,kind=8)

    !shift the data
    !call ifft_shift(cmplx(in,kind=8))

    !Make Plan
    call dfftw_plan_dft_c2r_1d(plan,size(in),in,out,FFTW_ESTIMATE)

    !Do the inverse transform
    call dfftw_execute_dft_c2r(plan,in,out)

    !Destroy the plan
    call dfftw_destroy_plan(plan)

    !shift the data
    !call ifft_shift(cmplx(out,kind=8))

    !Normalize
    out = real(out/sqrt(ndata*1.),kind=8)

    data_out = out
  end subroutine ifft_real_1d

  subroutine fft_real_3d(data_in,data_out)
    implicit none
    include "fftw3.f"
    real :: data_in(:,:,:)
    complex :: data_out(:,:,:)
    real(kind=8),allocatable :: in(:,:,:)
    complex(kind=8),allocatable :: out(:,:,:)
    integer*8 :: plan
    integer :: nx,ny,nz
    nx = size(data_in,1)
    ny = size(data_in,2)
    nz = size(data_in,3)

    allocate(in(nx,ny,nz))
    allocate(out(nx,ny,nz))

    in = real(0.,kind=8)
    out = cmplx(0.,kind=8)

    in = real(data_in,kind=8)

    !Make Plan
    call dfftw_plan_dft_r2c_3d(plan,nx,ny,nz,in,out,FFTW_ESTIMATE)

    !Do the inverse transform
    call dfftw_execute_dft_r2c(plan,in,out)

    !Destroy the plan
    call dfftw_destroy_plan(plan)

    !Normalize
    out = cmplx(out/sqrt(nx*ny*nz*1.),kind=8)

    data_out = out
  end subroutine fft_real_3d

  subroutine ifft_real_3d(data_in,data_out)
    implicit none
    include "fftw3.f"
    real :: data_out(:,:,:)
    complex :: data_in(:,:,:)
    real(kind=8),allocatable :: out(:,:,:)
    complex(kind=8),allocatable :: in(:,:,:)
    integer*8 :: plan
    integer :: nx,ny,nz
    nx = size(data_in,1)
    ny = size(data_in,2)
    nz = size(data_in,3)

    allocate(in(nx,ny,nz))
    allocate(out(nx,ny,nz))

    out = real(0.,kind=8)
    in = cmplx(0.,kind=8)

    in = cmplx(data_in,kind=8)

    !Make Plan
    call dfftw_plan_dft_c2r_3d(plan,nx,ny,nz,in,out,FFTW_ESTIMATE)

    !Do the inverse transform
    call dfftw_execute_dft_c2r(plan,in,out)

    !Destroy the plan
    call dfftw_destroy_plan(plan)

    !Normalize
    out = real(out/sqrt(nx*ny*nz*1.),kind=8)

    data_out = out
  end subroutine ifft_real_3d


end module Funcs
