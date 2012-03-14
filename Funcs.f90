module Funcs
  implicit none

  interface operator (.x.)
     module procedure cross
  end interface operator (.x.)

contains
  !Forward FFT algorithms
  !=======================================================!
  subroutine fft(data_in, data_out)
    !1D fft subroutine
    !This is really a wrapper for the intel library so that I don't have to deal with changing precision in the main program
    implicit none
    include "fftw3.f"
    complex(kind=8),intent(in) :: data_in(:)
    complex(kind=8) :: data_out(:)
    integer*8 :: plan
    integer :: ndata

    ndata=size(data_in)

    !shift the data
    !call fft_shift(in) 

    !Make Plan
    call dfftw_plan_dft_1d(plan,size(data_in),data_in,data_out,FFTW_FORWARD,FFTW_ESTIMATE)

    !Do the inverse transform
    call dfftw_execute_dft(plan,data_in,data_out)

    !Destroy the plan
    call dfftw_destroy_plan(plan)

    !shift again
    !call ifft_shift(out)

    !Normalize
    data_out = data_out/sqrt(ndata*1._8)
  end subroutine fft

  subroutine fft3d(data_in, data_out,shift)
    !This routine takes a (complex(kind=8)) 3d array as input and calls the 1D fft algorithm over each dimension of the array
    implicit none
    complex(kind=8), intent(in) :: data_in(:,:,:)
    complex(kind=8), intent(out) :: data_out(:,:,:)
    integer :: nx,ny,nz !size of the dimensions
    integer :: ix,iy,iz !counter variables for each dimension
    logical :: shift
    nx = size(data_in,1)
    ny = size(data_in,2)
    nz = size(data_in,3)

    if(shift)then
       !Shift the input array "along the half spaces"
       call ifftshift3d(data_in)
    end if
    !Call the fft's for all the dimensions
    do iz=1,nz
       do iy=1,ny
          call fft(data_in(:,iy,iz),data_out(:,iy,iz))
       end do
       do ix=1,nx
          call fft(data_in(ix,:,iz),data_out(ix,:,iz))
       end do
    end do
    do iy=1,ny
       do ix=1,nx
          call fft(data_in(ix,iy,:),data_out(ix,iy,:))
       end do
    end do

    if(shift)then
       !Shift back
       call fftshift3d(data_out)
    end if

  end subroutine fft3d

  subroutine fft_3d(data_in,data_out,shift)
    !fft routine which uses the 3d intel fft algorithm
    implicit none
    include "fftw3.f"
    !3d fftw
    complex(kind=8), intent(in) :: data_in(:,:,:)
    complex(kind=8) :: data_out(:,:,:)
    integer*8 :: plan
    integer :: nx,ny,nz
    logical :: shift

    nx = size(data_in,1)
    ny = size(data_in,2)
    nz = size(data_in,3)

    if(shift)then
       !Shift
       call ifftshift3d(data_in)
    end if
    
    call dfftw_plan_dft_3d(plan,nx,ny,nz,data_in,data_out,FFTW_FORWARD,FFTW_ESTIMATE)

    call dfftw_execute_dft(plan,data_in,data_out)

    call dfftw_destroy_plan(plan)
    
    if(shift)then
       !Shift back
       call fftshift3d(data_out)
    end if
    
    data_out = data_out/sqrt(nx*ny*nz*1._8)
  end subroutine fft_3d


  !==============================================================!
  !
  !Backwards Transforms
  !
  !==============================================================!
  subroutine ifft(data_in, data_out)
    !1D ifft subroutine
    !This is really a wrapper for the intel library so that I don't have to deal with changing precision in the main program
    implicit none
    include "fftw3.f"
    complex(kind=8),intent(in) :: data_in(:)
    complex(kind=8) :: data_out(:)
    integer*8 :: plan
    integer :: ndata

    ndata=size(data_in)

    !shift the data
    !call fft_shift(data_in) 

    !Make Plan
    call dfftw_plan_dft_1d(plan,size(data_in),data_in,data_out,FFTW_BACKWARD,FFTW_ESTIMATE)

    !Do the inverse transform
    call dfftw_execute_dft(plan,data_in,data_out)

    !Destroy the plan
    call dfftw_destroy_plan(plan)

    !shift again
    !call ifft_shift(data_out)

    !Normalize
    data_out = data_out/sqrt(ndata*1._8)
  end subroutine ifft

  subroutine ifft3d(data_in, data_out,shift)
    !This routine takes a (complex) 3d array as input and calls the 1D fft algorithm over each dimension of the array
    implicit none
    complex(kind=8), intent(in) :: data_in(:,:,:)
    complex(kind=8), intent(out) :: data_out(:,:,:)
    integer :: nx,ny,nz !size of the dimensions
    integer :: ix,iy,iz !counter variables for each dimension
    logical :: shift
    
    nx = size(data_in,1)
    ny = size(data_in,2)
    nz = size(data_in,3)

    if(shift)then 
       call ifftshift3d(data_in)
    end if

    !Call the fft's for all the dimensions
    do iz=1,nz
       do iy=1,ny
          call ifft(data_in(:,iy,iz),data_out(:,iy,iz))
       end do
       do ix=1,nx
          call ifft(data_in(ix,:,iz),data_out(ix,:,iz))
       end do
    end do
    do iy=1,ny
       do ix=1,nx
          call ifft(data_in(ix,iy,:),data_out(ix,iy,:))
       end do
    end do
    if(shift)then
       !Shift back
       call fftshift3d(data_out)
    end if
  end subroutine ifft3d

  subroutine ifft_3d(data_in,data_out,shift)
    implicit none
    include "fftw3.f"
    !3d fftw
    complex(kind=8), intent(in) :: data_in(:,:,:)
    complex(kind=8) :: data_out(:,:,:)
    integer*8 :: plan
    integer :: nx,ny,nz
    logical :: shift

    nx = size(data_in,1)
    ny = size(data_in,2)
    nz = size(data_in,3)

    if(shift)then
       !Shift the data along the "half spaces" (no idea what that really means)
       call fftshift3d(data_in)
    end if
    
    call dfftw_plan_dft_3d(plan,nx,ny,nz,data_in,data_out,FFTW_BACKWARD,FFTW_ESTIMATE)

    call dfftw_execute_dft(plan,data_in,data_out)

    call dfftw_destroy_plan(plan)

    if(shift)then
       !Shift back
       call ifftshift3d(data_out)
    end if
    
    data_out = data_out/sqrt(nx*ny*nz*1._8)
  end subroutine ifft_3d
  !===============================================================!
  !
  !Shift Functions
  !
  !===============================================================!
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
    middle = floor(fin/2._8)

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

!=============================================!
!
!PGPLOT plotting routine
!
!=============================================!
  subroutine plot(x,y,xaxis,yaxis,title)
    implicit none
    real(kind=4), intent(in) :: x(:), y(:)
    character(LEN=*) :: xaxis, yaxis, title
    call PGENV(minval(x),maxval(x),minval(y),maxval(y),0,0)
    call PGLAB(xaxis,yaxis,title)
    call PGLINE(size(x),x,y)
    return
  end subroutine plot

!==================================================!
!
!Vector products
!
!==================================================!
  function cross(A,B) result(C)
    implicit none
    complex(kind=8),intent(in), dimension(:,:,:,:) :: A,B
    complex(kind=8), allocatable, dimension(:,:,:,:) :: C
    integer :: ix,iy,iz
    ix = size(A,2)
    iy = size(A,3)
    iz = size(A,4)
    allocate(C(3,ix,iy,iz))
    !cx = ay*bz - az*by
    !cy = az*bx - ax*bz
    !cz = ax*by - ay*bx
    C(1,:,:,:) = A(2,:,:,:)*B(3,:,:,:) - A(3,:,:,:)*B(2,:,:,:)
    C(2,:,:,:) = A(3,:,:,:)*B(1,:,:,:) - A(1,:,:,:)*B(3,:,:,:)
    C(3,:,:,:) = A(1,:,:,:)*B(2,:,:,:) - A(2,:,:,:)*B(1,:,:,:)
  end function cross

  function dot(A,B) result(C)
    implicit none
    complex(kind=8),intent(in), dimension(:,:,:,:) :: A,B
    complex(kind=8), allocatable, dimension(:,:,:) :: C
    integer :: ix,iy,iz
    ix = size(A,2)
    iy = size(A,3)
    iz = size(A,4)
    allocate(C(ix,iy,iz))
    C(:,:,:) = A(1,:,:,:)*B(1,:,:,:) + A(2,:,:,:)*B(2,:,:,:) + A(3,:,:,:)*B(3,:,:,:)
  end function dot
  
  function grad(A,b) result(C)
    implicit none
    complex(kind=8),intent(in), dimension(:,:,:,:) :: A
    complex(kind=8),intent(in), dimension(:,:,:) :: b
    complex(kind=8), allocatable, dimension(:,:,:,:) :: C
    integer :: ix,iy,iz,count
    ix = size(A,2)
    iy = size(A,3)
    iz = size(A,4)
    allocate(C(3,ix,iy,iz))
    forall(count=1:3) C(count,:,:,:) = A(count,:,:,:)*b
  end function grad
  
  function average3(input) result(avg)
    !This function will take a 2 dimensional average
    !It will take in a 3d  array, take the average in each of the x-y planes and out a 3D array whose x-y planes are replaced by that average
    implicit none
    complex(kind=8),intent(in) :: input(:,:,:)
    complex(kind=8) :: avg(size(input,1),size(input,2),size(input,3))
    integer :: nx,ny,nz
    integer :: iz

    nx = size(input,1)
    ny = size(input,2)
    nz = size(input,3)
    
    do iz=1,nz
       avg(:,:,iz) = sum(input(:,:,iz))/real(nx,kind=8)/real(ny,kind=8)
    end do
  end function average3

    function average4(input) result(avg)
    !This function will take a 2 dimensional average
    !It will take in a 4d  array, take the average in each of the x-y planes and out a 4D array whose x-y planes are replaced by that average
    implicit none
    complex(kind=8),intent(in) :: input(:,:,:,:)
    complex(kind=8) :: avg(size(input,1),size(input,2),size(input,3),size(input,4))
    integer :: nx,ny,nz
    integer :: iz,id

    nx = size(input,1)
    ny = size(input,2)
    nz = size(input,3)
    do id=1,3
       do iz=1,nz
          avg(id,:,:,iz) = sum(input(id,:,:,iz))/real(nx,kind=8)/real(ny,kind=8)
       end do
    end do
  end function average4
  
end module Funcs
