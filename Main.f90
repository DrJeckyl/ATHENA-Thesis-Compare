!Program for taking in simulation data that has already been broken down into a column format text file
!The goal is to calculate:
!j . b = \bar{J . B} - \bar{J} . \bar{B}
!v x b = \bar{v x B} - \bar{v} x \bar{B}
!a . b = \bar{a . b} - \bar{a} . \bar{b}
!A x ( v x B + grad(phi)) - \bar{A} x (\bar{v x B} + grad{\bar{phi}})
!grad^2 phi = div (v x b)
!grad^2 a = -curl(B)
Program Main
IMPLICIT NONE
real, dimension(32*128*128) :: bx_raw,by_raw,bz_raw,ux_raw,uy_raw,uz_raw
real :: x(32),y(128),z(128)
real, dimension(32,128,128) :: bx,by,bz,ux,uy,uz

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
bx = reshape(bx_raw,(/32,128,128/))
by = reshape(by_raw,(/32,128,128/))
bz = reshape(bz_raw,(/32,128,128/))
ux = reshape(ux_raw,(/32,128,128/))
uy = reshape(uy_raw,(/32,128,128/))
uz = reshape(uz_raw,(/32,128,128/))

!Ok now everything is put into a big matrix
!Want to Fourier Transform these arrays

!Define a k-vector

!Calculate in k-space
!j = curl(B)
!j = k x B


!Inverse Fourier Transform
End Program Main
