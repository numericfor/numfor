!> test_wfitpack
program test_wfitpack
  USE utils, only: dp, Zero, M_PI
  USE grids, only: linspace
  USE fitpack
  implicit none
  real(dp), dimension(:), allocatable :: x
  real(dp), dimension(:), allocatable :: y
  integer, parameter :: N = 12
  integer :: i
  integer :: ier

  type(UnivSpline) :: tck

  allocate (x(N), y(N))
  x = linspace(Zero, M_PI, N)
  y = sin(x)
  do i = 1, N
    print *, x(i), y(i)
  end do

  print *, ''
  print *, ''
  call splrep(x, y, tck=tck)
  ier = -1000
  call splrep(x, y, tck=tck, ier=ier)
  print *, tck%c
  print *, ier
end program test_wfitpack
