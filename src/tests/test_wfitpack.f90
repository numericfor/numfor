!> test_wfitpack
program test_wfitpack
  USE utils, only: dp, Zero, M_PI
  USE grids, only: linspace
  USE fitpack
  implicit none
  real(dp), dimension(:), allocatable :: x
  real(dp), dimension(:), allocatable :: y
  integer, parameter :: N = 6
  integer :: i
  integer :: u

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

  ! call splrep(x, y, tck=tck, ier=ier)
  open (newunit=u, file='splrep1.dat')
  do i = 1, size(tck%c)
    write (u, "(2(g0.14, 1x))") tck%c(i), tck%t(i)
  end do

  close (u)
end program test_wfitpack
