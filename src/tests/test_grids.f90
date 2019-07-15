program test
  use basic
  use grids
  implicit none

  integer, parameter :: N = 10
  real(dp), dimension(2 * N) :: a
  real(dp), dimension(:), allocatable :: b
  real(dp) :: x1, x2

  x1 = -3._dp
  x2 = 3._dp

  a(:N) = logspace(x1, Zero, N, endpoint=.False.)
  a(N:) = linspace(1._dp, x2, N)
  print *, 'size(a)', size(a)
  print '(5(g8.2,1x))', a
  print '(A)', repeat('-', 65)

  print *, linspace(Zero, 3._dp, 4)
  print *, logspace(3._dp, Zero, 4)

  b = geomspace(-1000._dp, -1._dp, num=4)
  ! print *, 'size(b)', size(b)
  print '(4(g0,1x))', b
end program test

