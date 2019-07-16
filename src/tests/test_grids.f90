program test
  use basic
  use strings, only: center
  use grids
  implicit none

  integer, parameter :: N = 10
  real(dp), dimension(2 * N) :: a
  real(dp), dimension(:), allocatable :: b
  real(dp) :: x1, x2
  integer :: j

  x1 = -3._dp
  x2 = 3._dp

  print '(A)', center(' Test linspace and arange ', 65, '-')

  a(:N) = logspace(x1, Zero, N, endpoint=.False.)
  a(N + 1:) = linspace(1._dp, x2, N)

  print *, 'size(a)', size(a)
  print '(5(g8.2,1x))', a
  print '(A)', repeat('-', 65)
  print '(4(I3,1x))', arange(0, 122, 8)

  print '(A)', center(' Test logspace and geomspace ', 65, '-')
  print '(4(f8.2,1x))', logspace(3._dp, Zero, 4)
  print '(4(f8.2,1x))', geomspace(-1000._dp, -1._dp, num=4)

  print '(A)', repeat('-', 65)
  print '(A)', center(' Test searchsorted ', 65, '-')

  b = geomspace(1.e-5_dp, 10._dp, 200)
  x1 = 5.12
  j = searchsorted(b, x1)
  print '(A,f8.5,A,I3,A,f8.5,A,f8.5)', 'Found ', x1, ' in position ', j,&
    & ' between', b(j), ' and', b(j + 1)

  x1 = 10._dp
  j = searchsorted(b, x1)
  print '(A,f8.5,A,I3,A,f8.5)', 'Found ', x1, ' in position ', j, ' after ', b(j)

  x1 = 0.0013_dp
  j = searchsorted(b, x1)
  print '(A,f8.5,A,I3,A,f8.5,A,f8.5)', 'Found ', x1, ' in position ', j, &
    &' between', b(j), ' and', b(j + 1)

end program test

