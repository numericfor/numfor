program test
  use utils, only: dp, Zero, center
  use arrays, only: arange, linspace, geomspace, logspace, merge_sorted, mean, searchsorted, std
  implicit none

  integer, parameter :: N = 10
  ! real(dp), dimension(2 * N) :: a
  real(dp), dimension(:), allocatable :: a
  real(dp), dimension(:), allocatable :: b
  real(dp), dimension(:), allocatable :: c
  real(dp) :: x1, x2
  integer :: j

  x1 = -3._dp
  x2 = 3._dp

  print '(A)', center(' Test linspace and arange ', 65, '-')

  allocate (a(2 * N))
  a(:N) = linspace(0._dp, 1._dp, N, endpoint=.False.)
  a(N + 1:) = linspace(1._dp, 1.9_dp, N)

  print *, 'size(a)', size(a)
  print '(5(g8.3,1x))', a
  print '(A)', repeat('-', 65)
  print '(4(I3,1x))', arange(0, 122, 8)

  print *, 'mean, std', mean(a), std(a)

  print '(A)', center(' Test logspace and geomspace ', 65, '-')
  print '(4(f8.2,1x))', logspace(3._dp, Zero, 4)
  print '(4(f8.2,1x))', geomspace(-1000._dp, -1._dp, num=4)

  print '(A)', repeat('-', 65)

  b = geomspace(1.e-5_dp, 10._dp, 200)

  print *, 'mean, std', mean(b), std(b)

  print '(A)', repeat('-', 65)
  print '(A)', center(' Test searchsorted ', 65, '-')

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

  deallocate (b)
  b = geomspace(1.e-1_dp, 2._dp, 2 * N)
  print *, ""
  print '(A)', repeat('-', 30)//' a '//repeat('-', 30)
  print "(10(g0.3,1x))", a
  print '(A)', repeat('-', 30)//' b '//repeat('-', 30)
  print "(10(g0.3,1x))", b
  c = merge_sorted(a, b)
  print '(A)', repeat('-', 30)//' c '//repeat('-', 30)
  print "(10(g0.3,1x))", c
  print "(A,3(I2,1x))", 'Dimensiones: ', size(a), size(b), size(c)

  deallocate (b)
  b = geomspace(1.e-1_dp, 1.9_dp, 2 * N)
  print *, ""
  print '(A)', repeat('-', 30)//' a '//repeat('-', 30)
  print "(10(g0.3,1x))", a
  print '(A)', repeat('-', 30)//' b '//repeat('-', 30)
  print "(10(g0.3,1x))", b
  c = merge_sorted(a, b)
  print '(A)', repeat('-', 30)//' c '//repeat('-', 30)
  print "(10(g0.3,1x))", c
  print "(A,3(I2,1x))", 'Dimensiones: ', size(a), size(b), size(c)

  print '(A)', repeat('-', 28)//' tol<0 '//repeat('-', 28)
  c = merge_sorted(a, b, -1.e-300_dp)
  print '(A)', repeat('-', 30)//' c '//repeat('-', 30)
  print "(10(g0.3,1x))", c
  print "(A,3(I2,1x))", 'Dimensiones: ', size(a), size(b), size(c)

  deallocate (a, b, c)
  a = [1.9999999999_dp, 2._dp, 2._dp, 3._dp]
  b = [-1._dp, 2._dp, 3._dp, 4._dp]
  print *, ""
  print '(A)', repeat('-', 30)//' a '//repeat('-', 30)
  print "(4(g0.12,1x))", a
  print '(A)', repeat('-', 30)//' b '//repeat('-', 30)
  print "(4(g0.12,1x))", b
  c = merge_sorted(a, b)
  j = 30 - int(len("merge_sorted(a,b)") / 2._dp)
  print '(A)', repeat('-', 60)
  print '(A)', repeat('-', j)//' merge_sorted(a,b) '//repeat('-', j)
  print "(4(g0.12,1x))", c
  c = merge_sorted(a, b, 1.e-6_dp)
  j = 30 - int(len("merge_sorted(a,b,1.e-6)") / 2._dp)
  print '(A)', repeat('-', j)//' merge_sorted(a,b,1.e-6) '//repeat('-', j)
  print "(4(g0.12,1x))", c
  c = merge_sorted(a, b, -1.e-6_dp)
  j = 30 - int(len("merge_sorted(a,b,-1.e-6)") / 2._dp)
  print '(A)', repeat('-', j)//' merge_sorted(a,b,-1.e-6) '//repeat('-', j)
  print "(4(g0.12,1x))", c

end program test

