!> test_interpolate
program test_polinomial
  USE utils, only: dp, Zero, M_PI
  USE arrays, only: arange
  USE polynomial
  implicit none

  real(dp), dimension(6) :: p1
  ! real(dp), dimension(:), allocatable :: p2

  p1 = 1._dp
  print "(6(f7.4,1x))", p1

  ! print *, polyval(p1, -1._dp)
  ! p2 = polyder(p1)
  ! print *, polyder(p1, 1)
  print "(5(f7.3,1x))", polyder(p1, 1)
  print "(4(f7.3,1x))", polyder(p1, 2)
  print "(3(f7.3,1x))", polyder(p1, 3)
  print "(2(f7.3,1x))", polyder(p1, 4)
  print "(1(f7.3,1x))", polyder(p1, 5)

end program test_polinomial
