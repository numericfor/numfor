!> test
program test
  USE basic, only: dp
  USE sort, only: searchsorted
  implicit none
  real(dp), dimension(6) :: a = [1._dp, 3._dp, 5._dp, 9._dp, 32._dp, 124._dp]
  real(dp) :: elem
  integer :: n
  elem = 6.34_dp
  n = searchsorted(a, elem)

  print "(6(g0.4,1x))", a
  print *, 'n=', n
end program test
