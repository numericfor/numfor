!> test the number generator,
!>
!> Results for the coded values should be
!>@code
!> seed:            1447774167
!> 3.30885E-02  1.98364E-01  2.90381E-01
!> 4.75757E-01  5.08138E-01  1.34913E-01
!>
!> 5.52038E-01  1.98364E-01  2.90381E-01
!> 6.59989E-02  5.08138E-01  1.34913E-01
!>
!> 5.52038E-01  1.98364E-01  6.76710E-01
!> 6.59989E-02  5.08138E-01  7.77963E-01
!>
!> 8.92891E-01  6.01456E-02  9.18351E-01
!> 6.59989E-02  5.08138E-01  7.77963E-01
!> test finished OK:   0.153604707483576
!>@endcode
!>  Run Time is approximately: 0m12.509s
program test_random
  USE basic, only: dp, i8, str
  USE random, only: random_seed, random_sample
  integer, parameter :: TEST_SIZE = 400000000
  integer :: semilla
  ! integer(i8) :: i1
  real(dp), dimension(2, 3) :: r1
  real(dp) :: r

  semilla = 1234
  call random_seed(semilla)
  ! Single value
  call random_sample(r)
  print "(A)", str(r)

  ! Vector of random numbers
  call random_sample(r1(:, 1))
  print "(A)", str(r1(:, 1))

  ! Array of random numbers
  call random_sample(r1)
  print "(A)", str(r1)

  if (abs(r - 0.78682095486780212_dp) < epsilon(1._dp)) then
    print "(A,f20.17)", 'test finished OK: r =', r
  else
    print "(A,f20.17)", 'test failed, got r ='//str(r)//' instead of ', 0.78682095486780212_dp
  end if
end program test_random

