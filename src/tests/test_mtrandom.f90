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
program test
  integer, parameter :: TEST_SIZE = 400000000
  type(rng_state) :: state
  integer(4) :: semilla
  integer(I64) :: i1
  real(dp), dimension(2, 3) :: r1
  real(dp) :: r
  semilla = 1447774167

  call random(state, r1)

  ! open(UNIT=9, file='test_mtrandom.dat')
  write (*, '(A)') sep_2
  write (*, *) '# random seed:', state%mt(0)

  call seed(semilla, state)
  write (*, *) '# user-defined seed:', state%mt(0)
  call random(state, r1)
  write (*, '(3(1PE12.5,1x))') transpose(r1)
  write (*, *) ''
  call random(state, r1(:, 1))
  write (*, '(3(1PE12.5,1x))') transpose(r1)
  write (*, *) ''
  call random(state, r1(:, 3))
  write (*, '(3(1PE12.5,1x))') transpose(r1)
  write (*, *) ''
  call random(state, r1(1, :))
  write (*, '(3(1PE12.5,1x))') transpose(r1)
  do i1 = 1, TEST_SIZE
    call random(state, r)
  end do
  if (abs(r - 0.153604707483576) < 1.e-7) then
    write (*, *) 'test finished OK: ', r
  else
    write (*, *) 'test failed: ', r
  end if
  write (*, '(A)') sep_2
end program test

