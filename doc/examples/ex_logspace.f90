program ex_linspace
  use numfor, only: dp, Zero, logspace, geomspace, save_array
  implicit none

  integer, parameter :: N = 10
  real(dp), dimension(2 * N) :: a
  !< [logspace]
  print '(4(f8.2,1x))', logspace(2._dp, 3.0_dp, num=4)
  ! gives:  100.00   215.44   464.16  1000.00
  print '(4(f8.2,1x))', logspace(2.0_dp, 3.0_dp, num=4, base=2.0_dp)
  ! gives:    4.00     5.04     6.35     8.00
  print '(4(f8.2,1x))', logspace(3._dp, Zero, 4)
  ! gives: 1000.00   100.00    10.00     1.00
  !< [logspace]
  print "(A)", repeat('-', 20)
  !< [geomspace]
  print '(4(f8.2,1x))', geomspace(1._dp, 1000._dp, num=4)
  ! gives:     1.00    10.00   100.00  1000.00
  print '(4(f8.2,1x))', geomspace(-1000._dp, -1._dp, num=4)
  ! gives: -1000.00  -100.00   -10.00    -1.00
  !< [geomspace]
  print "(A)", repeat('-', 20)
  a = geomspace(1.e-5_dp, 10._dp, 20)
  call save_array(a, 5, fmt='f9.6')
end program ex_linspace

