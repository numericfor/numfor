program ex_linspace
  use numfor, only: dp, linspace, save_array
  implicit none

  integer, parameter :: N = 10
  real(dp), dimension(2 * N) :: a
  real(dp) :: step
  !< [linspace]
  print "(5(f4.2,1x))", linspace(2._dp, 3.0_dp, num=5)
  ! Gives: 2.00 2.25 2.50 2.75 3.00
  print "(5(f4.2,1x))", linspace(2._dp, 3.0_dp, num=5, endpoint=.False.)
  ! Gives: 2.00 2.20 2.40 2.60 2.80
  print "(5(f4.2,1x))", linspace(2._dp, 3.0_dp, num=5, retstep=step)
  ! Gives: 2.00 2.25 2.50 2.75 3.00
  !< [linspace]
  print "(A)", repeat('-', 20)
  a(:N) = linspace(0._dp, 0.5_dp, N, endpoint=.False.)
  a(N + 1:) = linspace(0.5_dp, 1.4_dp, N)
  call save_array(a, 5, fmt='f8.5')
end program ex_linspace

