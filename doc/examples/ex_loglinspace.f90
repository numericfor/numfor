program ex_loglinspace
  use numfor, only: dp, Zero, str, loglinspace, save_array
  implicit none

  character(len=:), allocatable :: header
  character(len=:), allocatable :: fname

  integer, parameter :: N = 100
  real(dp), dimension(N, 4) :: a
  real(dp), dimension(N) :: v
  real(dp), dimension(4) :: steps
  real(dp) :: step, ratio
  real(dp) :: xmin, xmax
  integer :: i
  header = ""

  !< [loglinspace]
  xmin = 1.e-4_dp; xmax = 10._dp
  ratio = 1.15_dp; steps = [0.125_dp, 0.25_dp, 0.5_dp, 1._dp]
  v = loglinspace(xmin, xmax, 50, step=0.4_dp, ratio=1.15_dp)
  !< [loglinspace]

  do i = 1, 4
    step = steps(i)
    a(:, i) = loglinspace(xmin, xmax, N, step=step, ratio=ratio)
    header = header//"     step = "//str(step)
  end do

  fname = "data/ex_loglinspace.dat"
  call save_array([a(:, 1), a(:, 2), a(:, 3), a(:, 4)], 4, fname=fname, fmt='f16.13', header=header)
end program ex_loglinspace

