  real(dp), intent(IN), dimension(2) :: xrange !< Limits in first axis
  real(dp), intent(IN), dimension(2) :: yrange !< Limits in second axis
  real(dp), intent(IN) :: eps !< Desired precision (estimated)
  integer, optional, intent(IN) :: limit !< Maximum number of function evaluations for each axis

  integer :: N
  integer :: icount
  integer, parameter :: NMAX = 10

  N = 21; IF (Present(limit)) N = limit
  icount = 1
  int0 = -huge(1._dp)
  int = NUMFOR_INTEG(func, xrange, yrange, N, N)
  do icount = 1, NMAX
    int0 = int
    N = 2 * N - 10
    int = NUMFOR_INTEG(func, xrange, yrange, N, N)
    IF (abs(int - int0) <= eps * abs(int)) exit
  enddo

