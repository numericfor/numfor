program cubic_splines_oo
  USE numfor, only: dp, linspace, save_array, str, M_PI
  USE numfor, only: CubicSpline
  implicit none
  integer, parameter :: Ndim = 20, Ndimnew = 4 * Ndim
  type(CubicSpline) :: csp
  real(dp), dimension(Ndim) :: x, y
  real(dp), dimension(Ndimnew) :: xnew, ynew, yder
  real(dp), dimension(:), allocatable :: zeros
  character(len=:), allocatable :: header
  character(len=:), allocatable :: fname

  ! Create the data
  x = linspace(1.e-6_dp, 20._dp, Ndim)
  y = sin(x)
  ! New grid
  xnew = linspace(1.e-6_dp, 19.99_dp, Ndimnew)
  !< [Using eval]
  ! Create the CubicSpline
  ! Notice that we know exactly the second derivative
  csp = CubicSpline(x, y, -y(1), -y(Ndim))
  ! Apply the interpolation to the new grid
  ynew = csp%evaluate(xnew)
  !< [Using eval]
  ! save to file fname
  fname = "data/ex_interp_cspline1.dat"
  header = "x                sin(x)             evaluate"
  call save_array([xnew, sin(xnew), ynew], 3, fname=fname, fmt="f17.14", header=header)
  !< [Using derivative]
  ! Given a CubicSpline object csp, evaluate the derivative
  yder = csp%derivative(xnew)
  !< [Using derivative]
  ! save to file (stdout really)
  fname = "data/ex_interp_cspline2.dat"
  header = "x                cos(x)             derivative"
  call save_array([xnew, cos(xnew), yder], 3, fname=fname, fmt="f17.14", header=header)

  print "(A,f8.6)", "Value of the function at 1.5     = ", csp%evaluate(1.5_dp)
  print "(A,f8.6)", "Value of first derivative at 1.5 = ", csp%derivative(1.5_dp)
  print "(A,f8.6)", "Value of integral from 0 to 1.5  = ", csp%integrate(0._dp, 1.5_dp)

  zeros = csp%roots()
  header = "There are "//str(size(zeros))//" roots. In units of π:"
  call save_array(zeros / M_PI, 1, fmt='f13.10', header=header)

end program cubic_splines_oo
