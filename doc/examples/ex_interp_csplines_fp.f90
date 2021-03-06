program cubic_splines_fp
  USE numfor, only: dp, linspace, save_array
  USE numfor, only: CubicSpline, csplrep, csplev, csplevder, csplint
  implicit none
  integer, parameter :: Ndim = 20, Ndimnew = 4 * Ndim
  type(CubicSpline) :: csp
  real(dp), dimension(Ndim) :: x, y
  real(dp), dimension(Ndimnew) :: xnew, ynew, yder
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
  call csplrep(x, y, -y(1), -y(Ndim), csp)
  ! Apply the interpolation to the new grid
  ynew = csplev(xnew, csp)
  !< [Using eval]
  ! save to file
  fname = "data/ex_interp_cspline3.dat"
  header = "x                sin(x)             evaluate"
  call save_array([xnew, sin(xnew), ynew], 3, fname=fname, fmt="f17.14", header=header)
  !< [Using derivative]
  ! Given a CubicSpline object csp evaluate the derivative
  yder = csplevder(xnew, csp)
  !< [Using derivative]
  ! save to file
  fname = "data/ex_interp_cspline4.dat"
  header = "x                cos(x)             derivative"
  call save_array([xnew, cos(xnew), yder], 3, fname=fname, fmt="f17.14", header=header)

  print "(A,f8.6)", "Value of the function at 1.5     = ", csplev(1.5_dp, csp)
  print "(A,f8.6)", "Value of first derivative at 1.5 = ", csplevder(1.5_dp, csp)
  print "(A,f8.6)", "Value of integral from 0 to 1.5  = ", csplint(0._dp, 1.5_dp, csp)
end program cubic_splines_fp
