program ex_splprep
  USE numfor, only: dp, Zero, M_PI, str, linspace, save_array
  USE numfor, only: UnivSpline, splrep, splev
  implicit none

  integer, parameter :: N = 6
  integer, parameter :: Nnew = 29
  real(dp), dimension(N) :: x
  real(dp), dimension(N) :: y
  real(dp), dimension(Nnew) :: xnew
  real(dp), dimension(Nnew) :: ynew
  character(len=:), allocatable :: header
  character(len=:), allocatable :: fname

  real(dp) :: s
  type(UnivSpline) :: tck

  ! Create data
  s = 0._dp
  x = linspace(Zero, M_PI, N)
  y = sin(x)
  xnew = linspace(Zero, M_PI, Nnew)
  !< [Using it]
  ! After setting data in x
  ! interpolate: (step one) Create interpolation and parameters u
  call splrep(x, y, tck=tck, s=s)
  ! interpolate: (step two) Apply to points u
  call splev(xnew, tck, ynew)
  !< [Using it]
  ! save to file
  fname = "data/ex_interp_splrep1.dat"
  header = "x                   y"
  call save_array([xnew, ynew], 2, fname=fname, fmt="f17.14", header=header)
end program ex_splprep
