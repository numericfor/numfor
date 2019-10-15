!> Example of use of splprep
program ex_splprep
  USE numfor, only: dp, Zero, M_PI, str, linspace, save_array
  USE numfor, only: UnivSpline, splprep, splev
  implicit none

  real(dp), dimension(:), allocatable :: r, phi, u
  real(dp), dimension(:, :), allocatable :: x, new_points
  type(UnivSpline) :: tck
  character(len=:), allocatable :: header
  character(len=:), allocatable :: fname

  integer :: Nd = 40          ! Number of points
  real(dp) :: s = 0._dp

  allocate (r(Nd), u(Nd), phi(Nd))
  allocate (x(2, Nd), new_points(2, Nd))

  ! Create data
  phi = linspace(Zero, 2 * M_PI, Nd)
  r = 0.5_dp + cos(phi)        ! polar coords
  x(1, :) = r * cos(phi)      ! convert to cartesian
  x(2, :) = r * sin(phi)      ! convert to cartesian
  !< [Using it]
  ! After setting data in x
  ! interpolate: (step one) Create interpolation and parameters u
  call splprep(x, u, tck, s=s)
  ! interpolate: (step two) Apply to points u
  call splev(u, tck, new_points)
  !< [Using it]
  ! save to file
  fname = "data/ex_interp_splprep1.dat"
  header = "u                    x                   y"
  call save_array([u, new_points(1, :), new_points(2, :)], 3, fname=fname, fmt="f17.14", header=header)
end program ex_splprep

