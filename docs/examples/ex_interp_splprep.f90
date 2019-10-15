!> Example of use of splprep
program ex_splprep
  USE numfor, only: dp, Zero, M_PI, str, linspace, save_array
  USE numfor, only: UnivSpline, splprep, splev
  implicit none

  real(dp), dimension(:), allocatable :: r, phi, u
  real(dp), dimension(:, :), allocatable :: x, new_points
  type(UnivSpline) :: tck
  character(len=:), allocatable :: header
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
  header = "u                    x                   y"
  call save_array([u, new_points(1, :), new_points(2, :)], 3, fmt="f17.14", header=header)
  ! Prints:
  !
  ! # u                    x                   y
  ! 0.00000000000000  1.50000000000000  0.00000000000000
  ! 0.03616230734577  1.46779335229253  0.23853963732962
  ! 0.07211603380768  1.37398960267532  0.45870512901973
  ! 0.10765440708122  1.22676038619218  0.64385351896871
  ! 0.14257428865429  1.03883011365998  0.78063018793468
  ! 0.17667803554905  0.82622920670009  0.86019572275769
  ! 0.20977542560675  0.60672992984431  0.87900005428954
  ! ....
  ! ....
  ! 0.79022457439325  0.60672992984431 -0.87900005428954
  ! 0.82332196445095  0.82622920670009 -0.86019572275769
  ! 0.85742571134571  1.03883011365998 -0.78063018793468
  ! 0.89234559291878  1.22676038619218 -0.64385351896871
  ! 0.92788396619232  1.37398960267532 -0.45870512901973
  ! 0.96383769265423  1.46779335229253 -0.23853963732962
  ! 1.00000000000000  1.50000000000000 -0.00000000000000
end program ex_splprep

