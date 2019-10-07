!> test_splprep
program test_splprep
  USE numfor, only: dp, Zero, M_PI, str, linspace, save_array
  USE numfor, only: UnivSpline, splprep, splev
  implicit none

  call test_parametric_splines()
contains
  subroutine test_parametric_splines()
    implicit none
    real(dp), dimension(:), allocatable :: r, phi, u
    real(dp), dimension(:, :), allocatable :: x, new_points
    type(UnivSpline) :: tck
    character(len=:), allocatable :: header
    integer :: Nd = 40          ! Number of points
    real(dp) :: s = 0._dp

    allocate (r(Nd), u(Nd), phi(Nd))
    allocate (x(2, Nd), new_points(2, Nd))

    phi = linspace(Zero, 2 * M_PI, Nd)
    r = 0.5_dp + cos(phi)        ! polar coords
    x(1, :) = r * cos(phi)      ! convert to cartesian
    x(2, :) = r * sin(phi)      ! convert to cartesian

    call splprep(x, u, tck, s=s)
    call splev(u, tck, new_points)

    header = "u                    x                   y"
    call save_array([u, new_points(1, :), new_points(2, :)], 3, fmt="f17.14", header=header)

  end subroutine test_parametric_splines
end program test_splprep

