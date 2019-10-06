!> test_wfitpack
program test_wfitpack
  USE utils, only: dp, Zero, M_PI, str
  USE arrays, only: linspace, save_array, savetxt
  USE strings, only: center

  USE fitpack
  implicit none

  real(dp), dimension(1) :: ss

  ! ss = [0._dp, 0.3_dp]
  ss = [0.3_dp]

  print *, "1D Splines"
  call test_splines1(ss)
  print *, "Parametric 1D Splines"
  call test_parametric_splines(ss)

contains
  !> test_splines1
  !!
  !! Examples:
  !!
  subroutine test_splines1(ss)
    implicit none

    real(dp), dimension(:), intent(IN) :: ss
    integer, parameter :: N = 6
    integer, parameter :: Nnew = 59
    real(dp), dimension(N) :: x
    real(dp), dimension(N) :: y
    real(dp), dimension(Nnew) :: xnew
    real(dp), dimension(Nnew) :: ynew
    character(len=:), allocatable :: fname, fofname
    character(len=:), allocatable :: fdata, fofdata

    integer :: i
    real(dp) :: s

    type(UnivSpline) :: tck

    fname = 'data/fosplrep_s'
    fdata = 'data/fosplev_s'

    x = linspace(Zero, M_PI, N)
    y = sin(x)
    xnew = linspace(Zero, M_PI, Nnew)

    ! Check x,y curves
    do i = 1, size(ss)
      s = ss(i)
      fofname = fname//str(s)//'.dat'
      call splrep(x, y, tck=tck, s=s)
      call save_array([tck%t, tck%c], [size(tck%t), 2], fname=fofname)
      call splev(xnew, tck, ynew)
      fofdata = fdata//str(s)//'.dat'
      call save_array([xnew, ynew], [size(xnew), 2], fname=fofdata)
    end do
  end subroutine test_splines1
  !> test_parametric_splines
  !!
  !! Examples:
  !!
  subroutine test_parametric_splines(ss)
    implicit none
    real(dp), dimension(:), intent(IN) :: ss
    real(dp), dimension(:), allocatable :: phi
    real(dp), dimension(:), allocatable :: r
    real(dp), dimension(:, :), allocatable :: x
    real(dp), dimension(:, :), allocatable :: new_points
    real(dp), dimension(:), allocatable :: u
    type(UnivSpline) :: tck
    real(dp) :: s = 0._8
    integer :: Nd = 40          ! Number of points
    integer :: i
    character(len=:), allocatable :: fname, fofname
    character(len=:), allocatable :: fdata, fofdata

    fname = 'data/fosplprep_s'
    fdata = 'data/fosplpev_s'

    allocate (r(Nd), u(Nd), phi(Nd))
    allocate (x(2, Nd), new_points(2, Nd))

    phi = linspace(Zero, 2.*M_PI, Nd)
    r = 0.5_8 + cos(phi)        ! polar coords
    x(1, :) = r * cos(phi)      ! convert to cartesian
    x(2, :) = r * sin(phi)      ! convert to cartesian

    do i = 1, size(ss)
      s = ss(i)
      fofname = fname//str(s)//'.dat'
      call splprep(x, u, tck, s=s)
      call save_array([tck%c, tck%c], [size(tck%c), 2], fname=fofname)

      fofdata = fdata//str(s)//'.dat'
      call splevp(u, tck, new_points)
      call save_array([u, new_points(1, :), new_points(2, :)], [size(u), 3], fofdata)
    end do
  end subroutine test_parametric_splines
end program test_wfitpack
