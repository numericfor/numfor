!> test_wfitpack
program test_wfitpack
  USE utils, only: dp, Zero, M_PI, str
  USE grids, only: linspace
  USE strings, only: center
  USE fitpack
  implicit none

  character(len=:), allocatable :: fname
  character(len=:), allocatable :: fdata
  character(len=:), allocatable :: header

  call test_parametric_splines()

  ! ds = 1._dp
  ! header = '#   t           c'
  ! do i = 1, 4
  !   s = (i - 1) * ds
  !   ! print *, center(' s = '//str(s)//' ', 120, '*')
  !   fname = 'data/testsplrep_per_s'//str(s)//'.dat'
  !   call splrep(x, y, tck=tck, s=s, per=.True.)
  !   call save_arrays(fname, tck%t, tck%c, head=header)
  !   ynew = splev(xnew, tck)
  ! end do

contains
  !> test_splines1
  !!
  !! Examples:
  !!
  subroutine test_splines1()
    implicit none

    real(dp), dimension(:), allocatable :: x
    real(dp), dimension(:), allocatable :: y

    integer, parameter :: N = 6
    integer, parameter :: Nnew = 59
    real(dp), dimension(Nnew) :: xnew, ynew, yd1, yd2
    integer :: i
    real(dp) :: s, ds

    type(UnivSpline) :: tck

    allocate (x(N), y(N))
    x = linspace(Zero, M_PI, N)
    xnew = linspace(-0.6_dp, M_PI, Nnew)
    y = sin(x)
    ds = 1._dp

    ! Check x,y curves
    do i = 1, 1
      s = (i - 1) * ds
      ! print *, center(' s = '//str(s)//' ', 120, '*')
      fname = 'data/testsplrep_s'//str(s)//'.dat'
      call splrep(x, y, tck=tck, s=s)
      ynew = splev(xnew, tck)
      yd1 = splev(xnew, tck, der=1)

      yd2 = splev(xnew, tck, der=1, ext=3)
      ! print *, ynew
      header = '#   t           c'
      call save_arrays(fname, tck%t, tck%c, head=header)
      header = "#   x           y"
      fdata = 'data/fosplev_s'//str(s)//'.dat'
      call save_arrays(fdata, xnew, ynew, head=header)
      fdata = 'data/fosplev1_s'//str(s)//'.dat'
      call save_arrays(fdata, xnew, yd1, head=header)
      fdata = 'data/fosplev2_s'//str(s)//'.dat'
      call save_arrays(fdata, xnew, yd2, head=header)

    end do
  end subroutine test_splines1
  !> test_parametric_splines
  !!
  !! Examples:
  !!
  subroutine test_parametric_splines()
    implicit none
    real(dp), dimension(:), allocatable :: phi
    real(dp), dimension(:), allocatable :: r
    real(dp), dimension(:, :), allocatable :: x
    ! real(dp), dimension(:, :), allocatable :: new_points
    real(dp), dimension(:), allocatable :: u
    type(UnivSpline) :: tck
    real(dp) :: s = 0._8
    integer :: Nd = 40          ! Number of points

    allocate (r(Nd), u(Nd), phi(Nd), x(2, Nd))
    phi = linspace(Zero, 2.*M_PI, Nd)
    r = 0.5_8 + cos(phi)        ! polar coords
    x(1, :) = r * cos(phi)      ! convert to cartesian
    x(2, :) = r * sin(phi)      ! convert to cartesian

    call splprep(x, u, tck, s=s)
    ! new_points = splev(u, tck)
    print *, "OUT"
    ! print *, tck%t(:6)

    ! fdata = 'data/fopsplev_s'//str(s)//'.dat'
    ! header = "#   x           y"
    ! call save_arrays(fdata, tck%t, tck%c, head=header)

  end subroutine test_parametric_splines
  !> save_arrays
  !!
  !! Examples:
  !!
  subroutine save_arrays(filename, x, y, head)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !<
    real(dp), dimension(size(x)), intent(IN) :: y !<
    character(len=:), allocatable, intent(IN) :: filename
    character(len=:), allocatable, intent(IN) :: head
    integer :: i, u
    open (newunit=u, file=filename)
    write (u, "(A)") head
    do i = 1, size(x)
      write (u, "(2(g0.14, 1x))") x(i), y(i)
    end do
    close (u)

  end subroutine save_arrays
end program test_wfitpack
