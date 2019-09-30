!> test_wfitpack
program test_wfitpack
  USE utils, only: dp, Zero, M_PI, str
  USE grids, only: linspace
  USE strings, only: center
  USE fitpack
  implicit none
  real(dp), dimension(:), allocatable :: x
  real(dp), dimension(:), allocatable :: y

  integer, parameter :: N = 6
  integer, parameter :: Nnew = 59
  real(dp), dimension(Nnew) :: xnew, ynew, yd1, yd2
  integer :: i
  real(dp) :: s, ds
  character(len=:), allocatable :: fname
  character(len=:), allocatable :: fdata
  character(len=:), allocatable :: header

  type(UnivSpline) :: tck

  allocate (x(N), y(N))
  x = linspace(Zero, M_PI, N)
  xnew = linspace(-0.6_dp, M_PI, Nnew)
  y = sin(x)
  ds = 1._dp

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
