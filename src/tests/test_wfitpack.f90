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
  integer :: i
  real(dp) :: s, ds
  character(len=:), allocatable :: fname
  character(len=:), allocatable :: header

  type(UnivSpline) :: tck

  allocate (x(N), y(N))
  x = linspace(Zero, M_PI, N)
  y = sin(x)
  ds = 1._dp
  header = '#   t           c'
  do i = 1, 4
    s = (i - 1) * ds
    ! print *, center(' s = '//str(s)//' ', 120, '*')
    fname = 'data/testsplrep_s'//str(s)//'.dat'
    call splrep(x, y, tck=tck, s=s)
    call save_arrays(fname, tck%t, tck%c, head=header)
  end do
  ! call splrep(x, y, tck=tck, ier=ier)
contains
  !> save_arrays
  !!
  !! Examples:
  !!
  subroutine save_arrays(fname, x, y, head)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !<
    real(dp), dimension(size(x)), intent(IN) :: y !<
    character(len=:), allocatable, intent(IN) :: fname
    character(len=:), allocatable, intent(IN) :: head
    integer :: i, u
    open (newunit=u, file=fname)
    write (u, "(A)") head
    do i = 1, size(x)
      write (u, "(2(g0.14, 1x))") x(i), y(i)
    end do
    close (u)

  end subroutine save_arrays
end program test_wfitpack
