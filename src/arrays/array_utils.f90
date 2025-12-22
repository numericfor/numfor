!> @file array_utils.f90
!! @date "2025-12-22 02:52:30"

!> This module provides convenience routines to operate or get information on arrays

module array_utils
  USE basic, only: dp, Zero, Small, stdout, print_msg
  USE strings, only: str

  implicit none

  character(len=1), parameter :: delimiter_default = " "
  character(len=1), parameter :: comment_default = "#"

  Private
  PUBLIC :: allclose, save_array, mean, std, merge_sorted
  ! PUBLIC :: savetxt

  !> save_array Stores an 1D or 2D array to file or stdout
  !!
  !! Examples:
  !!```
  !!  real(dp), dimension(20), allocatable :: x,y
  !!  real(dp), dimension(20,2) :: z
  !!  character(len=:), allocatable :: filename
  !!  filename = "output.dat"
  !!  x = linspace(0, 10, 20)
  !!  y = -x**2/10
  !!  save_array(x)  ! One array in one column to stdout
  !!  save_array(x, 1, filename) ! One array in one column to file
  !!  save_array([x,y], 2, filename) ! Two arrays in two columns to file
  !!  save_array([x,y], 2, fmt='g12.5', fout=filename) ! Two arrays in two columns to file with format
  !!  save_array(z)  ! A 2D array in two columns to stdout
  !!  save_array(z, fmt='g12.5', fout=filename) ! A 2D array in two columns to file with format
  !!```
  interface save_array
    module procedure save_array1D
    module procedure save_array2D
  end interface save_array

contains

  !> allclose returns True if two arrays are element-wise equal within a tolerance.
  !!
  !! Very similar to **Numpy** allclose
  !!
  !! The tolerance values are positive, typically very small numbers.  The
  !! relative difference (`rtol` * abs(`b`)) and the absolute difference
  !! `atol` are added together to compare against the absolute difference
  !! between `a` and `b`.\n
  !! If the following equation \n
  !!  abs(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`)) \n
  !!is element-wise True, then allclose returns `.True.`
  !!
  !! The above equation is not symmetric in `a` and `b`, so that
  !! ``allclose(a, b)`` might be different from ``allclose(b, a)`` in
  !! some rare cases.
  !!
  !! Examples:
  !! --------
  !! ```
  !! allclose([1e10,1e-7], [1.00001e10,1e-8])
  !! ! .False.
  !! allclose([1e10,1e-8], [1.00001e10,1e-9])
  !! ! .True.
  !! allclose([1e10,1e-8], [1.0001e10,1e-9])
  !! ! .False.
  !! !
  !! ```
  !!
  function allclose(a, b, rtol, atol) result(y)
    implicit none
    logical :: y !< True if the two arrays are equal within the given
    !!    tolerance; False otherwise.
    real(dp), dimension(:), intent(IN) :: a !< Array
    real(dp), dimension(size(a)), intent(IN) :: b !< Array
    real(dp), optional, intent(IN) :: rtol !< The relative tolerance parameter. Default = 1.e-5
    real(dp), optional, intent(IN) :: atol !< The absolute tolerance parameter. Default = 1.e-8
    real(dp) :: rtol_, atol_ !<
    atol_ = 1.e-8_dp; IF (Present(atol)) atol_ = atol
    rtol_ = 1.e-5_dp; IF (Present(atol)) rtol_ = rtol
    y = all(abs(a - b) < (atol_ + rtol_ * abs(b)))
  end function allclose

  !> This function creates a sorted array with values from two input sorted arrays
  !!
  !! Equal values (within tolerance) are only included once
  !!
  !! ### Examples:###
  !!```
  !!  real(dp), dimension(:), allocatable :: a
  !!  real(dp), dimension(:), allocatable :: b
  !!  real(dp), dimension(:), allocatable :: c
  !!  integer :: j
  !!  a = [1.9999999999_dp, 2._dp, 2._dp, 3._dp]
  !!  b = [-1._dp, 2._dp, 3._dp, 4._dp]
  !!  print *, ""
  !!  print '(A)', repeat('-', 30)//' a '//repeat('-', 30)
  !!  print "(4(g0.12,1x))", a
  !!  print '(A)', repeat('-', 30)//' b '//repeat('-', 30)
  !!  print "(4(g0.12,1x))", b
  !!  c = merge_sorted(a, b)
  !!  j = 30 - int(len("merge_sorted(a,b)") / 2._dp)
  !!  print '(A)', repeat('-', 60)
  !!  print '(A)', repeat('-', j)//' merge_sorted(a,b) '//repeat('-', j)
  !!  print "(4(g0.12,1x))", c
  !!  c = merge_sorted(a, b, 1.e-6_dp)
  !!  j = 30 - int(len("merge_sorted(a,b,1.e-6)") / 2._dp)
  !!  print '(A)', repeat('-', j)//' merge_sorted(a,b,1.e-6) '//repeat('-', j)
  !!  print "(4(g0.12,1x))", c
  !!  c = merge_sorted(a, b, -1.e-6_dp)
  !!  j = 30 - int(len("merge_sorted(a,b,-1.e-6)") / 2._dp)
  !!  print '(A)', repeat('-', j)//' merge_sorted(a,b,-1.e-6) '//repeat('-', j)
  !!  print "(4(g0.12,1x))", c
  !!  !
  !!  ! Outputs:
  !!  !
  !!  !  ------------------------------ a ------------------------------
  !!  !  1.99999999990 2.00000000000 2.00000000000 3.00000000000
  !!  !  ------------------------------ b ------------------------------
  !!  !  -1.00000000000 2.00000000000 3.00000000000 4.00000000000
  !!  !  ------------------------------------------------------------
  !!  !  ---------------------- merge_sorted(a,b) ----------------------
  !!  !  -1.00000000000 1.99999999990 2.00000000000 3.00000000000
  !!  !  4.00000000000
  !!  !  ------------------- merge_sorted(a,b,1.e-6) -------------------
  !!  !  -1.00000000000 1.99999999990 3.00000000000 4.00000000000
  !!  !  ------------------ merge_sorted(a,b,-1.e-6) ------------------
  !!  !  -1.00000000000 1.99999999990 2.00000000000 2.00000000000
  !!  !  2.00000000000 3.00000000000 3.00000000000 4.00000000000
  !!  !
  !! ```
  function merge_sorted(x1, x2, tolerance) result(y)
    implicit none
    real(dp), dimension(:), target, intent(IN) :: x1 !< First array
    real(dp), dimension(:), target, intent(IN) :: x2 !< Second array
    real(dp), optional, intent(IN) :: tolerance !< Defines the minimum value by which two numbers are considered different
    real(dp), dimension(:), allocatable :: y !< Output array with values from both x1 and x2
    real(dp), dimension(size(x1) + size(x2)) :: xo !< Output array with values from both x1 and x2
    real(dp) :: tol_
    integer :: i1               ! Used to visit the first array
    integer :: i2               ! Used to visit the second array
    integer :: io               ! Used to visit the output array
    integer :: n1, n2
    real(dp), dimension(:), pointer :: p1, p2

    tol_ = Small; IF (present(tolerance)) tol_ = tolerance
    n1 = size(x1); n2 = size(x2)
    if (x1(n1) <= x2(n2)) then  ! We choose the "smaller"
      p1 => x1
      p2 => x2
    else
      p1 => x2
      p2 => x1
    end if
    n1 = size(p1)
    n2 = size(p2)

    ! Now p1 is "smaller" than p2 (should exhaust first)
    i1 = 1; i2 = 1; io = 1

    ! First element
    if (p1(i1) < p2(i2)) then
      xo(io) = p1(i1)
      i1 = 2
    else
      xo(io) = p2(i2)
      i2 = 2
    end if
    io = 2

    smaller: do                 ! Loop over "smaller" array
      if (p1(i1) < p2(i2)) then
        if (p1(i1) - xo(io - 1) > tol_) then ! Only add them if different within tolerance
          xo(io) = p1(i1); io = io + 1
        end if
        i1 = i1 + 1
      else
        if (p2(i2) - xo(io - 1) > tol_) then
          xo(io) = p2(i2); io = io + 1
        end if
        i2 = i2 + 1
      end if
      IF (i1 > n1) exit smaller
    end do smaller

    remaining: do while (i2 <= n2)
      if (p2(i2) - xo(io - 1) > tol_) then
        xo(io) = p2(i2); io = io + 1
      end if
      i2 = i2 + 1
    end do remaining
    ! allocate (y(io - 1))
    y = xo(:io - 1)                ! Automatic allocation
  end function merge_sorted

  !> std Computes the standard deviation of the array.
  !!
  !! @note : Basically: `sqrt(mean(x - mean(x))* alfa )` with `alfa= (N/(N-1))`
  function std(x) result(y)
    implicit none
    real(dp) :: y !< Standard deviation
    real(dp), dimension(:), intent(IN) :: x !< Input array of real values
    integer :: N
    N = size(x)
    y = sqrt(mean((x - mean(x))**2) * (N / real(N - 1, kind=dp)))
  end function std

  !> mean Computes the arithmetic mean of the array.
  !!
  !! @note the mean is basically: `sum(x)/size(x)`
  function mean(x) result(y)
    implicit none
    real(dp) :: y !< Mean value
    real(dp), dimension(:), intent(IN) :: x !< Input array of real values
    y = sum(x) / size(x)
  end function mean

  ! save_array1D Stores an array to file or stdout
  subroutine save_array1D(X, ncols, fname, delimiter, fmt, header, footer, unit, comment)
    implicit none
    real(dp), dimension(:), intent(IN) :: X        !< Array to store
    integer, optional, intent(IN) :: ncols !< Number of columns to write. Default 1

    character(len=*), optional, intent(IN) :: fname  !< Filename. Default: stdout
    character(len=*), optional, intent(IN) :: fmt    !< String with format. Default 'g0.5' for each data
    character(len=1), optional, intent(IN) :: delimiter  !< Separator between columns. Default " "
    character(len=*), optional, intent(IN) :: header  !< If present, text to write before data.
    character(len=1), optional, intent(IN) :: comment  !< Separator between columns. Default " "
    character(len=*), optional, intent(IN) :: footer  !< If present, text to write after data.
    integer, optional, intent(IN) :: unit            !< If the file is already open, the unit at which it is associated.

    real(dp), dimension(:, :), allocatable :: b
    integer :: ncols_, nrows_
    integer :: unit_

    character(len=*), parameter :: default_fmt = "g0.5" ! Default
    character(len=:), allocatable :: fmt_
    character(len=:), allocatable :: fname_
    character(len=:), allocatable :: delimiter_
    character(len=:), allocatable :: comment_
    character(len=:), allocatable :: header_
    character(len=:), allocatable :: footer_

    ncols_ = 1; IF (Present(ncols)) ncols_ = ncols

    nrows_ = size(X) / ncols_
    IF (nrows_ * ncols_ /= size(X)) call print_msg("Incompatible&
    & number of columns"//str(ncols_)//" and size(X)="//str(size(X))&
    &//" for saving", "save_arrays", errcode=1)

    b = reshape(X, [ncols_, nrows_], order=[2, 1])

    IF (Present(fname) .and. Present(unit)) call print_msg("Only one of fname or unit &
    & may be present", "save_array", errcode=1)

    ! Check delimiter and comment
    delimiter_ = '"'//delimiter_default//'"'
    IF (Present(delimiter)) delimiter_ = "'"//delimiter//"'"
    comment_ = comment_default; IF (Present(comment)) comment_ = comment

    ! Check presence of filename
    fname_ = ''
    IF (Present(fname)) fname_ = trim(fname)

    ! Check presence of unit
    unit_ = stdout; if (Present(unit)) unit_ = unit

    ! Check presence of format
    fmt_ = default_fmt; IF (Present(fmt)) fmt_ = trim(fmt)

    ! Check presence of header and footer
    header_ = ''; if (Present(header)) header_ = trim(header)
    footer_ = ''; if (Present(footer)) footer_ = trim(footer)

    IF ((fname_ /= '')) open (newunit=unit_, file=fname_)

    call save_array(b, delimiter=delimiter_, fmt=fmt_,&
      & header=header_, footer=footer_, unit=unit_, comment=comment_)

    IF (unit_ /= stdout) close (unit_)
  end subroutine save_array1D

  subroutine save_array2D(X, fname, delimiter, fmt, header, footer, unit, comment)
    implicit none
    real(dp), dimension(:, :), intent(IN) :: X        !< Array to save to file
    !
    character(len=*), optional, intent(IN) :: fname  !< Filename. Default: stdout
    character(len=*), optional, intent(IN) :: fmt    !< String with format. Default 'g0.5' for each data
    character(len=1), optional, intent(IN) :: delimiter  !< Separator between columns. Default " "
    character(len=1), optional, intent(IN) :: comment  !< Separator between columns. Default " "
    character(len=*), optional, intent(IN) :: header  !< If present, text to write before data.
    character(len=*), optional, intent(IN) :: footer  !< If present, text to write after data.
    integer, optional, intent(IN) :: unit            !< If the file is already open, the unit at which it is associated.

    integer, dimension(2) :: theshape
    integer :: ncols_, nrows_
    character(len=3) :: delimiter_
    character(len=1) :: comment_

    integer :: i
    integer :: unit_

    character(len=*), parameter :: default_fmt = "g0.5" ! Default
    character(len=:), allocatable :: fmt_
    character(len=:), allocatable :: fname_
    character(len=:), allocatable :: header_
    character(len=:), allocatable :: footer_
    character(len=32) :: formato
    !
    IF (Present(fname) .and. Present(unit)) call print_msg("Only one of fname or unit &
    & may be present", "save_array", errcode=1)

    ! Check delimiter and comment
    delimiter_ = '"'//delimiter_default//'"'
    IF (Present(delimiter)) delimiter_ = "'"//delimiter//"'"
    comment_ = comment_default; IF (Present(comment)) comment_ = comment

    ! Check presence of filename
    fname_ = ''
    IF (Present(fname)) fname_ = trim(fname)

    ! Check presence of unit
    unit_ = stdout; if (Present(unit)) unit_ = unit

    ! Check presence of format
    fmt_ = default_fmt; IF (Present(fmt)) fmt_ = trim(fmt)

    ! Check presence of header and footer
    header_ = ''; if (Present(header)) header_ = trim(header)
    footer_ = ''; if (Present(footer)) footer_ = trim(footer)

    IF ((fname_ /= '')) open (newunit=unit_, file=fname_)
    !
    !
    theshape = shape(X)
    ncols_ = theshape(1)
    nrows_ = theshape(2)
    if (fmt_ /= 'default' .and. (fmt_ /= '')) then
      if (index(fmt_, '(') == 0) then
        write (formato, '(A,I1,A,A,A,A,A)') '(', ncols_, '(', trim(fmt_), ',', delimiter_, '))'
      else
        formato = fmt_
      end if
    else
      write (formato, '(A,I1,A,A,A,A,A)') '(', ncols_, '(', trim(default_fmt), ',', delimiter_, '))'
    end if

    if (header_ /= '') write (unit_, '(A)') comment_//" "//header_
    do i = 1, nrows_
      write (unit_, formato) X(:, i)
    end do
    if (footer_ /= '') write (unit_, '(A)') comment_//" "//footer_

    IF (unit_ /= stdout) close (unit_)
  end subroutine save_array2D
end module array_utils

! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:

