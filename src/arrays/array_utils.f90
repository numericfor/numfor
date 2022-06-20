!> @file array_utils.f90
!! @date "2022-06-20 11:40:10"

!> This module provides convenience routines to operate or get information on arrays

module array_utils
  USE basic, only: dp, Zero, Small, stdout, print_msg
  USE strings, only: str

  implicit none

  Private
  PUBLIC :: allclose, save_array, mean, std, merge_sorted
  ! PUBLIC :: savetxt

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

  ! !> savetxt Guarda un array 2D en un archivo de texto
  ! !!
  ! !! @note
  ! !! Si fname es "stdout" o " ", o no están presente ni fname ni unit,
  ! !! usa stdout
  ! !! Si se da fname el archivo se abre y cierra.
  ! !! Si se da unit, el archivo queda abierto
  ! subroutine savetxt(X, fname, fmt, unit)
  !   implicit none
  !   real(dp), dimension(:, :), intent(IN) :: X        !< Array a escribir a archivo de texto
  !   character(len=*), optional, intent(IN) :: fname  !< Nombre del archivo de salida
  !   character(len=*), optional, intent(IN) :: fmt    !< formato a usar para los datos. Default 'g0.5'
  !   integer, optional, intent(IN) :: unit            !< Unidad a escribir si el archivo está abierto

  !   real(dp), dimension(ubound(X, 2), ubound(X, 1)) :: b
  !   integer, dimension(2) :: sh
  !   integer :: i
  !   integer :: u

  !   character(len=32) :: form = "g0.5" ! Default
  !   character(len=32) :: formato
  !   logical :: closef

  !   ! Si fname está presente => Toma precedencia sobre unit. Si
  !   ! ninguna está usa stdout
  !   closef = .False.
  !   u = stdout

  !   if (Present(fname)) then
  !     if (trim(fname) /= '' .and. trim(fname) /= 'stdout') then
  !       open (newunit=u, file=trim(fname))
  !       closef = .True.
  !     end if
  !   else if (present(unit)) then ! The file was already open before
  !     ! invoking the function
  !     IF (unit >= 0 .and. unit <= 99) u = unit
  !   end if

  !   b = transpose(X)
  !   sh = shape(b)

  !   if (present(fmt) .and. (trim(fmt) /= 'default') .and. (trim(fmt) /= '')) then
  !     if (index('(', fmt) == 0) then
  !       write (formato, '(A,I1,A,A,A)') '(', sh(1), '(', trim(fmt), '&
  !         &,1x))'
  !     else
  !       formato = fmt
  !     end if
  !   else
  !     write (formato, '(A,I1,A,A,A)') '(', sh(1), '(', trim(form), '&
  !       &,1x))'
  !   end if

  !   do i = 1, sh(2)
  !     write (u, formato) b(:, i)
  !   end do

  !   IF (closef) close (u)
  ! end subroutine savetxt

  !> save_array Stores an array to file or stdout
  !!
  !! Examples:
  !!```
  !!  real(dp), dimension(20), allocatable :: x,y
  !!  character(len=:), allocatable :: filename
  !!  filename = "output.dat"
  !!  x = linspace(0, 10, 20)
  !!  y = -x**2/10
  !!  save_array(x)  ! One array in one column to stdout
  !!  save_array(x, 1, filename) ! One array in one column to file
  !!  save_array([x,y], 2, filename) ! One array in one column to file
  !!```
  subroutine save_array(X, ncols, fname, fmt, header, unit)
    implicit none
    real(dp), dimension(:), intent(IN) :: X        !< Array to store
    integer, optional, intent(IN) :: ncols !< Number of columns to write. Default 1

    character(len=*), optional, intent(IN) :: fname  !< Filename. Default: stdout
    character(len=*), optional, intent(IN) :: fmt    !< String with format. Default 'g0.5' for each data
    character(len=*), optional, intent(IN) :: header  !< If present, text to write before data.
    integer, optional, intent(IN) :: unit            !< If the file is already open, the unit at which it is associated.

    real(dp), dimension(:, :), allocatable :: b
    integer :: ncols_, nrows_
    integer :: i
    integer :: u

    character(len=32) :: form = "g0.5" ! Default
    character(len=32) :: formato
    logical :: closef

    ncols_ = 1; IF (Present(ncols)) ncols_ = ncols

    nrows_ = size(X) / ncols_
    IF (nrows_ * ncols_ /= size(X)) call print_msg("Incompatible&
      & number of columns"//str(ncols_)//" and size(X)="//str(size(X))&
      &//" for saving", "save_arrays", errcode=1)

    ! Si fname está presente => Toma precedencia sobre unit.
    closef = .False.
    u = stdout    ! Si no está ni fname ni unit está usa stdout

    if (Present(fname)) then
      if (trim(fname) /= '' .and. trim(fname) /= 'stdout') then
        open (newunit=u, file=trim(fname))
        closef = .True.
      end if
    else if (Present(unit)) then ! The file was already open before
      ! invoking the function
      IF (unit >= 0 .and. unit <= 99) u = unit
    end if

    b = reshape(X, [ncols_, nrows_], order=[2, 1])

    if (present(fmt) .and. (trim(fmt) /= 'default') .and. (trim(fmt) /= '')) then
      if (index(fmt, '(') == 0) then
        write (formato, '(A,I1,A,A,A)') '(', ncols_, '(', trim(fmt), '&
          &,1x))'
      else
        formato = fmt
      end if
    else
      write (formato, '(A,I1,A,A,A)') '(', ncols_, '(', trim(form), '&
        &,1x))'
    end if

    if (Present(header)) write (u, '(A)') "# "//trim(header)
    do i = 1, nrows_
      write (u, formato) b(:, i)
    end do

    IF (closef) close (u)
  end subroutine save_array
end module array_utils

! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:

