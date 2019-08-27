!> @file array_utils.f90
!! @date "2019-08-27 16:17:37"

!> This module provides convenience routines to operate or get information on arrays

module array_utils
  use basic, only: dp, Zero, Small, stdout, print_msg

  implicit none
  real(dp), parameter :: def_base = 10._dp

  private
  public :: savetxt, mean, std, merge_arrays

contains

  !> merge_sorted Creates a sorted array with the value from two sorted arrays
  !!
  !! Equal values are only added once
  !! Examples:
  !!
  function merge_sorted(x1, x2, tolerance) result(y)
    implicit none
    real(dp), dimension(:), intent(IN) :: x1 !< First array
    real(dp), dimension(:), intent(IN) :: x2 !< Second array
    real(dp), optional, intent(IN) :: tolerance
    real(dp), dimension(size(x1) + size(x2)) :: y !< Output array with values from both x1 and x2
    real(dp) :: tol_
    integer :: i1               ! Used to visit the first array
    integer :: i2               ! Used to visit the second array
    integer :: iy               ! Used to visit the output array
    real(dp), dimension(:), pointer :: p1, p2

    tol_ = 2 * Small; IF (present(tolerance)) tol_ = tolerance
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
    n = min(n1, n2)

    ! Now p1 is "smaller" than p2 (should exhaust first)
    i1 = 1; i2 = 1

    ! First element
    if (p1(i1) < p2(i2)) then
      y(iy) = p1(i1)
      i1 = 2
    else
      y(iy) = p2(i2)
      i2 = 2
    end if
    iy = 2

    do
      if (p1(i1) < p2(i2)) then
        if (p1(i1) - y(iy) > tol_) then
          y(iy) = p1(i1); iy = iy + 1
        end if
        i1 = i1 + 1
      else
        if (p2(i2) - y(iy) > tol_) then
          y(iy) = p2(i2); iy = iy + 1
        end if
        i2 = i2 + 1
      end if
      IF (i1 > n) exit
    end do

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

  !> savetxt Guarda un array 2D en un archivo de texto
  !!
  !! @note
  !! Si fname es "stdout" o " ", o no están presente ni fname ni unit,
  !! usa stdout
  !! Si se da fname el archivo se abre y cierra.
  !! Si se da unit, el archivo queda abierto
  subroutine savetxt(a, fmt, fname, unit)
    implicit none
    real(dp), dimension(:, :), intent(IN) :: a        !< Array a
    !escribir a archivo de texto
    character(len=*), optional, intent(in) :: fmt    !< formato a
    !usar para los datos. Default 'g0.5'
    character(len=*), optional, intent(in) :: fname  !< Nombre del
    !archivo de salida
    integer, optional, intent(in) :: unit            !< Unidad a
    !escribir si el archivo está abierto

    real(dp), dimension(ubound(a, 2), ubound(a, 1)) :: b
    integer, dimension(2) :: sh
    integer :: i
    integer :: u

    character(len=32) :: form = "g0.5" ! Default
    character(len=32) :: formato
    logical :: closef

    ! Si fname está presente => Toma precedencia sobre unit. Si
    ! ninguna está usa stdout
    closef = .False.
    u = stdout
    if (present(fname)) then
      if (trim(fname) /= '' .and. trim(fname) /= 'stdout') then
        open (newunit=u, file=trim(fname))
        closef = .True.
      end if
    else if (present(unit)) then ! The file was already open before
      ! invoking the function
      IF (unit >= 0 .and. unit <= 99) u = unit
    end if

    b = transpose(a)
    sh = shape(b)

    if (present(fmt) .and. (trim(fmt) /= 'default') .and. (trim(fmt) /= '')) then
      if (index('(', fmt) == 0) then
        write (formato, '(A,I1,A,A,A)') '(', sh(1), '(', trim(fmt), '&
          &,1x))'
      else
        formato = fmt
      end if
    else
      write (formato, '(A,I1,A,A,A)') '(', sh(1), '(', trim(form), '&
        &,1x))'
    end if
    do i = 1, sh(2)
      write (u, formato) b(:, i)
    end do

    ! write (u, formato) (b(:, i), i=1, sh(2))

    IF (closef) close (u)
  end subroutine savetxt

end module array_utils

! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:

