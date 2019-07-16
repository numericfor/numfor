!> This module provides convenience routines to work with grids and arrays
module grids
  use basic, only: dp, Zero, stdout, print_msg

  implicit none
  real(dp), parameter :: def_base = 10._dp

  private
  public :: linspace, logspace, geomspace, arange, searchsorted

contains

  !> Return evenly spaced numbers over a specified interval
  !!
  !! Returns `num` evenly spaced samples, calculated over the
  !! interval [`start`, `stop`].
  !!
  !!  Examples
  !!  ```
  !!  linspace(2.0, 3.0, num=5)
  !!  ! gives:  [ 2.  ,  2.25,  2.5 ,  2.75,  3.  ]
  !!  linspace(2.0, 3.0, num=5, endpoint=.False.)
  !!  ! gives: [ 2. ,  2.2,  2.4,  2.6,  2.8]
  !!  linspace(2.0, 3.0, num=5, retstep=step)
  !!  ! gives: [ 2.  ,  2.25,  2.5 ,  2.75,  3.  ]
  !!  ! and retstep= 0.25
  !!  !
  !! ```
  !!
  function linspace(start, end, num, endpoint, retstep) result(x)
    implicit none
    real(dp), intent(IN) :: start !< The starting value of the
    !sequence.
    real(dp), intent(IN) :: end   !< The end value of the sequence,
    integer, intent(IN) :: num !< Number of samples to generate. Must
    !be positive.
    logical, optional, intent(IN) :: endpoint !< If True, `end` is
    !the last sample. Otherwise, it is not included. Default is True
    real(dp), optional, intent(OUT) :: retstep !< If present, return
    !the step
    real(dp), dimension(num) :: x              !< An array of
    !uniformly spaced numbers
    real(dp) :: step
    integer :: i
    logical :: endpoint_
    x = Zero
    IF (num < 1) return

    x(1) = start
    IF (num == 1) return

    endpoint_ = .True.; IF (present(endpoint)) endpoint_ = endpoint
    if (endpoint_) then
      step = (end - start) / (num - 1._dp)
    else
      step = (end - start) / real(num, kind=dp)
    end if

    IF (present(retstep)) retstep = step

    forall (i=1:num - 1) x(i + 1) = start + step * i
    ! We make sure that the last point is the one desired
    IF (endpoint_ .and. (num > 1)) x(num) = end
  end function linspace

  !> Makes a grid with numbers spaced evenly on a log scale
  !!
  !! In linear space, the sequence starts at ``base**start``
  !! (`base` to the power of `start`) and ends with ``base**end``
  !!
  !! Examples
  !!
  !! ```
  !! logspace(2.0, 3.0, num=4)
  !! ! gives: [  100.        ,   215.443469  ,   464.15888336,  1000.
  !! ]
  !! logspace(2.0, 3.0, num=4, endpoint=False)
  !! ! gives: [ 100.        ,  177.827941  ,  316.22776602,
  !! 562.34132519]
  !! logspace(2.0, 3.0, num=4, base=2.0)
  !! ! gives: array([ 4.        ,  5.0396842 ,  6.34960421,  8.
  !! ])
  !! !
  !! ```
  !!
  function logspace(start, end, num, endpoint, base) result(x)
    implicit none
    real(dp), intent(IN) :: start !< ``base**start`` is the starting
    !value of the sequence.
    real(dp), intent(IN) :: end   !< ``base**end`` is the final value
    !of the sequence.
    integer, intent(IN) :: num !< Number of samples to generate. Must
    !be positive.
    logical, optional, intent(IN) :: endpoint !< If True, `end` is
    !the last sample. Otherwise, it is not included. Default is True
    real(dp), optional, intent(IN) :: base     !< The base of the log
    !space. Default is 10.
    real(dp), dimension(num) :: x              !< A sequence of
    !numbers spaced evenly on a log scale.

    real(dp) :: b_

    IF (num < 1) return

    b_ = def_base; IF (present(base)) b_ = base

    x = b_**(linspace(start, end, num, endpoint))
  end function logspace

  !> Makes a grid with numbers spaced evenly on a log scale
  !!
  !! @note: Is similar to logspace but with endpoints specified
  !! directly.
  !! Also accepts simultaneously negative `start` **and** `end`
  !! Examples
  !!
  !! ```
  !! geomspace(1, 1000.0, num=4)
  !! ! gives: [  1.        ,   10.0  ,   100.0,  1000.  ]
  !! geomspace(-1000, 1.0, num=4)
  !! ! gives: [ -1000. , -100.0  , -10.0 ,  -1.0 ]
  !! !
  !! ```
  function geomspace(start, end, num, endpoint) result(x)
    implicit none
    real(dp), intent(IN) :: start !< ``start`` is the starting value
    !of the sequence.
    real(dp), intent(IN) :: end   !< ``end`` is the final value of
    !the sequence.
    integer, intent(IN) :: num !< Number of samples to generate. Must
    !be positive.
    logical, optional, intent(IN) :: endpoint !< If True, `end` is
    !the last sample. Otherwise, it is not included. Default is True
    real(dp), dimension(num) :: x              !< A sequence of
    !numbers spaced evenly on a log scale.
    real(dp) :: sgout
    real(dp) :: lstart, lstop

    IF (num < 1) return

    IF (start * end <= Zero) call print_msg('Geometric sequence cannot include&
      & zero', 'geomspace')

    x(1) = start
    IF (num == 1) return

    lstart = log10(abs(start)); lstop = log10(abs(end))
    sgout = 1._dp; IF ((start < Zero) .and. (end < Zero)) sgout = &
 &-1._dp

    x = sgout * logspace(lstart, lstop, num=num, endpoint=endpoint,&
      & base=10._dp)

  end function geomspace

  !> arange: Return evenly spaced integer values within a given interval
  !!
  !! Values are generated within the half-open interval ``[start, end)``
  !! (in other words, the interval including `start` but excluding `end`).
  !! @todo function is failing if `step` is not present
  function arange(start, end, step) result(x)
    implicit none
    integer, intent(IN) :: start !< the starting value of the interval.
    integer, intent(IN) :: end   !< the final value of the interval (not included)
    integer, optional, intent(IN) :: step !< Spacing between values.
    integer, dimension(:), allocatable :: x !< A sequence of numbers spaced evenly
    integer :: num, i, step_
    !
    step_ = 1; IF (present(step)) step_ = step
    IF (step_ == 0) call print_msg('Step must be nonzero', 'arange')

    num = ceiling((end - start) / real(step, kind=dp))

    IF (allocated(x) .and. (size(x) /= num)) deallocate (x)
    IF (.not. allocated(x)) allocate (x(num))
    x(1) = start
    do i = 1, num
      x(i + 1) = x(i) + step_
    end do
  end function arange

  !> searchsorted: Find index where an element should be inserted to maintain order.
  !!
  !! Find the index into an ascending sorted array `x` such that, if `elem` was inserted
  !! after the index, the order of `x` would be preserved.
  !!
  !! @note Bisection is used to find the required insertion point
  function searchsorted(x, elem) result(n)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !< Array sorted in ascending order
    real(dp), intent(IN) :: elem !< element to insert
    integer :: n                 !< index of left
    integer :: upper, lower, mid

    lower = 1
    upper = size(x)
    IF (elem < x(lower)) call print_msg('Element below array', sub='searchsorted')
    IF (elem > x(upper)) call print_msg('Element above array', 'searchsorted')

    ! Instead of starting on (1, upper//2 , upper) we start with a linear guess
    n = nint(((elem - x(lower)) / (x(upper) - x(lower))) * upper)
    IF (n < 1) n = 1
    IF (n > size(x)) n = size(x)

    if ((x(n) == elem) .or. ((x(n) <= elem) .and. (x(n + 1) > elem))) then
      return                    ! Found the index
    else if (x(n + 1) == elem) then
      n = n + 1
      return
    else if (x(n) < elem) then
      lower = n
    else
      upper = n
    end if

    do while ((upper - lower > 1))
      mid = (upper + lower) / 2 ! integer arithmetic
      if (elem >= x(mid)) then
        lower = mid
      else
        upper = mid
      endif
    enddo
    n = lower

  end function searchsorted

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
    if (Present(fname)) then
      if (trim(fname) /= '' .and. trim(fname) /= 'stdout') then
        open (newunit=u, file=trim(fname))
        closef = .True.
      end if
    else if (Present(unit)) then ! The file was already open before
      ! invoking the function
      IF (unit >= 0 .and. unit <= 99) u = unit
    end if

    b = transpose(a)
    sh = shape(b)

    if (Present(fmt) .and. (trim(fmt) /= 'default') .and. (trim(fmt) /= '')) then
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

end module grids
