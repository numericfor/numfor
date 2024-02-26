!> @file grids.f90
!! @date "2024-02-26 17:16:17"

!> This module provides convenience routines to create grids
!! Description: @ref submodule-arrays
!!
module grids
  use basic, only: dp, Zero, Small, stdout, print_msg

  implicit none
  real(dp), parameter :: def_base = 10._dp

  private
  public :: linspace, logspace, geomspace, loglinspace, arange

contains

  !! Returns `num` evenly spaced samples, calculated over the
  !! interval [`start`, `stop`].
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
    real(dp), dimension(num) :: x !< An array of uniformly spaced numbers

    !> Examples:
    !! --------
    !!  @snippet ex_linspace.f90 linspace

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

    do concurrent(i=1:num - 1)
      x(i + 1) = start + step * i
    end do

    ! We make sure that the last point is the one desired
    IF (endpoint_ .and. (num > 1)) x(num) = end
  end function linspace

  !> Makes a grid with numbers spaced evenly on a log scale
  !!
  !! In linear space, the sequence starts at ``base**start``
  !! (`base` to the power of `start`) and ends with ``base**end``
  !!
  function logspace(start, end, num, endpoint, base) result(x)
    implicit none
    real(dp), intent(IN) :: start !< ``base**start`` is the starting value of the sequence.
    real(dp), intent(IN) :: end   !< ``base**end`` is the final value of the sequence.
    integer, intent(IN) :: num !< Number of samples to generate. Must be positive.
    logical, optional, intent(IN) :: endpoint !< If True, `end` is the last sample. Otherwise, it is not included. Default is True
    real(dp), optional, intent(IN) :: base    !< The base of the log space. Default is 10.
    real(dp), dimension(num) :: x !< A sequence of numbers spaced evenly on a log scale.

    !> Examples:
    !! --------
    !! @snippet ex_logspace.f90 logspace

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
    real(dp), dimension(num) :: x !< A sequence of numbers spaced evenly on a log scale.

    !> Examples:
    !! --------
    !! @snippet ex_logspace.f90 geomspace

    real(dp) :: sgout
    real(dp) :: lstart, lstop

    IF (num < 1) return

    IF (start * end <= Zero)&
      & call print_msg('Geometric sequence cannot include zero', 'geomspace')

    x(1) = start
    IF (num == 1) return

    lstart = log10(abs(start)); lstop = log10(abs(end))
    sgout = 1._dp; IF ((start < Zero) .and. (end < Zero)) sgout = -1._dp

    x = sgout * logspace(lstart, lstop, num=num, endpoint=endpoint, base=10._dp)

  end function geomspace

  !> loglinspace Computes a grid that may behave as linearly or logarithmically spaced
  !!
  !! @details From package RADIAL by Salvat et al 1995 Computer Physics Communications.
  !! The grid is such that:
  !!   -# R(1)=0, R(NP)=RN,
  !!   -# A*R(I)+B*DLOG(R(I))-C= I  (I > 0),
  !!    with  A=1.0/STEP  and   B=1.0/DLOG(RATIO).
  function loglinspace(start, end, num, step, ratio) result(x)
    implicit none
    real(dp), intent(IN) :: start !< Starting value
    real(dp), intent(IN) :: end !< Final value of the sequence
    integer, intent(IN) :: num !< Number of points
    real(dp), optional, intent(IN) :: step !< Approximated step in the linear region
    real(dp), optional, intent(IN) :: ratio !< quotient between consecutive points in the logarithmic region
    real(dp), dimension(num) :: x !< Array of points with the grid of dimension `num`
    !!
    !! Examples:
    !! -------

    real(dp) :: a, b, c
    real(dp) :: rr, ru, rl, r_N, c_i, fu, fr
    real(dp) :: step_, ratio_
    integer :: i
    real(dp), parameter :: accuracy = 1.e-10_dp

    r_N = end - start           ! Define the interval

    ! Default values
    step_ = r_N / (5 * num / 6._dp); IF (present(step)) step_ = step
    ratio_ = 1.15; IF (present(ratio)) ratio_ = ratio

    a = 1.0_dp / step_
    b = 1.0_dp / log(ratio_)
    c = num - a * r_N - b * log(r_N)
    x(1) = Zero

    rr = Small                  ! Initial guess for x(2)
    ! Solve the equation a*x(i) +b* ln(x(i)) + c - i = 0 (the root) by bisection
    do i = 2, num
      c_i = c - i
      rl = rr
      ru = rl
      search: do                     ! Search a point on the right of the root
        ru = 2 * ru
        fu = a * ru + b * log(ru) + c_i
        if (fu >= Zero) exit search
      end do search

      bisection: do
        rr = 0.5d0 * (ru + rl)
        fr = a * rr + b * log(rr) + c_i
        if (fr > Zero) then
          ru = rr
        else
          rl = rr
        end if
        if ((ru - rl) <= accuracy * rr) exit bisection
      end do bisection

      x(i) = rr
    end do

    x = x + start               ! Add the offset

  end function loglinspace

  !> arange: Return evenly spaced integer values within a given interval
  !!
  !! Values are generated within the half-open interval ``[start, end)``
  !! (in other words, the interval including `start` but excluding `end`).
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

    if (end == start) then
      allocate (x(1))
      x = start
      return
    end if

    num = ceiling(abs((end - start) / real(step_, kind=dp)))
    allocate (x(1:num))
    x(1) = start
    do concurrent(i=1:num)
      x(i + 1) = start + i * step_
    end do
  end function arange

end module grids
! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:

