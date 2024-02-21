!> @file polynomial.f90
!! @date "2024-02-21 15:00:50"

!> @ingroup interpolate
!> polynomials provides a framework for simple (and quite naive) work with polynomials
!! Further description in @ref submodule-interpolate
!!
!! It allows to easily evaluate, derivate, and integrate a polynomial
!! Examples:
!! -------
!! @include ex_polynomial.f90
!!
module polynomial

  USE basic, only: dp, Zero, print_msg

  !> Computes the value of the polynomial when applied to a number or list of numbers
  !!
  !! Examples:
  !! --------
  !! @snippet ex_polynomial.f90 evaluate
  interface polyval
    module procedure :: polyval_1, polyval_v
  end interface polyval

  private
  public polyval, polyder, polyint, bisect_pol

contains

  !> Evaluation of a polynomial in a number
  pure function polyval_1(p, x) result(y)
    implicit none
    real(dp) :: y               !< Value of polynomial evaluated in x
    real(dp), dimension(:), intent(IN) :: p !< Array of coefficients, from highest degree to constant term
    real(dp), intent(IN) :: x !< A number at which to evaluate the polynomial
    integer :: i

    y = p(1)
    do i = 2, size(p)
      y = y * x + p(i)
    end do

  end function polyval_1

  !> Evaluation of a polynomial in an array of numbers
  pure function polyval_v(p, x) result(y)
    implicit none
    real(dp), dimension(:), intent(IN) :: p !< Array of coefficients, from highest degree to constant term
    real(dp), dimension(:), intent(IN) :: x !< A number at which to evaluate the polynomial
    real(dp), dimension(size(x)) :: y !< Polynomial evaluated in x
    !! Examples:
    !! --------
    !! @snippet ex_polynomial.f90 derivative

    integer :: i, j

    y = p(1)
    do j = 1, size(x)
      do i = 2, size(p)
        y(j) = y(j) * x(j) + p(i)
      end do
    end do
  end function polyval_v

  !> polyder Computes the derivative of a polynomial. Returns an array with the coefficients
  function polyder(p, m) result(Pd)
    implicit none
    real(dp), dimension(:), intent(IN) :: p !< Array of coefficients, from highest degree to constant term
    integer, optional, intent(IN) :: m !< Order of derivation
    real(dp), dimension(:), allocatable :: Pd !< Derivative of polynomial
    !! Examples:
    !! --------
    !! @snippet ex_polynomial.f90 derivative
    integer :: m_
    integer :: i, k, n
    integer :: order

    m_ = 1; IF (present(m)) m_ = m
    IF (m_ < 0) call print_msg('Order of derivative must be positive', errcode=0)

    if (m_ == 0) then           ! Return the original polynomial
      Pd = p
      return
    end if

    if (size(p) - m_ < 0) then  ! Return the null polynomial
      allocate (Pd(1))
      Pd = Zero
      return
    end if

    order = size(p) - m_      ! Number of term of resulting polynomial
    allocate (Pd(order)); Pd = p(:order)

    morder: do k = 0, m_ - 1
      n = size(p) - k
      do i = 1, order
        Pd(i) = (n - i) * Pd(i)
      end do
    end do morder
  end function polyder

  !> polyint Computes m-esima antiderivative
  function polyint(p, m, k) result(p_I)
    implicit none
    real(dp), dimension(:), intent(IN) :: p !< Array of coefficients, from highest degree to constant term
    integer, optional, intent(IN) :: m !< Number of times that `p` must be integrated
    real(dp), optional, intent(IN) :: k !< Additive Constant
    real(dp), dimension(:), allocatable :: P_I !< Antiderivative polynomial
    !! Examples:
    !! --------
    !! @snippet ex_polynomial.f90 integrate
    integer :: m_
    real(dp) :: k_
    integer :: i, j, n
    integer :: order

    m_ = 1; IF (present(m)) m_ = m
    IF (m_ < 0) call print_msg('Order of antiderivative must be positive', errcode=0)

    k_ = Zero; IF (present(k)) k_ = k

    if (m_ == 0) then           ! Return the original polynomial
      p_I = p
      return
    end if

    order = size(p) + m_      ! Number of term of resulting polynomial
    allocate (p_I(order)); p_I(:size(p)) = p; p_I(size(p) + 1:) = k_

    morder: do j = 1, m_
      n = size(p) + j
      do i = 1, n - 1
        P_I(i) = P_I(i) / (n - i)
      end do
    end do morder
  end function polyint

  !> bisect_pol Classical bisection method for root finding on polynomials
  subroutine bisect_pol(x0, dx, p, toler, x)
    implicit none
    real(dp), intent(IN) :: x0 !< Initial value
    real(dp), intent(INOUT) :: dx !< range. It will probe in the range (x0-dx, x0+dx). On return it will have an estimation of error
    real(dp), dimension(:), intent(IN) :: p !< Array with coefficients of polynomial
    real(dp), intent(IN) :: toler           !< Tolerance in the root determination
    real(dp), intent(OUT) :: x              !< Value of the root

    real(dp) :: y0, yl, yu

    x = x0
    y0 = polyval(p, x)
    yl = polyval(p, x - dx)
    yu = polyval(p, x + dx)
    IF ((yl * y0 >= 0._dp) .and. (yu * y0 >= 0._dp)) return

    do while (dx > toler .and. abs(y0) > toler * 1.e-4)
      ! do while (dx > toler)
      dx = dx / 2
      if (yu * y0 < 0) then
        x = x + dx
      else
        x = x - dx
        yu = y0
      end if
      y0 = polyval(p, x)
    end do
  end subroutine bisect_pol

end module polynomial
! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:
