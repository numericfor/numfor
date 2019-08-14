!> polynomials provides a framework for simple (and quite naive) work with polynomials
!! It allows to easily evaluate, derivate, and integrate a polynomial
module polynomial

  USE basic, only: dp, Zero, print_msg

  private
  public polyval, polyder, polyint

  interface polyval
    module procedure polyval_1
    module procedure polyval_v
  end interface polyval

contains
  !> polyval Computes
  !!
  !! Examples:
  !!
  pure function polyval_1(p, x) result(y)
    implicit none
    real(dp) :: y !< Polynomial evaluated in x
    real(dp), dimension(:), intent(IN) :: p !< Array of coefficients, from highest degree to constant term
    real(dp), intent(IN) :: x !< A number at which to evaluate the polynomial
    integer :: i

    y = Zero
    do i = 1, size(p)
      y = y * x + p(i)
    end do

  end function polyval_1

  pure function polyval_v(p, x) result(y)
    implicit none
    real(dp), dimension(:), intent(IN) :: p !< Array of coefficients, from highest degree to constant term
    real(dp), dimension(:), intent(IN) :: x !< A number at which to evaluate the polynomial
    real(dp), dimension(size(x)) :: y !< Polynomial evaluated in x
    integer :: i, j

    y = p(1)
    do j = 2, size(x)
      do i = 1, size(p)
        y(j) = y(j) * x(j) + p(i)
      end do
    end do
  end function polyval_v

  !> polyder Computes the derivative of a polynomial. Returns an array with the coefficients
  !!
  !! Examples:
  !!```
  !! real(dp), dimension(5) :: p
  !! p = 1._dp;  p(3) = 2._dp  ! Set a polynomial
  !!
  !! print *, polyval(polyder(p), 1.5_dp)     ! First derivative evaluated in x = 1.5
  !! print *, polyval(polyder(p), 4, 0._dp)   ! Fourth derivative evaluated in x = 0
  !!
  function polyder(p, m) result(Pd)
    implicit none
    real(dp), dimension(:), intent(IN) :: p !<
    integer, optional, intent(IN) :: m !<
    real(dp), dimension(:), allocatable :: Pd !<
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
  !!
  !! Examples:
  !!
  function polyint(p, m, k) result(p_I)
    implicit none
    real(dp), dimension(:), intent(IN) :: p !<
    integer, optional, intent(IN) :: m !<
    real(dp), optional, intent(IN) :: k !<
    real(dp), dimension(:), allocatable :: p_I !<
    integer :: m_
    real(dp) :: k_
    integer :: i, j, n
    integer :: order

    m_ = 1; IF (present(m)) m_ = m
    IF (m_ < 0) call print_msg('Order of derivative must be positive', errcode=0)

    k_ = Zero; IF (present(k)) k_ = k

    if (m_ == 0) then           ! Return the original polynomial
      p_I = p
      return
    end if

    order = size(p) + m_      ! Number of term of resulting polynomial
    allocate (p_I(order)); p_I(:order) = p; p_I(order + 1:) = k_

    morder: do j = 0, m_ - 1
      n = size(p) - j
      do i = 1, order
        P_I(i) = P_I(i) / (n + i)
      end do
    end do morder
  end function polyint

end module polynomial
