!> polynomials provides a framework for
module polynomial

  USE basic, only: dp, Zero, print_msg

  private
  public polyval, polyder

  interface polyval
    module procedure polyval_1
    module procedure polyval_v
  end interface polyval

contains
  !> polyval Computes
  !!
  !! Examples:
  !!
  function polyval_1(p, x) result(y)
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

  function polyval_v(p, x) result(y)
    implicit none
    real(dp), dimension(:), intent(IN) :: p !< Array of coefficients, from highest degree to constant term
    real(dp), dimension(:), intent(IN) :: x !< A number at which to evaluate the polynomial
    real(dp), dimension(size(x)) :: y !< Polynomial evaluated in x
    integer :: i, j

    y = Zero
    do j = 1, size(x)
      do i = 1, size(p)
        y(j) = y(j) * x(j) + p(i)
      end do
    end do
  end function polyval_v

  !> polyder Computes
  !!
  !! Examples:
  !!
  function polyder(p, m) result(Pd)
    implicit none
    real(dp), dimension(:), intent(IN) :: p !<
    integer, optional, intent(IN) :: m !<
    real(dp), dimension(:), allocatable :: Pd !<
    integer :: m_
    integer :: i, k, n
    integer :: order, degree

    m_ = 1; IF (present(m)) m_ = m
    IF (m_ < 0) call print_msg('Order of derivative must be positive', errcode=0)

    if (size(p) - m_ < 0) then
      allocate (Pd(1))
      Pd = Zero
      return
    end if

    degree = size(p) - 1
    order = size(p) - m_         ! Number of term of resulting polynomial
    allocate (Pd(order)); Pd = p(:order)

    morder: do k = 0, m_ - 1
      n = degree + 1 - k
      do i = 1, order
        Pd(i) = (n - i) * Pd(i)
      end do
    end do morder
  end function polyder

  ! !> polyint Computes m-esima antiderivative
  ! !!
  ! !! Examples:
  ! !!
  ! function polyint(p, order=m) result(pI)
  !   implicit none
  !   real(dp) :: y !<
  !   real(dp), dimension(:), intent(IN) :: p !<
  !   real(dp), intent(IN) :: xL !<
  !   real(dp), intent(IN) :: xU !<
  !   real(dp), dimension(size(p) + m), intent(IN) :: pI
  !   pI = Zero

  !   IF (m=0) return

  ! end function polyint
end module polynomial