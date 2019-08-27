!> @file   splines.F90
!! @author Juan Fiol <juanfiol@gmail.com>
!! @date   Mon Jul 20 01:53:48 2009
!!
!! @brief implements several functions for simple use of cubic splines.
!! Several sources were reused.
!!
!! Most (all) of the routines were modified from the original.
!! See below for authors of original versions

!> csplines implements interpolation using cubic splines
!!
!! @note
!! It is recommended to use it through the standard interfaces
!! csplrep and splev
!!
!! Examples of use
!! @code
!! USE csplines, only: cspl_rep, csplrep, splev
!! type(cspl_rep) :: csp
!! integer, parameter :: N = 500
!! real(8), dimension(N) :: x, y, xnew, ynew
!! ! fill x,y and xnew ...
!! call csplrep(x, y, Zero, Zero, csp)
!! ynew= csplev(xnew, csp)
module csplines
  USE basic, only: dp, Zero, Small, print_msg
  USE sorting, only: searchsorted
  USE polynomial, only: polyder, polyval
  implicit none

  !> Type used to keep all information on splines
  !!
  !! For each interval, the value will be \f$f(t) \approx S(1) t^3 + S(2) t^2 + S(3) t + S(4)\f$
  type cspl_rep
    real(dp), dimension(:, :), allocatable :: S
    real(dp), dimension(:), allocatable :: x
  end type cspl_rep

  !> csplev Performs a spline interpolation in a point or in a table
  !! \see interp_spl and interp_spl_tab
  interface csplev
    module procedure :: interp_spl
    module procedure :: interp_spl_tab
    module procedure :: interp_devspl
    module procedure :: interp_devspl_tab
  end interface csplev

  private
  public :: cspl_rep, csplrep, csplev, spl_clean_rep, splint, splint_square, spleps
  public :: spline_coeff, csplder

contains

  !> cubic spline interpolation between tabulated data
  !! After calling this function the result may be used to evaluate the function as:
  !! `ynew= csplev(xnew, csp)`
  !!
  !! REF.: \ref M82 M.J. Maron, 'Numerical Analysis: A Practical Approach', Macmillan Publ. Co., New York 1982.
  subroutine csplrep(x, y, s1, sn, r)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !< independent grid points
    real(dp), dimension(:), intent(IN) :: y !< corresponding function values
    real(dp), intent(IN) :: s1              !< second derivative at x(1)
    real(dp), intent(IN) :: sn              !< second derivative at x(n) (Natural spline: s1=sn=0)
    type(cspl_rep), intent(OUT) :: r         !< Coefficients stored in type(cspl_rep)
    integer :: ierr
    integer :: Nd
    Nd = size(x)
    IF (allocated(r%x) .AND. size(r%x) /= Nd) deallocate (r%x, r%S)
    if (.NOT. allocated(r%x)) then
      allocate (r%x(Nd), r%S(4, Nd), STAT=ierr)
      IF (ierr /= 0) call print_msg('Allocation error', sub='csplrep', errcode=ierr)
    end if
    r%x = x
    r%S = spline_coeff(x, y, s1, sn, Nd)
  end subroutine csplrep

  !> csplder Computes the derivative of the cubic spline
  !!
  !! Examples:
  !!
  function csplder(csp, m) result(cspd)
    implicit none
    type(cspl_rep), intent(IN) :: csp !<
    integer, intent(IN) :: m !< order of derivation (must be 1 or 2)
    type(cspl_rep) :: cspd    !< Spline representation of derivative
    integer :: Nd
    integer :: i

    if ((m < 1) .or. (m > 2)) call print_msg('Not first or second derivative from cspline')

    Nd = size(csp%x)
    allocate (cspd%x(Nd), cspd%S(4 - m, Nd))
    cspd%x = csp%x
    do i = 1, Nd
      cspd%S(:, i) = polyder(csp%S(:, i), m)
    end do
  end function csplder

  !> Clean-up a spline representation
  subroutine spl_clean_rep(r)
    implicit none
    type(cspl_rep), intent(INOUT) :: r

    integer :: ierr
    if (allocated(r%x)) then
      deallocate (r%x, r%S, STAT=ierr)
      IF (ierr /= 0) call print_msg('Deallocation error', sub='csplrep')
    end if

  end subroutine spl_clean_rep

  !> Interpolates a function using previously calculated representation of splines
  !!
  !! @note Before calling this function, must be called `csplrep()`
  function interp_spl(xc, csp) result(y)
    real(dp), intent(IN) :: xc !< value where evaluate the interpolation
    type(cspl_rep), intent(IN) :: csp !<  spline coefficients
    real(dp) :: y                     !< the interpolated value
    integer :: ix
    ix = searchsorted(csp%x, xc)
    y = polyval(csp%S(:, ix), xc - csp%x(ix))
  end function interp_spl

  !> Interpolates the first derivative of a function
  !!
  !! @note Before calling this function, must be called `csplrep()`
  function interp_devspl(xc, csp, m) result(y)
    real(dp), intent(IN) :: xc !< value where evaluate the interpolation
    type(cspl_rep), intent(IN) :: csp !<  spline coefficients
    integer, intent(IN) :: m                      !< order of derivation
    real(dp) :: y                     !< the interpolated value
    integer :: ix
    ix = searchsorted(csp%x, xc)
    if (m /= 1 .or. m /= 2) call print_msg('Not first or second derivative from cspline')
    y = polyval(polyder(csp%S(:, ix), m), xc - csp%x(ix))
  end function interp_devspl

  function interp_devspl_tab(xnew, csp, m) result(y)
    type(cspl_rep), intent(IN) :: csp           !< spline coefficients
    real(dp), dimension(:), intent(IN) :: xnew !<  array of x values
    real(dp), dimension(size(xnew)) :: y
    integer, intent(IN) :: m    !< order of derivative
    type(cspl_rep) :: cspd
    real(dp) :: xc
    integer :: ix, i1

    if (m /= 1 .and. m /= 2) call print_msg('Not first or second derivative from cspline')
    cspd = csplder(csp, m)      ! Calculate the coefficients of derivative
    do concurrent(i1=1:size(xnew))
      xc = xnew(i1)
      ix = searchsorted(cspd%x, xc)
      y(i1) = polyval(cspd%S(:, ix), xc - csp%x(ix))
    end do
  end function interp_devspl_tab

  !> \copybrief interp_spl
  !! Works over an array of x values and returns an array of interpolated values
  function interp_spl_tab(xnew, csp) result(y)
    type(cspl_rep), intent(IN) :: csp           !< spline coefficients
    real(dp), dimension(:), intent(IN) :: xnew !<  array of x values
    real(dp), dimension(size(xnew)) :: y
    real(dp) :: xc
    integer :: ix, i1

    do concurrent(i1=1:size(xnew))
      xc = xnew(i1)
      ix = searchsorted(csp%x, xc)
      y(i1) = polyval(csp%S(:, ix), xc - csp%x(ix))
    end do
  end function interp_spl_tab

  !> Coefficients for cubic spline interpolation between tabulated data
  !!
  !! The interpolating polynomial in the i-th interval, from
  !! x(i) to x(i+1), is
  !!       \f$ P_i(X) = S(4,i)+h (S(3,i)+h (S(2,i)+h S(1,i))) \f$
  !!       \f$ P_i(X) = S(1,i) h^3 + S(2,i) h^2 S(3,i) h + S(4,i) \f$
  !!       with \f$ h = x-x(i) \f$
  !!  Ref: \ref M82 M.J. Maron, 'Numerical Analysis: A Practical Approach', Macmillan Publ. Co., New York 1982.
  function spline_coeff(x, y, s1, sN, nn) result(S)
    integer, intent(IN) :: nn                !< dimension of x, y
    real(dp), dimension(nn), intent(IN) :: x !< grid points
    real(dp), dimension(nn), intent(IN) :: y !< values of the function at grid points x
    real(dp), dimension(4, nn) :: S !< Coefficients, starting with the highest degree
    real(dp), intent(IN) :: s1      !< second derivative at x(1)
    real(dp), intent(IN) :: sN      !< second derivative at x(N)
    real(dp) :: R, SI1, SI, hi, ih
    integer :: i, k
    integer :: n1, n2

    IF (nn < 4) return
    n1 = nn - 1
    n2 = nn - 2

    associate (A=>S(4, :), B=>S(3, :), C=>S(2, :), D=>S(1, :))

      A = x(2:) - x(:n1)
      ! Check that the points of the grid are all different and ascending. Else abort
      ! Scipy leave that to the user (too expensive?)
      IF (any(A(:n1) < Small)) call print_msg('Points not in strict ascending order',&
        & sub='spline_coeff', errcode=1)

      D(:n1) = (y(2:) - y(:n1)) / A(:n1)

      ! Symmetric Coefficient Matrix (augmented).
      B(:n2) = 2._dp * (A(:n2) + A(2:n1))
      D(2:n1) = 6._dp * (D(2:n1) - D(1:n2))
      D(2) = D(2) - A(1) * s1
      D(n1) = D(n1) - A(n1) * sN

      do i = 2, n2   ! Gauss solution of the tridiagonal system.
        R = A(i) / B(i - 1)
        B(i) = B(i) - R * A(i)
        D(i + 1) = D(i + 1) - R * D(i)
      enddo
      D(n1) = D(n1) / B(n2)  ! Sigma Coefficients

      do k = n2, 2, -1
        D(k) = (D(k) - A(k) * D(k + 1)) / B(k - 1)
      end do
      D(nn) = sN

      ! Now we have:
      ! A = x_{i+1}-x_{i}
      ! D = sigma_{i}
      ! ---- Spline Coefficients.-------
      ! Bi= -(2*si/6+s1/6)*hi-(fi-f1)*ih
      ! Di= (s1/6-si/6)*ih
      ! Ci= 3*si/6
      ! Ai= fi
      SI1 = s1 / 6._dp
      do i = 1, n1
        SI = SI1
        SI1 = D(i + 1) / 6._dp
        hi = A(i)
        ih = 1._dp / hi
        A(i) = y(i)
        B(i) = -(SI1 + 2 * SI) * hi + (y(i + 1) - y(i)) * ih
        C(i) = 3 * SI
        D(i) = ih * (SI1 - SI)
      enddo
    end associate
  end function spline_coeff

  !> Estimates the error produced by using a spline approximation
  !!
  !! @details This subroutine estimates the error introduced by natural cubic spline
  !! interpolation in a table <code> x(i),y(i) (i=1,...,n)</code>.  the interpolation
  !! error in the vicinity of @a x(k) is approximated by the difference between y(k) and the
  !! value obtained from the spline that interpolates the table with the k-th point
  !! removed. err is the largest relative error along the table.
  !!
  !! @note Some tests seems to show that the error is about one order of magnitude better than
  !! estimated by this routine
  !! @note Modified from Salvat et al, Comp. Phys. Comm. (1995)
  !!
  !! @returns Err An array with the relative errors
  subroutine spleps(x, y, Err)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !< grid points
    real(dp), dimension(size(x)), intent(IN) :: y !< value of function at grid points
    real(dp), dimension(size(x)), intent(OUT):: Err !< Vector with error estimates
    real(dp), parameter :: Eps = 1.e-4_dp
    integer :: i, n1
    real(dp) :: YI, RC
    real(dp), dimension(size(x) - 1) :: R, F
    real(dp), dimension(4, size(x) - 1) :: S

    Err = Zero
    n1 = size(x) - 1
    do i = 2, n1                   ! Loop over skipped x-values
      R(1:i - 1) = x(1:i - 1); R(i:) = x(i + 1:)
      F(1:i - 1) = y(1:i - 1); F(i:) = y(i + 1:)

      S = spline_coeff(R, F, Zero, Zero, n1)
      RC = x(i) - x(i - 1)
      YI = polyval(S(:, i - 1), RC)
      if (abs(y(i)) > Eps) then
        Err(i) = 1.0_dp - YI / y(i) ! Relative error
      else
        Err(i) = YI - y(i)      ! Absolute error
      endif
    enddo
  end subroutine spleps

  !> Definite integral of a cubic spline function.
  !!
  function splint(xL, xU, csp) result(suma)
    real(dp), intent(IN) :: xL  !<   Lower limit in the integral.
    real(dp), intent(IN) :: xU  !<   Upper limit in the integral.
    type(cspl_rep), intent(IN) :: csp !< Interpolating object
    real(dp) :: suma                 !<  Value of integral

    real(dp) :: xll, xuu, x1, x2, sign
    integer :: iL, iU, i

    if (xU > xL) then            ! Set integration limits in increasing order.
      XLL = xL; XUU = xU; sign = 1._dp
    else
      XLL = xU; XUU = xL; sign = -1._dp
    endif

    x2 = csp%x(size(csp%x))
    IF (XLL < csp%x(1)) XLL = csp%x(1) + Small ! Check and correct integral limits.
    IF (XUU > x2) XUU = x2 - Small

    ! Find involved intervals.
    iL = searchsorted(csp%x, XLL)
    iU = searchsorted(csp%x, XUU)

    suma = Zero

    ! First interval (i = iL)
    i = iL                    ! Just to use always the same expression
    x1 = XLL - csp%x(i)
    x2 = min(csp%x(i + 1), XUU) - csp%x(i)
    suma = int_single(x2) - int_single(x1)
    if (iL == iU) then ! Both limits belong to the same interval
      suma = sign * suma
      return
    end if

    ! Add intermediate intervals
    do i = iL + 1, iU - 1
      x2 = csp%x(i + 1) - csp%x(i)
      suma = suma + int_single(x2)
    end do

    ! Last interval (i = iU)
    x2 = XUU - csp%x(i)
    suma = suma + int_single(x2)

    ! Consider the order of the limits of integration
    suma = sign * suma
  contains
    function int_single(xx) result(y)
      real(dp), intent(IN) :: xx
      real(dp) :: y
      associate (A=>csp%S(4, :), B=>csp%S(3, :), C=>csp%S(2, :), D=>csp%S(1, :))
        y = xx * (A(i) + xx * (B(i) / 2 + xx * (C(i) / 3 + xx * D(i) / 4)))
      end associate
    end function int_single

  end function splint

  !> Integral of the square of a function expressed as a cubic spline.
  function splint_square(xL, xU, csp) result(suma)
    real(dp), intent(IN) :: xL  !<   Lower limit in the integral.
    real(dp), intent(IN) :: xU  !<   Upper limit in the integral.
    type(cspl_rep), intent(IN) :: csp !< Interpolating object
    real(dp) :: suma                 !<  Value of integral

    real(dp) :: xll, xuu, x1, x2, sign
    integer :: iL, iU, i

    if (xU > xL) then            ! Set integration limits in increasing order.
      XLL = xL; XUU = xU; sign = 1._dp
    else
      XLL = xU; XUU = xL; sign = -1._dp
    endif

    x2 = csp%x(size(csp%x))
    IF (XLL < csp%x(1)) XLL = csp%x(1) + Small ! Check and correct integral limits.
    IF (XUU > x2) XUU = x2 - Small

    ! Find involved intervals.
    iL = searchsorted(csp%x, XLL)
    iU = searchsorted(csp%x, XUU)

    suma = Zero

    ! First interval (i = iL)
    i = iL                    ! Just to use always the same expression
    x1 = XLL - csp%x(i)
    x2 = min(csp%x(i + 1), XUU) - csp%x(i)
    suma = int_single(x2) - int_single(x1)
    if (iL == iU) then ! Both limits belong to the same interval
      suma = sign * suma
      return
    end if

    ! Add intermediate intervals
    do i = iL + 1, iU - 1
      x2 = csp%x(i + 1) - csp%x(i)
      suma = suma + int_single(x2)
    end do

    ! Last interval (i = iU)
    x2 = XUU - csp%x(i)
    suma = suma + int_single(x2)

    ! Consider the order of the limits of integration
    suma = sign * suma

  contains
    function int_single(xx) result(y)
      real(dp), intent(IN) :: xx
      real(dp) :: y
      associate (A=>csp%S(4, :), B=>csp%S(3, :), C=>csp%S(2, :), D=>csp%S(1, :))
        y = xx * (A(i) * (A(i) + xx * B(i))                           &
          &      + xx**2 * ((2 * A(i) * C(i) + B(i)**2) / 3._dp       &
          &      + xx * ((B(i) * C(i) + A(i) * D(i)) / 2._dp          &
          &      + xx * ((2 * B(i) * D(i) + C(i)**2) / 5._dp          &
          &      + xx * D(i) * (C(i) / 3._dp + xx * D(i) / 7._dp)))))
      end associate
      y = max(y, Zero)          ! If less than zero is an error
    end function int_single

  end function splint_square

end module csplines

