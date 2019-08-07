!> @file   splines.F90
!! @author Juan Fiol <juanfiol@gmail.com>
!! @date   Mon Jul 20 01:53:48 2009
!!
!! @brief implements several functions for simple use of cubic splines.
!! Several sources were reused.
!!
!! Most (all) of the routines were modified from the original.
!! See below for authors of original versions

!> @package splines
!! @brief implementation of interpolation using splines
!!
!! @note
!! It is recommended to use it through the standard interfaces
!! csplrep and splev
!!
!! Examples of use
!! @code
!! USE splines, only: cspl_rep, csplrep, splev
!! type(cspl_rep) :: tck
!! integer, parameter :: N = 500
!! real(8), dimension(N) :: x, y, xnew, ynew
!! ! fill x,y and xnew ...
!! call csplrep(x, y, Zero, Zero, tck)
!! ynew= csplev(xnew, tck)
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

  ! integer, parameter :: NMAX = 800

  !> csplev Performs a spline interpolation in a point or in a table
  !! \see interp_spl and interp_spl_tab
  interface csplev
    module procedure interp_spl
    module procedure interp_spl_tab
  end interface csplev

  private
  public :: cspl_rep, csplev, csplrep, spl_clean_rep, splint, splint_square, spleps
  public :: spline_coeff
contains

  !> cubic spline interpolation between tabulated data
  !! After calling this function the result may be used to evaluate the function as:
  !! `ynew= csplev(xnew, tck)`
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
  function interp_spl(xc, tck) result(y)
    real(dp), intent(IN) :: xc !< value where evaluate the interpolation
    type(cspl_rep), intent(IN) :: tck !<  spline coefficients
    real(dp) :: y                     !< the interpolated value
    integer :: ix
    ix = searchsorted(tck%x, xc)
    y = polyval(tck%S(:, ix), xc)
  end function interp_spl

  !> \copybrief interp_spl
  !! Works over an array of x values and returns an array of interpolated values
  function interp_spl_tab(xnew, tck) result(y)
    type(cspl_rep), intent(IN) :: tck           !< spline coefficients
    real(dp), dimension(:), intent(IN) :: xnew !<  array of x values
    real(dp), dimension(size(xnew)) :: y
    real(dp) :: xc
    integer :: ix, i1

    do concurrent(i1=1:size(xnew))
      xc = xnew(i1)
      ix = searchsorted(tck%x, xc)
      ! y(i1) = polyval(tck%S(:, ix), xc)
      y(i1) = tck%S(4, ix) + (tck%S(3, ix) + (tck%S(2, ix) + tck%S(1, ix) * xc) * xc) * xc
    end do
  end function interp_spl_tab

  !> Coefficients for cubic spline interpolation between tabulated data
  !!
  !! The interpolating polynomial in the i-th interval, from
  !! x(i) to x(i+1), is
  !!        P_i(X) = S(1,i)+x*(S(2,i)+x*(S(3,i)+x*S(4,i)))
  !!  REF.: \ref M82 M.J. Maron, 'Numerical Analysis: A Practical Approach', Macmillan Publ. Co., New York 1982.
  function spline_coeff(X, Y, s1, sN, nn) result(S)
    integer, intent(IN) :: nn                !< dimension of x, y
    real(dp), dimension(nn), intent(IN) :: X !< grid points
    real(dp), dimension(nn), intent(IN) :: Y !< values of the function at grid points x
    real(dp), dimension(4, nn) :: S !< Coefficients, starting with the highest degree
    real(dp), intent(IN) :: s1                   !< second derivative at x(1)
    real(dp), intent(IN) :: sN                   !< second derivative at x(N)
    real(dp) :: R, SI1, SI, H, HI
    integer :: i, k
    integer :: n1, n2

    IF (nn < 4) return
    n1 = nn - 1
    n2 = nn - 2

    ! Check that the points of the grid are all different and ascending. Else abort
    associate (A=>S(4, :), B=>S(3, :), C=>S(2, :), D=>S(1, :))
      A = X(2:) - X(:n1)
      IF (any(A(:n1) < Small)) call print_msg('Points not in strict ascending order',&
        & sub='spline_coeff', errcode=1)

      D(:n1) = (Y(2:) - Y(:nn - 1)) / A(:n1)

      do i = 1, n2   ! Symmetric Coefficient Matrix (augmented).
        B(i) = 2._dp * (A(i) + A(i + 1))
        k = nn - i
        D(k) = 6._dp * (D(k) - D(k - 1))
      enddo
      D(2) = D(2) - A(1) * s1
      D(n1) = D(n1) - A(n1) * sN

      do i = 2, n2   ! Gauss solution of the tridiagonal SYSTEM.
        R = A(i) / B(i - 1)
        B(i) = B(i) - R * A(i)
        D(i + 1) = D(i + 1) - R * D(i)
      enddo
      D(n1) = D(n1) / B(n2)  ! The Sigma Coefficients are stored in array D.

      do i = 2, n2
        k = nn - i
        D(k) = (D(k) - A(k) * D(k + 1)) / B(k - 1)
      enddo
      D(nn) = sN

      !   Spline_coeff Coefficients.
      SI1 = s1
      do i = 1, n1
        SI = SI1
        SI1 = D(i + 1)
        H = A(i)
        HI = 1._dp / H
        A(i) = (HI / 6._dp) * (SI * X(i + 1)**3 - SI1 * X(i)**3)       &
          &      + HI * (Y(i) * X(i + 1) - Y(i + 1) * X(i))            &
          &      + (H / 6._dp) * (SI1 * X(i) - SI * X(i + 1))
        B(i) = (HI / 2._dp) * (SI1 * X(i)**2 - SI * X(i + 1)**2)       &
          &      + HI * (Y(i + 1) - Y(i)) + (H / 6._dp) * (SI - SI1)
        C(i) = (HI / 2._dp) * (SI * X(i + 1) - SI1 * X(i))
        D(i) = (HI / 6._dp) * (SI1 - SI)
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
      RC = x(i)
      YI = polyval(S(:, i - 1), RC)
      if (abs(y(i)) > Eps) then
        Err(i) = 1.0_dp - YI / y(i) ! Relative error
      else
        Err(i) = YI - y(i)      ! Absolute error
      endif
    enddo
  end subroutine spleps

  !> Integral of a cubic spline function.
  !!
  function splint(xL, xU, tck) result(suma)
    real(dp), intent(IN) :: xL  !<   Lower limit in the integral.
    real(dp), intent(IN) :: xU  !<   Upper limit in the integral.
    type(cspl_rep), intent(IN) :: tck !< Interpolating object
    real(dp) :: suma                 !<  Value of integral

    real(dp) :: xll, xuu, x1, x2, sign
    integer :: iL, iU, i

    if (xU > xL) then            ! Set integration limits in increasing order.
      XLL = xL; XUU = xU; sign = 1._dp
    else
      XLL = xU; XUU = xL; sign = -1._dp
    endif

    x1 = tck%x(1); x2 = tck%x(size(tck%x))
    IF (XLL < x1) XLL = x1 + Small ! Check integral limits.
    IF (XUU > x2) XUU = x2 - Small

    ! Find involved intervals.
    iL = searchsorted(tck%x, XLL)
    iU = searchsorted(tck%x, XUU)

    i = iL                    ! Just to use the same expression always
    suma = Zero

    ! First interval
    x1 = XLL; x2 = min(tck%x(i + 1), XUU)
    suma = int_single()
    if (iL == iU) then ! Both limits belong to the same interval
      suma = sign * suma
      return
    end if

    ! Add intermediate intervals
    do i = iL + 1, iU - 1
      x1 = tck%x(i)
      x2 = tck%x(i + 1)
      suma = suma + int_single()
    end do

    ! Last interval
    x1 = x2
    x2 = XUU
    suma = suma + int_single()

    ! Consider the limits of integration
    suma = sign * suma
  contains
    function int_single() result(y)
      real(dp) :: y
      associate (A=>tck%S(4, :), B=>tck%S(3, :), C=>tck%S(2, :), D=>tck%S(1, :))
        y = x2 * (A(i) + x2 * (B(i) / 2 + x2 * (C(i) / 3 + x2 * D(i) / 4)))&
          &  - x1 * (A(i) + x1 * (B(i) / 2 + x1 * (C(i) / 3 + x1 * D(i) / 4)))
      end associate
    end function int_single

  end function splint

  !> Integral of the square of a function expressed as a cubic spline.
  function splint_square(xL, xU, tck) result(suma)
    real(dp), intent(IN) :: xL  !<   Lower limit in the integral.
    real(dp), intent(IN) :: xU  !<   Upper limit in the integral.
    type(cspl_rep), intent(IN) :: tck !< Interpolating object
    real(dp) :: suma                 !<  Value of integral

    real(dp) :: xll, xuu, x1, x2, sign
    integer :: iL, iU, i

    if (xU > xL) then            ! Set integration limits in increasing order.
      XLL = xL; XUU = xU; sign = 1._dp
    else
      XLL = xU; XUU = xL; sign = -1._dp
    endif

    x1 = tck%x(1); x2 = tck%x(size(tck%x))
    IF (XLL < x1) XLL = x1 + Small ! Check integral limits.
    IF (XUU > x2) XUU = x2 - Small

    ! Find involved intervals.
    iL = searchsorted(tck%x, XLL)
    iU = searchsorted(tck%x, XUU)

    i = iL                    ! Just to use the same expression always
    suma = Zero

    ! First interval
    x1 = XLL; x2 = min(tck%x(i + 1), XUU)
    suma = int_single()
    if (iL == iU) then ! Both limits belong to the same interval
      suma = sign * suma
      return
    end if

    ! Add intermediate intervals
    do i = iL + 1, iU - 1
      x1 = tck%x(i)
      x2 = tck%x(i + 1)
      suma = suma + int_single()
    end do

    ! Last interval
    x1 = x2
    x2 = XUU
    suma = suma + int_single()

    ! Consider the limits of integration
    suma = sign * suma
  contains
    function int_single() result(y)
      real(dp) :: y
      associate (A=>tck%S(4, :), B=>tck%S(3, :), C=>tck%S(2, :), D=>tck%S(1, :))
        y = x2 * (A(i) * (A(i) + x2 * B(i))                           &
          &      + x2**2 * ((2 * A(i) * C(i) + B(i)**2) / 3._dp       &
          &      + x2 * ((B(i) * C(i) + A(i) * D(i)) / 2._dp          &
          &      + x2 * ((2 * B(i) * D(i) + C(i)**2) / 5._dp          &
          &      + x2 * D(i) * (C(i) / 3._dp + x2 * D(i) / 7._dp))))) &
          &      - x1 * (A(i) * (A(i) + x1 * B(i))                    &
          &      + x1**2 * ((2 * A(i) * C(i) + B(i)**2) / 3._dp       &
          &      + x1 * ((B(i) * C(i) + A(i) * D(i)) / 2._dp          &
          &      + x1 * ((2 * B(i) * D(i) + C(i)**2) / 5._dp          &
          &      + x1 * D(i) * (C(i) / 3._dp + x1 * D(i) / 7._dp)))))
      end associate
      y = max(y, Zero)
    end function int_single

  end function splint_square

end module csplines

