!> @file csplines.f90
!! @author Juan Fiol <juanfiol@gmail.com>
!! @date "2019-12-19 08:35:06"
!!
!! @brief implements several functions for simple use of cubic splines.
!!
!! @note Several sources were reused.
!! All routines were modified from the original (and some heavily modified).
!! See below for authors of original and earlier versions

!> csplines implements interpolation using cubic splines
!! Description: @ref docinterpolate
module csplines
  USE basic, only: dp, Zero, Small, print_msg
  USE sorting, only: searchsorted
  USE polynomial, only: polyder, polyval, polyint, bisect_pol
  USE fitpack, only: fpcuro
  implicit none

  !> Type used to keep all information on splines
  !!
  !! For each interval, the value will be \f$f(t) \approx S(1) t^3 + S(2) t^2 + S(3) t + S(4)\f$ where \f$t = x - x_{i}\f$ is the distance to the nearest point on the left.
  type CubicSpline
    real(dp), dimension(:, :), allocatable :: S !<  Coefficients of the polynomial  `S(:,i)` valid for each interval i, `x(i) <= x < x(i+1)`
    real(dp), dimension(:), allocatable :: x    !< Limits of the intervals
  contains
    ! Evaluation
    procedure, pass(csp) :: cspl_interp
    procedure, pass(csp) :: cspl_interp_tab
    generic :: evaluate => cspl_interp, cspl_interp_tab
    ! Evaluation of derivative
    procedure, pass(csp) :: cspl_interpdev
    procedure, pass(csp) :: cspl_interpdev_tab
    generic :: derivative => cspl_interpdev, cspl_interpdev_tab
    ! Evaluation of integral
    procedure, pass(csp) :: integrate => csplint

    ! Evaluation of root
    procedure, pass(csp) :: roots => csplroots

  end type CubicSpline

  !> csplev Performs a spline interpolation in a point or in a table
  !! \see cspl_interp and cspl_interp_tab
  !!
  !! Example:
  !! --------
  !! @snippet ex_interp_csplines_fp.f90 Using eval
  !!
  interface csplev
    module procedure :: cspl_interp
    module procedure :: cspl_interp_tab
  end interface csplev

  !> csplev Performs a spline interpolation in a point or in a table
  !! \see cspl_interpdev and cspl_interpdev_tab
  !!
  !! Example:
  !! --------
  !! @snippet ex_interp_csplines_fp.f90 Using derivative
  !!
  interface csplevder
    module procedure :: cspl_interpdev
    module procedure :: cspl_interpdev_tab
  end interface csplevder

  !> CubicSpline is the OO interface to cubic splines interpolation
  !!
  !! Example:
  !! --------
  !! @snippet ex_interp_csplines_oo.f90 Using eval
  !!
  !! And to evaluate the derivative
  !!
  !! @snippet ex_interp_csplines_oo.f90 Using derivative
  !!
  interface CubicSpline
    module procedure :: init
  end interface CubicSpline

  Private
  Public :: CubicSpline, cspl_clean, csplrep, csplev, csplevder, csplint, csplint_square, cspleps, csplroots
  Public :: csplder, csplantider

contains

  function init(x, y, s1, sn) result(csp)
    implicit none
    type(CubicSpline) :: csp         !< Coefficients stored in type(CubicSpline)
    real(dp), dimension(:), intent(IN) :: x !< independent grid points
    real(dp), dimension(:), intent(IN) :: y !< corresponding function values
    real(dp), intent(IN) :: s1              !< second derivative at x(1)
    real(dp), intent(IN) :: sn              !< second derivative at x(n) (Natural spline: s1=sn=0)

    call csplrep(x, y, s1, sn, csp)
  end function init

  !> cubic spline interpolation between tabulated data
  !! After calling this function the result may be used to evaluate the function as:\n
  !! `ynew= csplev(xnew, csp)`\n
  !!
  !! REF.: M.J. Maron, 'Numerical Analysis: A Practical Approach', Macmillan Publ. Co., New York 1982.
  subroutine csplrep(x, y, s1, sn, csp)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !< independent grid points
    real(dp), dimension(:), intent(IN) :: y !< corresponding function values
    real(dp), intent(IN) :: s1              !< second derivative at x(1)
    real(dp), intent(IN) :: sn              !< second derivative at x(n) (Natural spline: s1=sn=0)
    type(CubicSpline), intent(OUT) :: csp   !< Coefficients stored in

    integer :: ierr
    integer :: Nd
    Nd = size(x)
    IF (allocated(csp%x) .AND. size(csp%x) /= Nd) deallocate (csp%x, csp%S)
    if (.NOT. allocated(csp%x)) then
      allocate (csp%x(Nd), csp%S(4, Nd - 1), STAT=ierr)
      IF (ierr /= 0) call print_msg('Allocation error', sub='csplrep', errcode=ierr)
    end if
    csp%x = x
    csp%S = spline_coeff(x, y, s1, sn, Nd)
  end subroutine csplrep

  !> csplder Computes the derivative of the cubic spline
  function csplder(csp, m) result(cspd)
    implicit none
    type(CubicSpline), intent(IN) :: csp !< Interpolating object
    integer, intent(IN) :: m   !< order of derivation (must be 1 or 2)
    type(CubicSpline) :: cspd  !< Spline representation of derivative
    integer :: Nd
    integer :: i

    if ((m < 1) .or. (m > 2)) call print_msg('Not first or second derivative from cspline', &
      &"csplder", errcode=m)

    Nd = size(csp%x)
    allocate (cspd%x(Nd), cspd%S(4 - m, Nd))
    cspd%x = csp%x
    do i = 1, Nd
      cspd%S(:, i) = polyder(csp%S(:, i), m)
    end do
  end function csplder

  !> Computes the antiderivative of the CubicSpline approximation
  !!
  !! The result is a CubicSpline (not exactly, it is a polynomial of order m+4)
  function csplantider(csp, m) result(cspa)
    implicit none
    type(CubicSpline), intent(IN) :: csp !< CubicSpline object holding the spline
    integer, intent(IN) :: m             !< order of integration
    type(CubicSpline) :: cspa            !< spline holding the antiderivative

    integer :: Nd
    integer :: i

    Nd = size(csp%x)
    allocate (cspa%x(Nd), cspa%S(4 + m, Nd))
    cspa%x = csp%x
    do i = 1, Nd
      cspa%S(:, i) = polyint(csp%S(:, i), m)
    end do
  end function csplantider

  !> Clean-up a spline representation
  subroutine cspl_clean(r)
    implicit none
    type(CubicSpline), intent(INOUT) :: r

    integer :: ierr
    if (allocated(r%x)) then
      deallocate (r%x, r%S, STAT=ierr)
      IF (ierr /= 0) call print_msg('Deallocation error', sub='csplrep')
    end if
  end subroutine cspl_clean

  !> Interpolates a function using previously calculated representation of splines
  !!
  !! @note Before calling this function, must be called `csplrep()`
  function cspl_interp(xc, csp) result(y)
    real(dp), intent(IN) :: xc !< value where evaluate the interpolation
    class(CubicSpline), intent(IN) :: csp !<  spline coefficients
    real(dp) :: y                     !< the interpolated value
    integer :: ix
    ix = searchsorted(csp%x, xc)
    y = polyval(csp%S(:, ix), xc - csp%x(ix))
  end function cspl_interp

  !> \copybrief cspl_interp
  !! Works over an array of x values and returns an array of interpolated values
  function cspl_interp_tab(xnew, csp) result(y)
    class(CubicSpline), intent(IN) :: csp           !< spline coefficients
    real(dp), dimension(:), intent(IN) :: xnew !<  array of x values
    real(dp), dimension(size(xnew)) :: y       !< Interpolated values at xnew positions
    real(dp) :: xc
    integer :: ix, i1

    do concurrent(i1=1:size(xnew))
      xc = xnew(i1)
      ix = searchsorted(csp%x, xc)
      y(i1) = polyval(csp%S(:, ix), xc - csp%x(ix))
    end do
  end function cspl_interp_tab

  !> Interpolates the first derivative of a function
  !!
  !! @note Before calling this function, must be called `csplrep()`
  function cspl_interpdev(xc, csp, m) result(y)
    real(dp), intent(IN) :: xc !< value where evaluate the interpolation
    class(CubicSpline), intent(IN) :: csp !<  spline coefficients
    integer, optional, intent(IN) :: m                      !< order of derivation
    integer :: m_                      !< order of derivation
    real(dp) :: y                     !< the interpolated value
    integer :: ix
    m_ = 1; IF (present(m)) m_ = m
    ix = searchsorted(csp%x, xc)
    if (m_ /= 1 .and. m_ /= 2) call print_msg('Not first or second derivative from cspline', &
      & "cspl_interpdev", errcode=m_)
    y = polyval(polyder(csp%S(:, ix), m_), xc - csp%x(ix))
  end function cspl_interpdev

  function cspl_interpdev_tab(xnew, csp, m) result(y)
    class(CubicSpline), intent(IN) :: csp           !< spline coefficients
    real(dp), dimension(:), intent(IN) :: xnew !<  array of x values
    real(dp), dimension(size(xnew)) :: y       !< Interpolated values
    integer, optional, intent(IN) :: m    !< order of derivative

    real(dp) :: xc
    integer :: ix, i1
    integer :: m_
    ! Register which polynomial were already calculated
    logical, dimension(size(csp%S(1, :))) :: tofill
    real(dp), dimension(3, size(csp%x)) :: pol

    m_ = 1; IF (present(m)) m_ = m
    if (m_ /= 1 .and. m_ /= 2) call print_msg('Not first or second derivative from cspline', &
      &"cspl_interpdev_tab", errcode=m_)

    tofill = .True.
    do i1 = 1, size(xnew)
      xc = xnew(i1)
      ix = searchsorted(csp%x, xc)
      IF (ix < 0 .or. ix > size(csp%x)) call print_msg('Error in indice', 'cspl_interpdev_tab', abs(ix))
      if (tofill(ix)) then
        pol(:, ix) = polyder(csp%S(:, ix), m_)      ! Calculate the coefficients of derivative
        tofill(ix) = .False.
      end if
      y(i1) = polyval(pol(:, ix), xc - csp%x(ix))
    end do
  end function cspl_interpdev_tab

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
    real(dp), dimension(4, nn - 1) :: S !< Coefficients, starting with the highest degree
    real(dp), intent(IN) :: s1      !< second derivative at x(1)
    real(dp), intent(IN) :: sN      !< second derivative at x(N)
    real(dp) :: R, SI1, SI, hi, ih
    integer :: i, k
    integer :: n1, n2

    IF (nn < 4) return
    n1 = nn - 1
    n2 = nn - 2

    associate (A=>S(4, :), B=>S(3, :), C=>S(2, :), D=>S(1, :))

      A(:n1) = x(2:) - x(:n1)
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
        B(i) = -(SI1 + 2 * SI) * hi + (y(i + 1) - y(i)) * ih
        C(i) = 3 * SI
        D(i) = ih * (SI1 - SI)
        A(i) = y(i)
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
  subroutine cspleps(x, y, Err)
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
  end subroutine cspleps

  !> Definite integral of a cubic spline function.
  !! @todo Testing when extrapolate is True
  function csplint(xL, xU, csp, extrapolate) result(suma)
    real(dp), intent(IN) :: xL  !<   Lower limit in the integral.
    real(dp), intent(IN) :: xU  !<   Upper limit in the integral.
    class(CubicSpline), intent(IN) :: csp !< Interpolating object
    logical, optional, intent(IN) :: extrapolate !< Flag signaling if we extrapolate outside interval
    real(dp) :: suma                 !<  Value of integral

    real(dp) :: xll, xuu, x1, x2, sign
    integer :: iL, iU, i
    logical :: ext_

    ext_ = (Present(extrapolate) .and. (extrapolate))

    if (xU > xL) then            ! Set integration limits in increasing order.
      XLL = xL; XUU = xU; sign = 1._dp
    else
      XLL = xU; XUU = xL; sign = -1._dp
    endif

    x1 = csp%x(1)     ! Last value
    x2 = csp%x(size(csp%x))     ! Last value

    suma = Zero

    ! Here fix the limits according to argument extrapolate
    if (ext_) then
      if (XLL < x1) then  ! Add interval before the first point
        i = 1
        suma = suma + int_single(x1) - int_single(XLL)
        XLL = x1
      end if
      if (XUU > x2) then        ! Add interval after the last point
        suma = suma + int_single(XUU) - int_single(x2)
        XUU = x2
      end if
    end if

    ! Find involved intervals.
    iL = searchsorted(csp%x, XLL)
    iU = searchsorted(csp%x, XUU)

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
      real(dp), intent(IN) :: xx  ! Position of interval
      real(dp) :: y               ! Integrate a single interval
      associate (A=>csp%S(4, :), B=>csp%S(3, :), C=>csp%S(2, :), D=>csp%S(1, :))
        y = xx * (A(i) + xx * (B(i) / 2 + xx * (C(i) / 3 + xx * D(i) / 4)))
      end associate
    end function int_single

  end function csplint

  !> Integral of the square of a function expressed as a cubic spline.
  !! @todo Testing
  function csplint_square(xL, xU, csp) result(suma)
    real(dp), intent(IN) :: xL  !<   Lower limit in the integral.
    real(dp), intent(IN) :: xU  !<   Upper limit in the integral.
    type(CubicSpline), intent(IN) :: csp !< Interpolating object
    real(dp) :: suma                 !< Value of integral

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
    suma = int_single_a(x2) - int_single_a(x1)
    if (iL == iU) then ! Both limits belong to the same interval
      suma = sign * suma
      return
    end if

    ! Add intermediate intervals
    do i = iL + 1, iU - 1
      x2 = csp%x(i + 1) - csp%x(i)
      suma = suma + int_single_a(x2)
    end do

    ! Last interval (i = iU)
    x2 = XUU - csp%x(i)
    suma = suma + int_single_a(x2)

    ! Consider the order of the limits of integration
    suma = sign * suma

  contains

    ! Integrate a single interval
    function int_single_a(xx) result(y)
      real(dp), intent(IN) :: xx ! Position of the given interval
      real(dp) :: y
      associate (A=>csp%S(4, :), B=>csp%S(3, :), C=>csp%S(2, :), D=>csp%S(1, :))
        y = xx * (A(i) * (A(i) + xx * B(i))                           &
          &      + xx**2 * ((2 * A(i) * C(i) + B(i)**2) / 3._dp       &
          &      + xx * ((B(i) * C(i) + A(i) * D(i)) / 2._dp          &
          &      + xx * ((2 * B(i) * D(i) + C(i)**2) / 5._dp          &
          &      + xx * D(i) * (C(i) / 3._dp + xx * D(i) / 7._dp)))))
      end associate
      y = max(y, Zero)          ! If less than zero is an error
    end function int_single_a

  end function csplint_square

  !> csplroots Computes the roots of the Spline approximation
  function csplroots(csp) result(z)
    implicit none
    class(CubicSpline), intent(IN) :: csp !< Spline approximation to consider
    real(dp), dimension(:), allocatable :: z !< Roots (zeros) of the spline function
    !!
    !! Examples:
    !! --------
    !! ```
    !!   real(dp), dimension(:), allocatable :: zeros
    !!   zeros = csp%roots()
    !!```
    real(dp), dimension(3 * size(csp%x)) :: c
    integer :: i, nc, n, j, k
    real(dp), parameter :: toler = 1.e-5_8

    nc = 0
    c = 0._dp
    do i = 1, size(csp%x) - 1
      call fpcuro(csp%S(1, i), csp%S(2, i), csp%S(3, i), csp%S(4, i), c(nc + 1:nc + 3), n)
      k = 0
      do j = 1, n
        IF (c(nc + j) >= -toler .and. c(nc + j) - toler < csp%x(i + 1) - csp%x(i)) then
          k = k + 1
          c(nc + k) = c(nc + j) + csp%x(i)
        end IF
        nc = nc + k
      end do

    end do
    if (nc > 0) then
      allocate (z(nc))
      z(:nc) = c(:nc)
    end if
  end function csplroots

end module csplines

! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:

