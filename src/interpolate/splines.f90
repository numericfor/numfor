!> @file   splines.F90
!! @author Juan Fiol <juanfiol@gmail.com>
!! @date   Mon Jul 20 01:53:48 2009
!!
!! @brief implements several functions for simple use of splines. Several sources were
!! reused.
!!
!! Most (all) of the routines were modified from the original. See below for authors of
!! original versions

!> @package splines
!! @brief implementation of interpolation using splines
!!
!! @note
!! It is recommended to use it through the standard interfaces
!! splrep and splev
!!
!! Examples of use
!! @code
!! USE splines, only: spl_rep, splrep, splev
!! type(spl_rep) :: tck
!! integer, parameter :: N = 500
!! real(8), dimension(N) :: x, y, xnew, ynew
!! ! fill x,y and xnew ...
!! call splrep(x, y, Zero, Zero, tck)
!! ynew= splev(xnew, tck)
module splines
  USE basic, only: dp, Zero, Small, print_msg
  USE sorting, only: searchsorted
  implicit none

  !> Type used to keep all information on splines
  type spl_rep
    real(dp), allocatable :: A(:)
    real(dp), allocatable :: B(:)
    real(dp), allocatable :: C(:)
    real(dp), allocatable :: D(:)
    real(dp), allocatable :: x(:)
  end type spl_rep

  integer, parameter :: NMAX = 800

  !> splev Performs a spline interpolation in a point or in a table
  !! \see interp_spl and interp_spl_tab
  interface splev
    module procedure interp_spl
    module procedure interp_spl_tab
  end interface

  private
  public :: spl_rep, splev, splrep, spl_clean_rep, splint, splint_square, spleps
  public :: spline_coeff, INTEG2
contains

  !> cubic spline interpolation between tabulated data
  !! After calling this function the result may be used to evaluate the function as:
  !! `ynew= splev(xnew, tck)`
  !!
  !! REF.: \ref M82 M.J. Maron, 'Numerical Analysis: A Practical Approach', Macmillan Publ. Co., New York 1982.
  subroutine splrep(x, y, s1, sn, r)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !< independent grid points
    real(dp), dimension(:), intent(IN) :: y !< corresponding function values
    real(dp), intent(IN) :: s1              !< second derivative at x(1)
    real(dp), intent(IN) :: sn              !< second derivative at x(n) (Natural spline: s1=sn=0)
    type(spl_rep), intent(OUT) :: r         !< Coefficients stored in type(spl_rep)
    integer :: ierr
    integer :: Nd
    Nd = size(x)
    IF (allocated(r%x) .AND. size(r%x) /= Nd) deallocate (r%x, r%A, r%B, r%C, r%D)
    if (.NOT. allocated(r%x)) then
      allocate (r%x(Nd), r%A(Nd), r%B(Nd), r%C(Nd), r%D(Nd), STAT=ierr)
      IF (ierr /= 0) call print_msg('Allocation error in splrep', errcode=ierr)
    end if
    r%x = x
    ierr = spline_coeff(x, y, r%A, r%B, r%C, r%D, s1, sn, Nd)
    IF (ierr /= 0) call print_msg('Calculation error in splrep', 'splrep', ierr)
  end subroutine splrep

  subroutine spl_clean_rep(r)
    implicit none
    type(spl_rep), intent(INOUT) :: r

    integer :: ierr
    if (allocated(r%x)) then
      deallocate (r%x, r%A, r%B, r%C, r%D, STAT=ierr)
      IF (ierr /= 0) call print_msg('Deallocation error in splrep')
    end if

  end subroutine spl_clean_rep

  !> Interpolates a function using previously calculated representation of splines
  !!
  !! @note Before calling this function, must be called `splrep()`
  function interp_spl(xc, tck) result(y)
    real(dp), intent(IN) :: xc !< value where evaluate the interpolation
    type(spl_rep), intent(IN) :: tck !<  spline coefficients
    real(dp) :: y                    !< the interpolated value
    integer :: ix
    ix = searchsorted(tck%x, xc)
    y = tck%A(ix) + (tck%B(ix) + (tck%C(ix) + tck%D(ix) * xc) * xc) * xc
  end function interp_spl

  !> \copybrief interp_spl
  !! Works over an array of x values and returns an array of interpolated values
  function interp_spl_tab(xnew, tck) result(y)
    type(spl_rep), intent(IN) :: tck           !< spline coefficients
    real(dp), dimension(:), intent(IN) :: xnew !<  array of x values
    real(dp), dimension(size(xnew)) :: y
    real(dp) :: xc
    integer :: ix, i1

    ix = 1
    do i1 = 1, size(xnew)
      xc = xnew(i1)
      ix = searchsorted(tck%x, xc)
      y(i1) = tck%A(ix) + (tck%B(ix) + (tck%C(ix) + tck%D(ix) * xc) * xc) * xc
    end do
  end function interp_spl_tab

  !> cubic spline interpolation between tabulated data
  !! @param X X(I) (I=1, ...,N) grid points.
  !!                (The X values must be in INCREASING ORDER).
  !! @param Y Y(I) (I=1, ...,N) corresponding function values.
  !! @param S1 second derivatives at X(1)
  !! @param SN second derivatives at X(N)
  !!        (THE NATURAL SPLINE CORRESPONDS TO TAKING S1=SN=0).
  !! @param nn Number of grid points.
  !!
  !! The interpolating polynomial in the i-th interval, from
  !! x(i) to x(i+1), is
  !!        P_i(X) = A(i)+x*(B(i)+x*(C(i)+x*D(i)))
  !! @param[out] A  Spline Coefficients (constant).
  !! @param[out] B Spline Coefficients
  !! @param[out] C Spline Coefficients
  !! @param[out] D Spline Coefficients
  !!
  !!  REF.: \ref M82 M.J. Maron, 'Numerical Analysis: A Practical Approach', Macmillan Publ. Co., New York 1982.
  function spline_coeff(X, Y, A, B, C, D, S1, SN, nn) result(spl)
    integer, intent(IN) :: nn
    real(dp), dimension(nn), intent(IN) :: X
    real(dp), dimension(nn), intent(IN) :: Y
    real(dp), dimension(nn), intent(OUT) :: A
    real(dp), dimension(nn), intent(OUT) :: B
    real(dp), dimension(nn), intent(OUT) :: C
    real(dp), dimension(nn), intent(OUT) :: D
    real(dp), intent(IN) :: S1
    real(dp), intent(IN) :: SN
    real(dp) :: R, SI1, SI, H, HI
    integer :: i, k
    integer :: n1, n2
    integer :: spl

    spl = 1
    IF (nn < 4) return
    n1 = nn - 1
    n2 = nn - 2

    ! Check that the points of the grid are all different and ascending
    A = X(2:) - X(:n1)

    IF (any(A(:n1) < Small)) call print_msg('Points not in strict ascending order',&
        & sub='spline_coeff', errcode=1)

    D(:n1) = (Y(2:) - Y(:nn - 1)) / A(:n1)

    do i = 1, n2   ! Symmetric Coefficient Matrix (augmented).
      B(i) = 2._dp * (A(i) + A(i + 1))
      k = nn - i
      D(k) = 6._dp * (D(k) - D(k - 1))
    enddo
    D(2) = D(2) - A(1) * S1
    D(n1) = D(n1) - A(n1) * SN

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
    D(nn) = SN

    !   Spline_coeff Coefficients.
    SI1 = S1
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
    spl = 0
  end function spline_coeff

  !> Estimates the error produced by using a spline approximation
  !!
  !! @param x
  !! @param y value of function at grid points
  !! @param N dimension
  !! @param[out] Err Vector with errors estimates
  !!
  !! @details This subroutine estimates the error introduced by natural cubic spline
  !! interpolation in a table <code> x(i),y(i) (i=1,...,n)</code>.  the interpolation
  !! error in the vicinity of @a x(k) is approximated by the difference between y(k) and the
  !! value obtained from the spline that interpolates the table with the k-th point
  !! removed. err is the largest relative error along the table.
  !!
  !! @note that some tests show that the error is about one order of magnitude better than
  !! estimated by this routine
  !!
  !! @returns Err An array with the relative errors
  subroutine spleps(x, y, Err)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !< grid points
    real(dp), dimension(size(x)), intent(IN) :: y
    real(dp), dimension(size(x)), intent(OUT):: Err
    real(dp), parameter :: Epsilon = 1.e-4_dp
    integer :: i, n1, ierr
    real(dp) :: YI, RC
    real(dp), dimension(size(x) - 1) :: R, F, A, B, C, D

    Err = Zero
    n1 = size(x) - 1
    do i = 2, n1                   ! loop over skipped x-values
      R(1:i - 1) = x(1:i - 1); R(i:) = x(i + 1:)
      F(1:i - 1) = y(1:i - 1); F(i:) = y(i + 1:)
      ierr = spline_coeff(R, F, A, B, C, D, Zero, Zero, n1)
      RC = x(i)
      YI = A(i - 1) + RC * (B(i - 1) + RC * (C(i - 1) + RC * D(i - 1))) ! splev
      if (abs(y(i)) > Epsilon) then
        Err(i) = 1.0_dp - YI / y(i)
      else
        Err(i) = YI - y(i)
      endif
    enddo
  end subroutine spleps

  !> Integral of a cubic spline function.
  !! @param xL  Lower limit in the integral.
  !! @param xU  Upper limit in the integral.
  !! @param tck spline Coefficients
  !! @remarks (slightly modified) From RADIAL package
  !! @return error code (0 in success)
  function splint(xL, xU, tck) result(suma)
    real(dp), intent(IN) :: xL
    real(dp), intent(IN) :: xU
    type(spl_rep), intent(IN) :: tck
    real(dp) :: suma
    real(dp) :: xll, xuu, x1, x2, sumaP, sign
    integer :: N, iL, iU, i

    if (xU > xL) then            ! Set integration limits in increasing order.
      XLL = xL; XUU = xU; sign = 1._dp
    else
      XLL = xU; XUU = xL; sign = -1._dp
    endif

    N = size(tck%x)

    x1 = tck%x(1); x2 = tck%x(N)
    IF (XLL < x1) XLL = x1 + Small ! Check integral limits.
    IF (XUU > x2) XUU = x2 - Small

    ! Find involved intervals.
    iL = searchsorted(tck%x, XLL)
    iU = searchsorted(tck%x, XUU)

    suma = Zero
    if (iL == iU) then    ! Only a single interval involved.
      x1 = XLL
      x2 = XUU
      suma = x2 * (tck%A(iL) + x2 * ((tck%B(iL) / 2) + x2 * ((tck%C(iL) / 3) + x2 * tck%D(iL) / 4))) - &
        & x1 * (tck%A(iL) + x1 * ((tck%B(iL) / 2) + x1 * ((tck%C(iL) / 3) + x1 * tck%D(iL) / 4)))
    else      ! Contributions from different intervals.
      x1 = XLL
      x2 = tck%x(iL + 1)
      suma = x2 * (tck%A(iL) + x2 * ((tck%B(iL) / 2) + x2 * ((tck%C(iL) / 3) + x2 * tck%D(iL) / 4)))&
        & - x1 * (tck%A(iL) + x1 * ((tck%B(iL) / 2) + x1 * ((tck%C(iL) / 3) + x1 * tck%D(iL) / 4)))
      iL = iL + 1
      do i = iL, iU
        x1 = tck%x(i)
        x2 = tck%x(i + 1)
        if (i == iU) x2 = XUU
        sumaP = x2 * (tck%A(i) + x2 * ((tck%B(i) / 2) + x2 * ((tck%C(i) / 3) + x2 * tck%D(i) / 4))) - &
          &  x1 * (tck%A(i) + x1 * ((tck%B(i) / 2) + x1 * ((tck%C(i) / 3) + x1 * tck%D(i) / 4)))
        suma = suma + sumaP
      enddo
    endif
    suma = sign * suma
  end function splint

  !! @param XL  Lower limit in the integral.
  !! @param XU  Upper limit in the integral.
  !! @param[out] SUM value of integral
  !! @param N Number of grid points.
  !! @remarks (slightly modified) From RADIAL package
  !! @return error code (0 in success)
  !! @remarks We should convert it to work with splrep and splev

  !> Integral of the square of a cubic spline function.
  !! @param xL  Lower limit in the integral.
  !! @param xU  Upper limit in the integral.
  !! @param tck spline Coefficients
  !! @return
  function splint_square(xL, xU, tck) result(suma)
    !     INTEGRAL OF A SQUARED CUBIC SPLINE FUNCTION.
    !     INPUT:
    !     X(I) (I=1, ...,N) ........ GRID POINTS.
    !                    (THE X VALUES MUST BE IN INCREASING ORDER).
    !     A(I),B(I),C(I),D(I) ...... SPLINE COEFFICIENTS.
    !     N ........................ NUMBER OF GRID POINTS.
    !     XL ....................... LOWER LIMIT IN THE INTEGRAL.
    !     XU ....................... UPPER LIMIT IN THE INTEGRAL.
    !     OUTPUT:
    !     SUM ...................... VALUE OF THE INTEGRAL.
    !     -----------------------------------------------------------------
    real(dp), intent(IN) :: xL
    real(dp), intent(IN) :: xU
    type(spl_rep), intent(IN) :: tck
    real(dp) :: suma
    real(dp) :: SUMP
    real(dp) :: x1, x2, xll, xuu, sign
    integer :: i, IL, IU, N

    if (xU > xL) then            ! Set integration limits in increasing order.
      XLL = xL; XUU = xU; sign = 1._dp
    else
      XLL = xU; XUU = xL; sign = -1._dp
    endif

    N = size(tck%x)

    x1 = tck%x(1); x2 = tck%x(N)
    IF (XLL < x1) XLL = x1 + Small ! Check integral limits.
    IF (XUU > x2) XUU = x2 - Small

    ! Find involved intervals.
    iL = searchsorted(tck%x, XLL)
    iU = searchsorted(tck%x, XUU)

    suma = Zero
    if (IL == IU) then     ! Only a single interval involved.
      x1 = xLL
      x2 = xUU
      suma = x2 * (tck%A(IL) * (tck%A(IL) + x2 * tck%B(IL))                            &
        &     + x2 * x2 * (((2 * tck%A(IL) * tck%C(IL) + tck%B(IL)**2) / 3._dp)          &
        &     + x2 * (((tck%B(IL) * tck%C(IL) + tck%A(IL) * tck%D(IL)) / 2._dp)        &
        &     + x2 * (((2 * tck%B(IL) * tck%D(IL) + tck%C(IL)**2) / 5._dp)             &
        &     + x2 * tck%D(IL) * ((tck%C(IL) / 3._dp) + x2 * tck%D(IL) / 7._dp)))))      &
        & - x1 * (tck%A(IL) * (tck%A(IL) + x1 * tck%B(IL))                            &
        &     + x1 * x1 * (((2 * tck%A(IL) * tck%C(IL) + tck%B(IL)**2) / 3._dp)          &
        &     + x1 * (((tck%B(IL) * tck%C(IL) + tck%A(IL) * tck%D(IL)) / 2._dp)        &
        &     + x1 * (((2 * tck%B(IL) * tck%D(IL) + tck%C(IL)**2) / 5._dp)             &
        &     + x1 * tck%D(IL) * ((tck%C(IL) / 3._dp) + x1 * tck%D(IL) / 7._dp)))))
    else             ! Contributions from different intervals.
      x1 = xLL
      x2 = tck%x(IL + 1)
      suma = x2 * (tck%A(IL) * (tck%A(IL) + x2 * tck%B(IL))                          &
        &     + x2 * x2 * (((2 * tck%A(IL) * tck%C(IL) + tck%B(IL)**2) / 3._dp)          &
        &     + x2 * (((tck%B(IL) * tck%C(IL) + tck%A(IL) * tck%D(IL)) / 2._dp)        &
        &     + x2 * (((2 * tck%B(IL) * tck%D(IL) + tck%C(IL)**2) / 5._dp)             &
        &     + x2 * tck%D(IL) * ((tck%C(IL) / 3._dp) + x2 * tck%D(IL) / 7._dp)))))      &
        &  - x1 * (tck%A(IL) * (tck%A(IL) + x1 * tck%B(IL))                          &
        &     + x1 * x1 * (((2 * tck%A(IL) * tck%C(IL) + tck%B(IL)**2) / 3._dp)          &
        &     + x1 * (((tck%B(IL) * tck%C(IL) + tck%A(IL) * tck%D(IL)) / 2._dp)        &
        &     + x1 * (((2 * tck%B(IL) * tck%D(IL) + tck%C(IL)**2) / 5._dp)             &
        &     + x1 * tck%D(IL) * ((tck%C(IL) / 3._dp) + x1 * tck%D(IL) / 7._dp)))))
      IL = IL + 1
      do i = IL, IU
        x1 = tck%x(i)
        x2 = tck%x(i + 1)
        if (i == IU) x2 = xUU
        SUMP = x2 * (tck%A(i) * (tck%A(i) + x2 * tck%B(i))                          &
          &      + x2 * x2 * (((2 * tck%A(i) * tck%C(i) + tck%B(i)**2) / 3._dp)         &
          &      + x2 * (((tck%B(i) * tck%C(i) + tck%A(i) * tck%D(i)) / 2._dp)        &
          &      + x2 * (((2 * tck%B(i) * tck%D(i) + tck%C(i)**2) / 5._dp)            &
          &      + x2 * tck%D(i) * ((tck%C(i) / 3._dp) + x2 * tck%D(i) / 7._dp)))))     &
          &  - x1 * (tck%A(i) * (tck%A(i) + x1 * tck%B(i))                          &
          &      + x1 * x1 * (((2 * tck%A(i) * tck%C(i) + tck%B(i)**2) / 3._dp)         &
          &      + x1 * (((tck%B(i) * tck%C(i) + tck%A(i) * tck%D(i)) / 2._dp)        &
          &      + x1 * (((2 * tck%B(i) * tck%D(i) + tck%C(i)**2) / 5._dp)            &
          &      + x1 * tck%D(i) * ((tck%C(i) / 3._dp) + x1 * tck%D(i) / 7._dp)))))
        if (SUMP < Zero) SUMP = Zero ! if negative comes from numerical error
        suma = suma + SUMP
      enddo
    endif
    suma = SIGN * suma
  end function splint_square

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Integral of a squared cubic spline function.
  !! @param[in] X grid points (must be in increasing order). Dimension N
  !! @param[in] A Spline Coefficients
  !! @param[in] B Spline Coefficients
  !! @param[in] C Spline Coefficients
  !! @param[in] D Spline Coefficients
  !! @param[in] XL  Lower limit in the integral.
  !! @param[in] XU  Upper limit in the integral.
  !! @param[out] SUM value of integral
  !! @param N Number of grid points.
  !! @remarks (slightly modified) From RADIAL package
  !! @remarks It is still needed by the subroutine Sch_Bound for Normalization
  !! @return error code (0 in success)
  integer function INTEG2(X, A, B, C, D, XL, XU, SUM, N)
    integer :: N
    real(dp), dimension(N) :: X, A, B, C, D
    real(dp) :: XL, XU, SUM, SUMP
    real(dp) :: X1, X2, XLL, XUU, SIGN
    integer :: i, IL, IU

    integ2 = 0
    if (XU > XL) then    ! Set integration limits in increasing order.
      XLL = XL; XUU = XU; SIGN = 1._dp
    else
      XLL = XU; XUU = XL; SIGN = -1._dp
    endif
    ! Check integral limits.
    if (XLL < X(1) .or. XUU > X(N)) then
      integ2 = -1
      return
    end if

    ! Find involved intervals.
    SUM = Zero
    IL = searchsorted(X, XLL)
    IU = searchsorted(X, XUU)

    if (IL == IU) then     ! Only a single interval involved.
      X1 = XLL
      X2 = XUU
      SUM = X2 * (A(IL) * (A(IL) + X2 * B(IL))                          &
        &   + X2 * X2 * (((2 * A(IL) * C(IL) + B(IL)**2) / 3._dp)          &
        &   + X2 * (((B(IL) * C(IL) + A(IL) * D(IL)) / 2._dp)            &
        &   + X2 * (((2 * B(IL) * D(IL) + C(IL)**2) / 5._dp)             &
        &   + X2 * D(IL) * ((C(IL) / 3._dp) + X2 * D(IL) / 7._dp)))))      &
        &   - X1 * (A(IL) * (A(IL) + X1 * B(IL))                       &
        &   + X1 * X1 * (((2 * A(IL) * C(IL) + B(IL)**2) / 3._dp)          &
        &   + X1 * (((B(IL) * C(IL) + A(IL) * D(IL)) / 2._dp)            &
        &   + X1 * (((2 * B(IL) * D(IL) + C(IL)**2) / 5._dp)             &
        &   + X1 * D(IL) * ((C(IL) / 3._dp) + X1 * D(IL) / 7._dp)))))
    else             ! Contributions from different intervals.
      X1 = XLL
      X2 = X(IL + 1)
      SUM = X2 * (A(IL) * (A(IL) + X2 * B(IL))                          &
        &   + X2 * X2 * (((2 * A(IL) * C(IL) + B(IL)**2) / 3._dp)          &
        &   + X2 * (((B(IL) * C(IL) + A(IL) * D(IL)) / 2._dp)            &
        &   + X2 * (((2 * B(IL) * D(IL) + C(IL)**2) / 5._dp)             &
        &   + X2 * D(IL) * ((C(IL) / 3._dp) + X2 * D(IL) / 7._dp)))))      &
        &   - X1 * (A(IL) * (A(IL) + X1 * B(IL))                       &
        &   + X1 * X1 * (((2 * A(IL) * C(IL) + B(IL)**2) / 3._dp)          &
        &   + X1 * (((B(IL) * C(IL) + A(IL) * D(IL)) / 2._dp)            &
        &   + X1 * (((2 * B(IL) * D(IL) + C(IL)**2) / 5._dp)             &
        &   + X1 * D(IL) * ((C(IL) / 3._dp) + X1 * D(IL) / 7._dp)))))
      IL = IL + 1
      do i = IL, IU
        X1 = X(i)
        X2 = X(i + 1)
        if (i == IU) X2 = XUU
        SUMP = X2 * (A(i) * (A(i) + X2 * B(i))                          &
          &      + X2 * X2 * (((2 * A(i) * C(i) + B(i)**2) / 3._dp)        &
          &      + X2 * (((B(i) * C(i) + A(i) * D(i)) / 2._dp)           &
          &      + X2 * (((2 * B(i) * D(i) + C(i)**2) / 5._dp)           &
          &      + X2 * D(i) * ((C(i) / 3._dp) + X2 * D(i) / 7._dp)))))    &
          &      - X1 * (A(i) * (A(i) + X1 * B(i))                     &
          &      + X1 * X1 * (((2 * A(i) * C(i) + B(i)**2) / 3._dp)        &
          &      + X1 * (((B(i) * C(i) + A(i) * D(i)) / 2._dp)           &
          &      + X1 * (((2 * B(i) * D(i) + C(i)**2) / 5._dp)           &
          &      + X1 * D(i) * ((C(i) / 3._dp) + X1 * D(i) / 7._dp)))))
        if (SUMP < Zero) SUMP = Zero
        SUM = SUM + SUMP
      enddo
    endif
    SUM = SIGN * SUM
  end function INTEG2
end module splines
