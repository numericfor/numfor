!> fitpack provides a framework for fitting and interpolation using B-Splines
module fitpack

  USE utils, only: dp, str, print_msg

  !> Type used to keep all information on spline fitting
  !!
  !!
  type UnivSpline
    real(dp), dimension(:), allocatable :: t !< knots
    real(dp), dimension(:), allocatable :: c !< coefficients
    integer :: k
    real(dp) :: fp !< weighted sum of squared residuals of the spline approximation returned.
    real(dp), dimension(:), allocatable :: wrk
    integer, dimension(:), allocatable :: iwrk
    integer :: ier
  end type UnivSpline

  Public :: UnivSpline, splrep, splrep_msg

  character, private, parameter :: nl = new_line('a')

contains

  !> splrep_msg
  !!
  !! Examples:
  !!
  subroutine splrep_msg(flag)
    implicit none
    integer, intent(IN) :: flag !< Value of i
    character(len=:), allocatable :: msg
    select case (flag)
    case (0)
      msg = "Normal return. The spline returned has a residual sum of squares fp such that abs(fp - s)/s <= 1.e-3"
    case (-1)
      msg = "Normal return. The spline interpolates the data (fp=0)."
    case (-2)
      msg = "Normal return. The spline is the weighted least-squares polynomial of degree k."//nl//" &
        &    fp gives the upper bound fp0 for the smoothing factor s."
    case (1)
      msg = "Error. The required storage space exceeds the available storage space."//nl//"&
        & Probable causes: size(x) too small or smoothing parameter s is too small (fp > s) ."
    case (2)
      msg = "Error. A theoretically impossible result when finding a smoothing spline with fp = s."//nl//&
        &"Probable cause: s too small. abs(fp - s) / s > 0.001"
    case (3)
      msg = "Error. The maximal number of iterations (maxit=20) allowed for finding a smoothing spline &
        &with fp = s has been reached."//nl//"Probable cause: s too small abs(fp - s) / s > 0.001."
    case (10)
      msg = "Error. On entry the input data are controlled on validity the following restrictions must be satisfied."//nl//&
        &"  -1 <= task <= 1,"//nl//"   1 <= k <= 5,"//nl//"   m > k,"//nl//" w(i) > 0, i = 1, 2, ..., m "//nl//&
        &"  xb <= x(1) < x(2) < ... < x(m) <= xe,"//nl// &
        &"  if task == -1: then 2 * k + 2 <= n "//nl//&
        &"  xb < t(k + 2) < t(k + 3) < ... < t(n - k - 1) < xe "//&
        &"  The Schoenberg-Whitney conditions, i.e.there must be a subset of data points xx(j) such that "//nl//&
        &"     t(j) < xx(j) < t(j + k + 1), j = 1, 2, ..., n - k - 1"//nl//&
        &"  if task >= 0:  s >= 0."
    end select

    print *, msg
  end subroutine splrep_msg

  !> splprep Computes the B-spline representation of an N-dimensional curve.
  !!
  !! Given a list of N rank-1 arrays, `x`, which represent a curve in N-dimensional space parametrized by `u`, find a smooth
  !! approximating spline curve g(`u`). Uses the routine parcur from (slightly modified) FITPACK.
  !! Examples:
  !!
  subroutine splprep(x, w, u, ub, ue, k, task, s, t, per, tck, ier)
    ! subroutine splprep(x, y, w, ub, ue, k, task, s, t, per, tck)
    implicit none
    real(dp), dimension(:, :), intent(IN) :: x !<
    real(dp), dimension(:), intent(INOUT) :: u !< An array with parameters. If
    real(dp), optional, dimension(size(x(1))), target, intent(IN) :: w !< Strictly positive rank-1 array of weights the same size as u.

    !!The weights are used in computing the weighted least-squares spline fit. If the errors in the y values have standard-deviation
    !!given by the vector d, then w should be 1/d.
    real(dp), optional, intent(INOUT) :: ub !< Lower limit of the interval to approximate (ub <= x(1))
    real(dp), optional, intent(INOUT) :: ue !< Upper limit of the interval to approximate (ue >= x(size(x))).
    integer, optional, intent(IN) :: k !< The degree of the spline fit. It is recommended to use cubic splines.
    !! Even values of k should be avoided especially with small s values. 1 <= k <= 5
    integer, optional, intent(IN) :: task !< {1, 0, -1},
    !! If task==0 find t and c for a given smoothing factor, s.
    !!
    !! If task==1 find t and c for another value of the smoothing factor, s.
    !! There must have been a previous call with task=0 or task=1 for the same set of data (it will be stored an used internally)
    !!
    !!If task==-1 find the weighted least square spline for a given set of knots t.
    !!    These should be interior knots as knots on the ends will be added automatically
    real(dp), optional, intent(IN) :: s !< A smoothing condition. The amount of smoothness is determined by satisfying the conditions:
    !!
    !! `sum((w * (y - g))**2) <= s` where `g(x)` is the smoothed interpolation of `(x,y)`.
    !!
    !!The user can use `s` to control the tradeoff between closeness and smoothness of fit. Larger `s` means more smoothing while
    !! smaller values of `s` indicate less smoothing.  Recommended values of s depend on the weights, w.
    !!
    !!If the weights represent the inverse of the standard-deviation of y, then a good s value should be found in the range
    !!(m-sqrt(2*m),m+sqrt(2*m)) where m is the number of datapoints in x, y, and w. default : \f$ s= m-\sqrt(2 m)\f$ if weights are
    !!supplied. `s = 0.0` (interpolating) if no weights are supplied.

    real(dp), optional, dimension(:), intent(IN) :: t !< Input knots (interior knots only). if task = -1. If given then task is automatically set to -1.

    logical, optional, intent(IN) :: per  !< Flag indicating if data are considered periodic
    type(UnivSpline), intent(OUT) :: tck !<  Coefficients, knots, and additional information for interpolation/fitting.
    !! On output the following values will be set:
    !!    c: coefficients
    !!    t: knots
    !!    k: order of splines.
    !!    wrk, iwrk: workspace for internal use only
    !!    fp: weighted sum of squared residuals
    !!
    !! On input the user may provide values for:
    !!    wrk, iwrk: For tasks -1, 1. (usually set by a previous call)
    ! integer, optional, intent(OUT) :: ier !<

    real(dp), dimension(size(w)), pointer :: w_
    real(dp), dimension(:), allocatable :: c_
    real(dp), dimension(:), allocatable :: t_
    real(dp), dimension(:), allocatable :: wrk_
    integer, dimension(:), allocatable :: iwrk_
    real(dp) :: ub_, ue_, s_
    integer :: task_, k_, ier_
    logical :: per_

    integer :: m
    integer :: nest
    integer :: nknots, k1, nwrk, n

    ! Is it periodic data
    per_ = .FALSE.; if (Present(per)) per_ = per

    m = size(x)
    IF (size(y) /= m) call print_msg('size(y) different than size(x)='//str(m), 'splprep')

    ! Weights
    if (Present(w)) then
      IF (size(w) /= m) call print_msg('size(w) different than size(x)='//str(m), 'splprep')
      w_ = w
      IF (.not. Present(s)) s_ = m - sqrt(2._dp * m)
    else
      ! allocate (w_(m))
      w_ = 1._dp
      IF (.not. Present(s)) s_ = 0._dp
    end if

    if (Present(s)) s_ = s

    ! Order of spline
    k_ = 3; IF (Present(k)) k_ = k
    k1 = k_ + 1

    ! Limits
    ub_ = x(1); IF (Present(ub)) ub_ = ub
    ue_ = x(m); IF (Present(ue)) ue_ = ue

    task_ = 0; IF (Present(task)) task_ = task

    IF (Present(t)) then        ! overwrite task input
      task_ = -1  ! If knots given, then task = -1
    end IF

    if (task_ == -1) then
      ! Interior knots t are given. Copy to work array
      IF (.not. Present(t)) call print_msg('Knots are required for task = -1')
      nknots = size(t)
      nest = nknots + 2 * k_ + 2
      allocate (t_(nest))
      t_(k1 + 1:nest - k1) = t
    else if (task_ == 0) then
      ! Reserve memory for knots t. Not initialized
      if (per_) then
        nest = max(m + 2 * k_, 2 * k1 + 1)
      else
        nest = max(m + k1, 2 * k1 + 1)
      end if
      allocate (t_(nest))
    end if

    allocate (c_(nest))

    ! Set workspace
    if (task_ <= 0) then
      if (per_) then
        nwrk = m * k1 + nest * (3 + 5 * k1)
      else
        nwrk = m * k1 + nest * (4 + 3 * k1)
      end if
      allocate (wrk_(nwrk))
      allocate (iwrk_(nest))

    else
      wrk_ = tck%wrk
      nwrk = size(wrk_)
      iwrk_ = tck%iwrk
    end if

    if (per_) then
      call percur(task_, m, x, y, w_, k_, s_, nest, n, t_, c_, tck%fp, wrk_, nwrk, iwrk_, ier_)
    else
      call curfit(task_, m, x, y, w_, ub_, ue_, k_, s_, nest, n, t_, c_, tck%fp, wrk_, nwrk, iwrk_, ier_)
    end if

    tck%ier = ier_
    ! if (Present(ier)) ier = ier_

    tck%k = k_

    ! Set the first n (number of knots) values to the output
    allocate (tck%wrk(n), tck%iwrk(n))
    tck%wrk = wrk_(:n)
    tck%iwrk = iwrk_(:n)

    allocate (tck%t(n), tck%c(n))
    tck%c = c_(:n)
    tck%t = t_(:n)

    deallocate (c_, t_, wrk_, iwrk_)

  end subroutine splprep

  !> splrep Computes
  !!
  !! Examples:
  !!
  ! subroutine splrep(x, y, w, xb, xe, k, task, s, t, per, tck, ier)
  subroutine splrep(x, y, w, xb, xe, k, task, s, t, per, tck)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !< Values of independent variable
    real(dp), dimension(size(x)), intent(IN) :: y !< The data points y=f(x)
    real(dp), optional, dimension(size(x)), intent(IN) :: w !< Strictly positive rank-1 array of weights the same size as x and y.
    !!The weights are used in computing the weighted least-squares spline fit. If the errors in the y values have standard-deviation
    !!given by the vector d, then w should be 1/d.
    real(dp), optional, intent(IN) :: xb !< Lower limit of the interval to approximate (xb <= x(1))
    real(dp), optional, intent(IN) :: xe !< Upper limit of the interval to approximate (xe >= x(size(x))).
    integer, optional, intent(IN) :: k !< The degree of the spline fit. It is recommended to use cubic splines.
    !! Even values of k should be avoided especially with small s values. 1 <= k <= 5
    integer, optional, intent(IN) :: task !< {1, 0, -1},
    !! If task==0 find t and c for a given smoothing factor, s.
    !!
    !! If task==1 find t and c for another value of the smoothing factor, s.
    !! There must have been a previous call with task=0 or task=1 for the same set of data (it will be stored an used internally)
    !!
    !!If task==-1 find the weighted least square spline for a given set of knots t.
    !!    These should be interior knots as knots on the ends will be added automatically
    real(dp), optional, intent(IN) :: s !< A smoothing condition. The amount of smoothness is determined by satisfying the conditions:
    !!
    !! `sum((w * (y - g))**2) <= s` where `g(x)` is the smoothed interpolation of `(x,y)`.
    !!
    !!The user can use `s` to control the tradeoff between closeness and smoothness of fit. Larger `s` means more smoothing while
    !! smaller values of `s` indicate less smoothing.  Recommended values of s depend on the weights, w.
    !!
    !!If the weights represent the inverse of the standard-deviation of y, then a good s value should be found in the range
    !!(m-sqrt(2*m),m+sqrt(2*m)) where m is the number of datapoints in x, y, and w. default : \f$ s= m-\sqrt(2 m)\f$ if weights are
    !!supplied. `s = 0.0` (interpolating) if no weights are supplied.

    real(dp), optional, dimension(:), intent(IN) :: t !< Input knots (interior knots only). if task = -1. If given then task is automatically set to -1.

    logical, optional, intent(IN) :: per  !< Flag indicating if data are considered periodic
    type(UnivSpline), intent(OUT) :: tck !<  Coefficients, knots, and additional information for interpolation/fitting.
    !! On output the following values will be set:
    !!    c: coefficients
    !!    t: knots
    !!    k: order of splines.
    !!    wrk, iwrk: workspace for internal use only
    !!    fp: weighted sum of squared residuals
    !!
    !! On input the user may provide values for:
    !!    wrk, iwrk: For tasks -1, 1. (usually set by a previous call)
    ! integer, optional, intent(OUT) :: ier !<

    real(dp), dimension(size(x)) :: w_
    real(dp), dimension(:), allocatable :: c_
    real(dp), dimension(:), allocatable :: t_
    real(dp), dimension(:), allocatable :: wrk_
    integer, dimension(:), allocatable :: iwrk_
    real(dp) :: xb_, xe_, s_
    integer :: task_, k_, ier_
    logical :: per_

    integer :: m
    integer :: nest
    integer :: nknots, k1, nwrk, n

    ! Is it periodic data
    per_ = .FALSE.; if (Present(per)) per_ = per

    m = size(x)
    IF (size(y) /= m) call print_msg('size(y) different than size(x)='//str(m), 'splrep')

    ! Weights
    if (Present(w)) then
      IF (size(w) /= m) call print_msg('size(w) different than size(x)='//str(m), 'splrep')
      w_ = w
      IF (.not. Present(s)) s_ = m - sqrt(2._dp * m)
    else
      ! allocate (w_(m))
      w_ = 1._dp
      IF (.not. Present(s)) s_ = 0._dp
    end if

    if (Present(s)) s_ = s

    ! Order of spline
    k_ = 3; IF (Present(k)) k_ = k
    k1 = k_ + 1

    ! Limits
    xb_ = x(1); IF (Present(xb)) xb_ = xb
    xe_ = x(m); IF (Present(xe)) xe_ = xe

    task_ = 0; IF (Present(task)) task_ = task

    IF (Present(t)) then        ! overwrite task input
      task_ = -1  ! If knots given, then task = -1
    end IF

    if (task_ == -1) then
      ! Interior knots t are given. Copy to work array
      IF (.not. Present(t)) call print_msg('Knots are required for task = -1')
      nknots = size(t)
      nest = nknots + 2 * k_ + 2
      allocate (t_(nest))
      t_(k1 + 1:nest - k1) = t
    else if (task_ == 0) then
      ! Reserve memory for knots t. Not initialized
      if (per_) then
        nest = max(m + 2 * k_, 2 * k1 + 1)
      else
        nest = max(m + k1, 2 * k1 + 1)
      end if
      allocate (t_(nest))
    end if

    allocate (c_(nest))

    ! Set workspace
    if (task_ <= 0) then
      if (per_) then
        nwrk = m * k1 + nest * (3 + 5 * k1)
      else
        nwrk = m * k1 + nest * (4 + 3 * k1)
      end if
      allocate (wrk_(nwrk))
      allocate (iwrk_(nest))

    else
      wrk_ = tck%wrk
      nwrk = size(wrk_)
      iwrk_ = tck%iwrk
    end if

    if (per_) then
      call percur(task_, m, x, y, w_, k_, s_, nest, n, t_, c_, tck%fp, wrk_, nwrk, iwrk_, ier_)
    else
      call curfit(task_, m, x, y, w_, xb_, xe_, k_, s_, nest, n, t_, c_, tck%fp, wrk_, nwrk, iwrk_, ier_)
    end if

    tck%ier = ier_
    ! if (Present(ier)) ier = ier_

    tck%k = k_

    ! Set the first n (number of knots) values to the output
    allocate (tck%wrk(n), tck%iwrk(n))
    tck%wrk = wrk_(:n)
    tck%iwrk = iwrk_(:n)

    allocate (tck%t(n), tck%c(n))
    tck%c = c_(:n)
    tck%t = t_(:n)

    deallocate (c_, t_, wrk_, iwrk_)

  end subroutine splrep

  !> splev Computes a B-spline or its derivatives.
  !! Given the knots and coefficients of a B-spline representation, evaluate
  !! the value of the smoothing polynomial and its derivatives.  This is a
  !! wrapper around the FORTRAN routines splev and splder of FITPACK.
  !!
  !! Examples:
  !!
  function splev(x, tck, der, ext, ier) result(y)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !< Points at which to return the value of the smoothed spline or its derivative
    type(UnivSpline), intent(IN) :: tck !< A spline representation returned by splrep
    integer, optional, intent(IN) :: der !< The order of the derivative of the spline to compute (must be less than k).
    integer, optional, intent(IN) :: ext !< Flag controling the result for  ``x`` outside the
    integer, optional, intent(OUT) :: ier !< Flag given output status
    !! interval defined by the knot sequence.
    real(dp), dimension(size(x)) :: y !< Smoothed or interpolated spline values

    integer :: k, n, m, e, ier_, nu
    k = tck%k
    n = size(tck%t)
    m = size(x)
    e = 0; IF (Present(ext)) e = ext
    nu = 0; IF (Present(der)) nu = der

    IF ((nu < 0) .or. (nu > k)) call print_msg('Must be 0 <= der='//str(nu)//' <= k='//str(k))

    IF ((e < 0) .or. (e > 3)) call print_msg('ext = '//str(e)//' not one of (0, 1, 2, 3)')
    call splder(tck%t, n, tck%c, k, nu, x, y, m, e, tck%wrk, ier_)

    IF (ier_ == 10) call print_msg('Invalid input data', 'splev', errcode=0)
    IF (ier_ == 1) call print_msg('x value out of bounds and e == 2', 'splev', errcode=0)
    IF (Present(ier)) ier = ier_
  end function splev

end module fitpack
