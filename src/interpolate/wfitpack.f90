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
  end type UnivSpline

  Private
  Public :: UnivSpline, splrep, splprep, splev

  character, private, parameter :: nl = new_line('a')

  interface splev
    module procedure :: splevp, splevc
  end interface splev

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
  !!
  !! Examples:
  !! --------
  !!
  !! ```
  !!  real(dp), dimension(:), allocatable :: phi
  !!  real(dp), dimension(:), allocatable :: r
  !!  real(dp), dimension(:, :), allocatable :: x
  !!  real(dp), dimension(:, :), allocatable :: new_points
  !!  real(dp), dimension(:), allocatable :: u
  !!  type(UnivSpline) :: tck
  !!  character(len=:), allocatable :: header
  !!  real(dp) :: s = 0._dp
  !!  ! Generate a discretization of a limacon curve in the polar coordinates:
  !!  phi = linspace(Zero, 2 * M_PI, Nd)
  !!  r = 0.5_8 + cos(phi)        ! polar coords
  !!  allocate (r(Nd), u(Nd), phi(Nd))
  !!  allocate (x(2, Nd), new_points(2, Nd))
  !!  x(1, :) = r * cos(phi)      ! convert to cartesian
  !!  x(2, :) = r * sin(phi)      ! convert to cartesian
  !!  ! interpolate
  !!  call splprep(x, u, tck, s=s)
  !!  call splevp(u, tck, new_points)
  !!  ! and write to stdout
  !!  header = "u                    x                   y"
  !!  call save_array([u, new_points(1, :), new_points(2, :)], 3, fmt="f12.8", header=header)
  !!  ! Prints:
  !!  !
  !!  ! # u                    x                   y
  !!  ! 0.00000000000000  1.50000000000000  0.00000000000000
  !!  ! 0.03616230734577  1.46779335229253  0.23853963732962
  !!  ! 0.07211603380768  1.37398960267532  0.45870512901973
  !!  ! 0.10765440708122  1.22676038619218  0.64385351896871
  !!  ! 0.14257428865429  1.03883011365998  0.78063018793468
  !!  ! 0.17667803554905  0.82622920670009  0.86019572275769
  !!  ! 0.20977542560675  0.60672992984431  0.87900005428954
  !!  ! ....
  !!  ! ....
  !!  ! 0.79022457439325  0.60672992984431 -0.87900005428954
  !!  ! 0.82332196445095  0.82622920670009 -0.86019572275769
  !!  ! 0.85742571134571  1.03883011365998 -0.78063018793468
  !!  ! 0.89234559291878  1.22676038619218 -0.64385351896871
  !!  ! 0.92788396619232  1.37398960267532 -0.45870512901973
  !!  ! 0.96383769265423  1.46779335229253 -0.23853963732962
  !!  ! 1.00000000000000  1.50000000000000 -0.00000000000000
  !!  !
  !!```
  !!
  !!  Notice that (i) we force interpolation by using `s=0`,
  !!  (ii) the parameterization, ``u``, is generated automatically.
  subroutine splprep(x, u, tck, w, ulim, k, task, upar, s, t, per, ier)
    implicit none
    real(dp), dimension(:, :), intent(IN) :: x !< 2D-Array representing the curve in an n-dimensional space
    real(dp), dimension(:), intent(INOUT) :: u !< An array with parameters. The routine will fill it if upar is False or not present.
    real(dp), optional, dimension(size(x(1, :))), target, intent(IN) :: w !< Strictly positive rank-1 array of weights the same size as u.

    !!The weights are used in computing the weighted least-squares spline fit. If the errors in the y values have standard-deviation
    !!given by the vector d, then w should be 1/d.
    real(dp), dimension(2), optional, intent(INOUT) :: ulim !< Lower and upper limits of the interval to approximate (ulim(1) <= u(1) and ulim(2) >= u(m) )
    integer, optional, intent(IN) :: k !< The degree of the spline fit. It is recommended to use cubic splines.
    !! Even values of k should be avoided especially with small s values. 1 <= k <= 5
    integer, optional, intent(IN) :: task !< {1, 0, -1},\\
    !! - If `task==0` find t and c for a given smoothing factor, s.
    !! - If `task==1` find t and c for another value of the smoothing factor, s.
    !! There must have been a previous call with task=0 or task=1 for the same set of data (it will be stored an used internally)
    !! - If `task==-1` find the weighted least square spline for a given set of knots t.
    !!    These should be interior knots as knots on the ends will be added automatically
    real(dp), optional, intent(IN) :: s !< A smoothing condition. The amount of smoothness is determined by satisfying the conditions:\n
    !! `sum((w * (y - g))**2) <= s` where `g(x)` is the smoothed interpolation of `(x,y)`.\n
    !!The user can use `s` to control the tradeoff between closeness and smoothness of fit. Larger `s` means more smoothing while
    !! smaller values of `s` indicate less smoothing.  Recommended values of s depend on the weights, w.\n
    !!If the weights represent the inverse of the standard-deviation of y, then a good s value should be found in the range
    !!\f$ (m-\sqrt{2 m},m+\sqrt{2 m})\f$ where m is the number of datapoints in x, y, and w. default : \f$ s= m-\sqrt{2 m}\f$ if weights are
    !!supplied. `s = 0.0` (interpolating) if no weights are supplied.

    real(dp), optional, dimension(:), intent(IN) :: t !< Input knots (interior knots only) needed for task = -1. If given then task is automatically set to -1.

    logical, optional, intent(IN) :: upar  !< Flag indicating if u is given or must be automatically calculated. Default `.False.`
    logical, optional, intent(IN) :: per  !< Flag indicating if data are considered periodic
    type(UnivSpline), intent(OUT) :: tck !<  Coefficients, knots, and additional information for interpolation/fitting.\n
    !! On output the following values of tck will be set:
    !! - c: coefficients
    !! - t: knots
    !! - k: order of splines.
    !! - wrk, iwrk: workspace for internal use only
    !! - fp: weighted sum of squared residuals
    !!
    !! For tasks -1, or 1, on input the user may provide values for:
    !!    wrk and iwrk: (but are usually set by a previous call)
    integer, optional, intent(OUT) :: ier !< Error code

    real(dp) :: tol
    integer :: i, ia1, ia2, ib, ifp, ig1, ig2, iq, iz, maxit, ncc

    real(dp), dimension(:), allocatable :: w_
    real(dp), dimension(:), allocatable :: c_
    real(dp), dimension(:), allocatable :: t_
    real(dp), dimension(:), allocatable :: wrk_
    integer, dimension(:), allocatable :: iwrk_
    real(dp) :: ub, ue, s_, Du
    integer :: task_, k_, ier_
    logical :: per_, upar_
    integer :: m, idim
    integer :: nest
    integer :: mx, k1, k2, nwrk, n
    integer, dimension(2) :: shx
    ! integer :: pos

    ! Set up the parameters tol and maxit
    tol = 0.001_8
    maxit = 20

    shx = shape(x)
    idim = shx(1)
    m = shx(2)
    mx = idim * m

    ! !!!!!!!!!! Checks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check dimension of space where lives the curve
    IF (idim < 1 .or. idim > 10) call print_msg('0 < idim < 11 must hold', errcode=1)

    ! Check periodicity
    per_ = .FALSE.; if (Present(per)) per_ = per
    if (per_) then
      if (any(x(:, 1) /= x(:, m))) then ! Check that it really is periodic, otherwise fix it
        call print_msg('Warning: Setting data periodic', 'splprep', errcode=0)
        ! x(:, m) = x(:, 1)
      end if
    end if

    ! Check parameter array u
    IF (size(u) < m) call print_msg('size(u) < size(m)', 'splprep', errcode=1)
    upar_ = .False.; IF (Present(upar)) upar_ = upar

    ! Check weights and smoothing factor
    if (Present(w)) then
      IF (size(w) /= m) call print_msg('size(w) different than size(x)='//str(m), 'splprep', 10)
      IF (any(w <= 0._8)) call print_msg('Weights values must be non-negative', 'splprep', 10)
      w_ = w
      IF (.not. Present(s)) s_ = m - sqrt(2._dp * m)
    else
      allocate (w_(m))
      w_ = 1._dp
      IF (.not. Present(s)) s_ = 0._dp
    end if
    IF (Present(s)) s_ = s
    IF (s_ < 0._8) call print_msg('Smooth factor s='//str(s_)//'must be non-negative', 'splprep', 10)

    ! Check order of spline
    k_ = 3; IF (Present(k)) k_ = k
    IF (k_ < 0 .or. k_ > 5) call print_msg('1 <= k= '//str(k_)//' <= 5 must hold', 'splprep', 10)
    IF (m <= k_) call print_msg('k= '//str(k_)//' < m ='//str(m)//' must hold', 'splprep', 10)
    k1 = k_ + 1
    k2 = k1 + 1

    ! Check the task
    task_ = 0; IF (Present(task)) task_ = task
    !! Scipy version do not force the task if t is present in this routine. We do (for now at least)
    ! Check knots
    IF (Present(t)) task_ = -1 ! If knots given, then task = -1
    IF (abs(task_) > 1) call print_msg('task = '//str(task_)//', should be -1, 0 or 1', 'splprep', 10)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!      Set working array     !!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check knots (according to task desired)
    if (task_ == -1) then      ! Interior knots t are given. Copy to work array
      IF (.not. Present(t)) call print_msg('Knots are required for task = -1', errcode=10)
      nest = size(t) + 2 * k_ + 2
      allocate (t_(nest))
      t_(k1 + 1:nest - k1) = t
    else if (task_ == 0) then      ! Reserve memory for knots t. Not initialized
      if (per_) then
        nest = m + 2 * k_
      else
        nest = m + k_ + 1
      end if
      allocate (t_(nest))
      t_ = 0._8
    else                        ! All knots are given
      nest = size(t)
      allocate (t_(nest))
      t_ = t
    end if

    ncc = nest * idim
    allocate (c_(ncc))

    ! Set workspace
    if (task_ <= 0) then        ! Not previously defined
      if (per_) then
        nwrk = m * k1 + nest * (2 + idim + 5 * k1)
      else
        nwrk = m * k1 + nest * (3 + idim + 3 * k1)
      end if

      IF ((allocated(wrk_)) .and. (size(wrk_) /= nwrk)) deallocate (wrk_)
      IF (.not. allocated(wrk_)) allocate (wrk_(nwrk))

      IF ((allocated(iwrk_)) .and. (size(iwrk_) /= nest)) deallocate (iwrk_)
      IF (.not. allocated(iwrk_)) allocate (iwrk_(nest))
    else                        ! Use from previous call
      wrk_ = tck%wrk
      iwrk_ = tck%iwrk
    end if

    ! Set the parameter of the curve: "u", and the limits
    if (upar_) then
      if (Present(ulim)) then
        ub = ulim(1); ue = ulim(2)
      else
        ub = u(1); ue = u(m)
      end if
    else if (task_ <= 0) then    ! If upar_ is False and task_=-1,0  we calculate the values of u automatically
      u(1) = 0._8
      do i = 2, m
        u(i) = u(i - 1) + norm2(x(:, i) - x(:, i - 1))
      end do
      IF (u(m) <= 0._8) call print_msg("u(m)="//str(u(m))//" <= 0", errcode=1)
      u = u / u(m)
      ub = 0._8
      ue = 1._8
    end if

    IF (ub > u(1)) call print_msg('ulim(1)='//str(ub)//' <= u(1)'//str(u(1))//' must hold', errcode=10)
    IF (ue < u(m)) call print_msg('ulim(2)='//str(ue)//' >= u(m)'//str(u(m))//' must hold', errcode=10)
    IF (any(u(1:m - 1) >= u(2:m))) call print_msg('u must be in ascending order', errcode=10)

    n = nest
    if (task_ == -1) then      ! t(k+2:k-n-1) are provided by the user. We set the exterior knots
      if (per_) then
        Du = u(m) - u(1)
        t_(k1) = u(1); t_(n - k) = u(m) ! First two values
        t_(:k_) = t_(n - 2 * k_:n - k1) - Du
        t_(n - k_ + 1:n) = t_(k1 + 1:k1 + k_) + Du
        call fpchep(u, m, t_, n, k_, ier_)
      else
        t_(1:k1) = ub
        t_(n - k:n) = ue
        call fpchec(u, m, t_, n, k_, ier_)
      end if

      IF (ier /= 0) &
        & call print_msg('Knots not positioned correctly for task=-1', errcode=ier)
    else
      IF ((per_ .and. nest < (m + 2 * k_)) .or. ((.not. per_) .and. nest < (m + k1)))&
        & call print_msg('Too few knots provided', 'splprep', errcode=10)
    end if

    ! We partition the working space and determine the spline curve.
    ier_ = 0

    ifp = 1
    iz = ifp + nest
    ia1 = iz + ncc
    ia2 = ia1 + nest * k1

    if (per_) then              ! periodic data
      ib = ia2 + nest * k
      ig1 = ib + nest * k2
      ig2 = ig1 + nest * k2
      iq = ig2 + nest * k1
      call fpclos(task_, idim, m, u, mx, x, w_, k_, s_, nest, tol, maxit, k1, k2, n, t_,&
        & ncc, c_, tck%fp, wrk_(ifp), wrk_(iz), wrk_(ia1), wrk_(ia2), wrk_(ib), wrk_(ig1),&
        &  wrk_(ig2), wrk_(iq), iwrk_, ier_)
    else
      ig1 = ia2 + nest * k2
      iq = ig1 + nest * k2
      call fppara(task_, idim, m, u, mx, x, w_, ub, ue, k_, s_, nest, tol, maxit, k1, k2,&
        & n, t_, ncc, c_, tck%fp, wrk_(ifp:iz - 1), wrk_(iz), wrk_(ia1:ia2 - 1), &
        &wrk_(ia2:ig1 - 1), wrk_(ig1:iq - 1), wrk_(iq:), iwrk_, ier_)
    end if

    if (Present(ier)) ier = ier_
    tck%k = k_

    ! Save workspaces for possible future use
    tck%wrk = wrk_
    tck%iwrk = iwrk_

    ! Set the first n (number of knots) values to the output
    tck%c = c_(:idim * n)
    tck%t = t_(:n)

    deallocate (c_)
    deallocate (w_)
    deallocate (t_)
    deallocate (wrk_)
    deallocate (iwrk_)
  end subroutine splprep

  !> splrep Finds the B-spline representation of 1-D curve.
  !! Given the set of data points ``(x[i], y[i])`` determine a smooth spline
  !! approximation of degree k on the interval ``xb <= x <= xe``.
  !!
  !!
  !! Examples:
  !! --------
  !!```
  !!  integer, parameter :: N = 6
  !!  integer, parameter :: Nnew = 29
  !!  real(dp), dimension(N) :: x
  !!  real(dp), dimension(N) :: y
  !!  real(dp), dimension(Nnew) :: xnew
  !!  real(dp), dimension(Nnew) :: ynew
  !!  character(len=:), allocatable :: header
  !!  real(dp) :: s
  !!  type(UnivSpline) :: tck
  !!  ! Generate data
  !!  x = linspace(Zero, M_PI, N)
  !!  y = sin(x)
  !!  ! The new array where evaluate the spline
  !!  xnew = linspace(Zero, M_PI, Nnew)
  !!  ! interpolate
  !!  call splrep(x, y, tck=tck, s=s)
  !!  call splev(xnew, tck, ynew)
  !!  ! and write to stdout
  !!  header = "x                   y"
  !!  call save_array([xnew, ynew], 2, fmt="f17.14", header=header)
  !!  ! Prints:
  !!  !
  !!  ! 0.00000000000000  0.00000000000000
  !!  ! 0.11219973762821  0.11401345780438
  !!  ! 0.22439947525641  0.22522911715938
  !!  ! 0.33659921288462  0.33269883705930
  !!  ! 0.44879895051283  0.43547447649845
  !!  ! 0.56099868814103  0.53260789447113
  !!  !  ...
  !!  !  ...
  !!  ! 2.46839422782055  0.62315094997166
  !!  ! 2.58059396544876  0.53260789447113
  !!  ! 2.69279370307697  0.43547447649845
  !!  ! 2.80499344070517  0.33269883705930
  !!  ! 2.91719317833338  0.22522911715938
  !!  ! 3.02939291596159  0.11401345780438
  !!  ! 3.14159265358979  0.00000000000000
  !!  !
  !!```
  !!
  subroutine splrep(x, y, w, xb, xe, k, task, s, t, per, tck, ier)
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
    !! - If task==0 find t and c for a given smoothing factor, s.
    !! - If task==1 find t and c for another value of the smoothing factor, s.
    !! There must have been a previous call with task=0 or task=1 for the same set of data (it will be stored an used internally)
    !! - If task==-1 find the weighted least square spline for a given set of knots t.
    !! These should be interior knots, as knots on the ends will be added automatically
    real(dp), optional, intent(IN) :: s !< A smoothing condition. The amount of smoothness is determined by satisfying the conditions:\n
    !! `sum((w * (y - g))**2) <= s` where `g(x)` is the smoothed interpolation of `(x,y)`.
    !!The user can use `s` to control the tradeoff between closeness and smoothness of fit. Larger `s` means more smoothing while
    !! smaller values of `s` indicate less smoothing.  Recommended values of s depend on the weights, w.\n
    !!If the weights represent the inverse of the standard-deviation of y, then a good s value should be found in the range
    !!\f$ (m-\sqrt{2 m},m+\sqrt{2 m}) where m is the number of datapoints in x, y, and w. default : \f$ s= m-\sqrt{2 m}\f$ if weights are supplied. It will use `s = 0.0` (interpolating) if no weights are supplied.
    real(dp), optional, dimension(:), intent(IN) :: t !< Input knots (interior knots only) for task = -1. If given then task is automatically set to -1.
    logical, optional, intent(IN) :: per  !< Flag indicating if data are considered periodic
    type(UnivSpline), intent(OUT) :: tck !<  Coefficients, knots, and additional information for interpolation/fitting.
    !! On output the following values will be set:
    !!   - c: coefficients
    !!   - t: knots
    !!   - k: order of splines.
    !!   - wrk, iwrk: workspace arrays, for internal use only
    !!   - fp: weighted sum of squared residuals
    !! On input, for tasks -1, 1 the user may provide values for:
    !!    wrk, iwrk: For tasks. (but usually are set by a previous call)
    integer, optional, intent(OUT) :: ier !< Flag of status at output

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
      ! allocate (w_(m))
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

    IF (Present(ier)) ier = ier_

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
  subroutine splevc(x, tck, y, der, ext, ier)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !< Points at which to return the value of the smoothed spline or its derivative
    type(UnivSpline), intent(IN) :: tck !< A spline representation returned by splrep
    real(dp), dimension(size(x)), intent(OUT) :: y !< Smoothed or interpolated spline values
    integer, optional, intent(IN) :: der !< The order of the derivative of the spline to compute (must be less than k).
    integer, optional, intent(IN) :: ext !< Flag controling the result for  ``x`` outside the
    integer, optional, intent(OUT) :: ier !< Flag given output status
    !! interval defined by the knot sequence.
    real(dp), dimension(:), allocatable :: wrk_

    integer :: k, n, m, e, ier_, nu

    k = tck%k
    n = size(tck%t)
    m = size(x)
    e = 0; IF (Present(ext)) e = ext
    nu = 0; IF (Present(der)) nu = der

    allocate (wrk_(size(tck%wrk)))
    wrk_ = tck%wrk
    IF ((nu < 0) .or. (nu > k)) call print_msg('Must be 0 <= der='//str(nu)//' <= k='//str(k))
    IF ((e < 0) .or. (e > 3)) call print_msg('ext = '//str(e)//' not one of (0, 1, 2, 3)')
    !
    call splder(tck%t, n, tck%c, k, nu, x, y, m, e, wrk_, ier_)
    IF (ier_ == 10) call print_msg('Invalid input data', 'splev', errcode=0)
    IF (ier_ == 1) call print_msg('x value out of bounds and e == 2', 'splev', errcode=0)
    IF (Present(ier)) ier = ier_
  end subroutine splevc

  !> splev Computes a B-spline or its derivatives.
  !! Given the knots and coefficients of a B-spline representation, evaluate
  !! the value of the smoothing polynomial and its derivatives.  This is a
  !! wrapper around the FORTRAN routines splev and splder of FITPACK.
  !!
  !! Examples:
  !!
  subroutine splevp(u, tck, y, der, ext, ier)
    implicit none
    real(dp), dimension(:), intent(IN) :: u !< Points at which to return the value of the smoothed spline or its derivative
    type(UnivSpline), intent(IN) :: tck !< A spline representation returned by splrep
    integer, optional, intent(IN) :: der !< The order of the derivative of the spline to compute (must be less than k).
    integer, optional, intent(IN) :: ext !< Flag controling the result for  ``x`` outside the
    ! integer, optional, intent(IN) :: idim !< Number of dimensions for parametric curve

    integer, optional, intent(OUT) :: ier !< Flag given output status
    real(dp), dimension(:), allocatable :: wrk_
    ! Result
    real(dp), dimension(:, :), intent(OUT) :: y !< Smoothed or interpolated spline values

    integer :: k, n, m, e, ier_, nu
    integer :: idim_             ! If parametric, dimension of curve
    integer :: shy(2)
    integer :: i, pos

    ! Dimension of the curve

    k = tck%k
    n = size(tck%t)
    m = size(u)
    e = 0; IF (Present(ext)) e = ext
    nu = 0; IF (Present(der)) nu = der
    shy = shape(y)

    idim_ = shy(1)
    IF (shy(2) /= m) call print_msg('Number of points in y='//str(shy(2))//' different from number of parameters u'//str(m))

    allocate (wrk_(size(tck%wrk)))
    wrk_ = tck%wrk
    IF ((nu < 0) .or. (nu > k)) call print_msg('Must be 0 <= der='//str(nu)//' <= k='//str(k))

    IF ((e < 0) .or. (e > 3)) call print_msg('ext = '//str(e)//' not one of (0, 1, 2, 3)')
    !

    do i = 1, idim_
      pos = (i - 1) * n
      call splder(tck%t, n, tck%c(pos + 1:pos + n), k, nu, u, y(i, :), m, e, wrk_, ier_)
      IF (ier_ == 10) call print_msg('Invalid input data', 'splev', errcode=0)
      IF (ier_ == 1) call print_msg('u value out of bounds and e == 2', 'splev', errcode=0)
    end do
    IF (Present(ier)) ier = ier_

  end subroutine splevp

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We include the necessary fitpack routines
  include "fitp.inc"

end module fitpack
