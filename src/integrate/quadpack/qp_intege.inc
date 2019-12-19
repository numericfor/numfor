  ! PRNM(qage) (f, a, b, epsabs, epsrel, key, result, abserr, neval, ier, &
  !  & limit, alist, blist, rlist, elist, iord, last)
  PRNM(qage) (f, a, b, epsabs, epsrel, rule, result, abserr, neval, ier, &
   & limit, alist, blist, rlist, elist, iord, last)
!*****************************************************************************80
!! QAGE estimates a definite integral.
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Parameters:
!
!    Input, external real(dp) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(dp) :: f
!      real(dp) :: x
!    which evaluates the integrand function.
!
!    Input, real(dp) :: A, B, the limits of integration.
!
!    Input, real(dp) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer KEY, chooses the order of the local integration rule:
!    1,  7 Gauss points, 15 Gauss-Kronrod points,
!    2, 10 Gauss points, 21 Gauss-Kronrod points,
!    3, 15 Gauss points, 31 Gauss-Kronrod points,
!    4, 20 Gauss points, 41 Gauss-Kronrod points,
!    5, 25 Gauss points, 51 Gauss-Kronrod points,
!    6, 30 Gauss points, 61 Gauss-Kronrod points.
!
!    Input, integer LIMIT, the maximum number of subintervals that
!    can be used.
!
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!
!    Output, real(dp) :: ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!    Output, integer IER, return code.
!    0, normal and reliable termination of the routine.  It is assumed that the
!      requested accuracy has been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can
!      allow more subdivisions by increasing the value of LIMIT in QAG.
!      However, if this yields no improvement it is advised to analyze the
!      integrand to determine the integration difficulties.  If the position
!      of a local difficulty can be determined, such as a singularity or
!      discontinuity within the interval) one will probably gain from
!      splitting up the interval at this point and calling the integrator
!      on the subranges.  If possible, an appropriate special-purpose
!      integrator should be used which is designed for handling the type
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
!
!    Workspace, real(dp) :: ALIST(LIMIT), BLIST(LIMIT), contains in entries 1
!    through LAST the left and right ends of the partition subintervals.
!
!    Workspace, real(dp) :: RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.
!
!    Workspace, real(dp) :: ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.
!
!    Output, integer IORD(LIMIT), the first K elements of which are pointers
!    to the error estimates over the subintervals, such that
!    elist(iord(1)), ..., elist(iord(k)) form a decreasing sequence, with
!    k = last if last <= (limit/2+2), and k = limit+1-last otherwise.
!
!    Output, integer LAST, the number of subintervals actually produced
!    in the subdivision process.
!
!  Local parameters:
!    alist     - list of left end points of all subintervals
!                       considered up to now
!    blist     - list of right end points of all subintervals
!                       considered up to now
!    elist(i)  - error estimate applying to rlist(i)
!    maxerr    - pointer to the interval with largest error estimate
!    errmax    - elist(maxerr)
!    area      - sum of the integrals over the subintervals
!    errsum    - sum of the errors over the subintervals
!    errbnd    - requested accuracy max(epsabs,epsrel*abs(result))
!    *****1    - variable for the left subinterval
!    *****2    - variable for the right subinterval
!    last      - index for subdivision
  NUMFOR_KINDF :: f          !< function to integrate
  NUMFOR_KINDR, intent(OUT) :: result

  real(dp), intent(IN) :: a !<
  real(dp), intent(IN) :: b !<
  real(dp), intent(IN) :: epsabs !<
  real(dp), intent(IN) :: epsrel !<
  character(len=4), optional, intent(IN) :: rule      !< rule to use
  ! integer, intent(IN), optional :: key !<
  real(dp), intent(OUT) :: abserr !<
  integer, intent(OUT) :: neval !<
  integer, intent(OUT) :: ier !<
  integer, intent(IN) :: limit !<
  real(dp), dimension(limit), intent(OUT) :: alist !<
  real(dp), dimension(limit), intent(OUT) :: blist !<
  NUMFOR_KINDR, dimension(limit), intent(OUT) :: rlist !<
  real(dp), dimension(limit), intent(OUT) :: elist !<
  integer, dimension(limit), intent(OUT) :: iord !<
  integer, intent(OUT) :: last !<

!! Examples:
!!
  NUMFOR_KINDR :: area
  NUMFOR_KINDR :: area1
  NUMFOR_KINDR :: area12
  NUMFOR_KINDR :: area2
  real(dp) :: a1
  real(dp) :: a2
  real(dp) :: b1
  real(dp) :: b2
  real(dp) :: c
  real(dp) :: defabs
  real(dp) :: defab1
  real(dp) :: defab2
  real(dp) :: errbnd
  real(dp) :: errmax
  real(dp) :: error1
  real(dp) :: error2
  real(dp) :: erro12
  real(dp) :: errsum

  character(len=4) :: rule_
  integer :: Nrule

  integer iroff1
  integer iroff2
  integer maxerr
  integer nrmax
  real(dp) :: resabs

!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0_8
  abserr = 0.0_8
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0_8
  elist(1) = 0.0_8
  iord(1) = 0
  if (epsabs < 0.0_8 .and. epsrel < 0.0_8) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  rule_ = 'qk21'; IF (Present(rule)) rule_ = rule
  read (rule(3:4), *) Nrule

  ! keyf = max(1, min(key, 6))
  ! if (keyf == 1) then
  !   Nrule = 15
  ! else
  !   Nrule = keyf * 10 + 1
  ! end if

  neval = 0
  ! call qgk(f, a, b, result, qgkrule(keyf), abserr, defabs, resabs)
  call qgk(f, a, b, result, rule_, abserr, defabs, resabs)

  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
!
!  Test on accuracy.
!
  errbnd = max(epsabs, epsrel * abs(result))
  IF (abserr <= 50 * epsilon(defabs) * defabs .and. errbnd < abserr) ier = 2
  IF (limit == 1) ier = 1

  if (ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs) .or. abserr == 0) then
    neval = Nrule * (2 * neval + 1)
    return
  end if
!
!  Initialization.
!
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  iroff1 = 0
  iroff2 = 0
  do last = 2, limit
    !
    !  Bisect the subinterval with the largest error estimate.
    !
    a1 = alist(maxerr)
    b1 = 0.5_8 * (alist(maxerr) + blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)

    ! call qgk(f, a1, b1, area1, qgkrule(keyf), error1, resabs, defab1)
    ! call qgk(f, a2, b2, area2, qgkrule(keyf), error2, resabs, defab2)
    call qgk(f, a1, b1, area1, rule_, error1, resabs, defab1)
    call qgk(f, a2, b2, area2, rule_, error2, resabs, defab2)

    !
    !  Improve previous approximations to integral and error and
    !  test for accuracy.
    !
    neval = neval + 1
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)
    if (defab1 /= error1 .and. defab2 /= error2) then
      IF (abs(rlist(maxerr) - area12) <= 1.e-5_8 * abs(area12) .and. 0.99_8 * errmax <= erro12)&
        & iroff1 = iroff1 + 1
      IF (10 < last .and. errmax < erro12) iroff2 = iroff2 + 1
    end if
    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max(epsabs, epsrel * abs(area))
    !
    !  Test for roundoff error and eventually set error flag.
    !
    if (errbnd < errsum) then
      IF (6 <= iroff1 .or. 20 <= iroff2) ier = 2
      !
      !  Set error flag in the case that the number of subintervals
      !  equals limit.
      IF (last == limit) ier = 1
      !
      !  Set error flag in the case of bad integrand behavior
      !  at a point of the integration range.
      IF (max(abs(a1), abs(b2)) <= (1.0_8 + c * 1.0e3_8 * epsilon(a1)) * &
        & (abs(a2) + 10000 * tiny(a2))) ier = 3
    end if
    !
    !  Append the newly-created intervals to the list.
    !
    if (error2 <= error1) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with the largest error estimate (to be bisected next).
    call qsort(limit, last, maxerr, errmax, elist, iord, nrmax)
    if (ier /= 0 .or. errsum <= errbnd) then
      exit
    end if
  end do
!
!  Compute final result.
  result = sum(rlist(1:last))
  abserr = errsum
  neval = Nrule * (2 * neval + 1)
  end PRNM(qage)

  ! ========================================================================
!*****************************************************************************80
!! QAGS estimates the integral of a function.
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!
!  Parameters:
!
!    Input, external real(dp) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(dp) :: f
!      real(dp) :: x
!    which evaluates the integrand function.
!
!    Input, real(dp) :: A, B, the limits of integration.
!
!    Input, real(dp) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!
!    Output, real(dp) :: ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!    Output, integer IER, error flag.
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the data value of
!                             limit in qags (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement
!                             it is advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a
!                             local difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the sub-
!                             ranges. if possible, an appropriate
!                             special-purpose integrator should be used,
!                             which is designed for handling the type
!                             of difficulty involved.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some  points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved, and that the
!                             returned result is the best which can be
!                             obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             result, abserr and neval are set to zero.
!
!  Local Parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2 containing
!                       the part of the epsilon table which is still
!                       needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left interval
!           *****2    - variable for the right interval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true value)
!
  PRNM(qagse) (f, a, b, epsabs, epsrel, rule, result, abserr, neval, ier, &
    & limit, alist, blist, rlist, elist, iord, last)
  implicit none
  NUMFOR_KINDF :: f         !< function to integrate
  real(dp), intent(IN) :: a       !< lower limit
  real(dp), intent(IN) :: b       !< upper limit
  real(dp), intent(IN) :: epsabs  !< Desired absolute error
  real(dp), intent(IN) :: epsrel  !< Desired relative error
  character(len=4), optional, intent(IN) :: rule      !< rule to use
  NUMFOR_KINDR, intent(OUT) :: result !< Estimation of integral
  real(dp), intent(OUT) :: abserr !< Estimation of error
  integer, intent(OUT) :: neval   !< Number of function evaluations
  integer, intent(OUT) :: ier     !< Error code
!
  integer, INTENT(IN) :: limit
  real(dp), intent(OUT) :: alist(limit)
  real(dp), intent(OUT) :: blist(limit)
  real(dp), intent(OUT) :: elist(limit)
  NUMFOR_KINDR, intent(OUT) :: rlist(limit)
  integer, intent(OUT) :: iord(limit)
  integer, intent(OUT)::  last

  real(dp) :: abseps
  NUMFOR_KINDR :: area, area1, area12, area2
  real(dp) :: a1, a2
  real(dp) :: b1, b2
  real(dp) :: correc
  real(dp) :: defabs, defab1, defab2, dres

  real(dp) :: erlarg, erlast, errbnd, errmax
  real(dp) :: error1, error2, erro12, errsum, ertest
  logical extrap, noext
  integer id, k
  integer ierro, iroff1, iroff2, iroff3
  integer jupbnd, ksgn, ktmin, maxerr, nres, nrmax
  integer numrl2

  real(dp) :: resabs
  NUMFOR_KINDR :: reseps
  NUMFOR_KINDR :: res3la(3)
  NUMFOR_KINDR :: rlist2(52)
  real(dp) :: small
  character(len=4) :: rule_
  ! integer :: Nrule, keyf
  integer :: Nrule
!
!  The dimension of rlist2 is determined by the value of
!  limexp in QEXTR (rlist2 should be of dimension
!  (limexp+2) at least).
!
!  Test on validity of parameters.
  ier = 0; neval = 0; last = 0; result = 0._8; abserr = 0._8
  alist(1) = a; blist(1) = b; rlist(1) = 0._8; elist(1) = 0._8
  if (epsabs < epsilon(0._8) .and. (epsrel < 50 * epsilon(1._8))) then
    ier = 6
    return
  end if

  rule_ = 'qk21'; IF (Present(rule)) rule_ = rule
  read (rule(3:4), *) Nrule
!
!  First approximation to the integral.
  ierro = 0

  ! call qgk(f, a, b, result, qgkrule(keyf), abserr, defabs, resabs)
  call qgk(f, a, b, result, rule_, abserr, defabs, resabs)
!
!  Test on accuracy.
  dres = abs(result)
  errbnd = max(epsabs, epsrel * dres)
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  IF (abserr <= 100 * epsilon(defabs) * defabs .and. abserr > errbnd) ier = 2
  IF (limit == 1) ier = 1
  IF (ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs) .or. abserr == 0._8) go to 140
!
!  Initialization.
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge(abserr)
  nrmax = 1; nres = 0
  numrl2 = 2; ktmin = 0
  extrap = .false.; noext = .false.
  iroff1 = 0; iroff2 = 0; iroff3 = 0
  if (dres >= (1._8 - 50 * epsilon(defabs)) * defabs) then
    ksgn = 1
  else
    ksgn = -1
  end if
  do last = 2, limit
    !  Bisect the subinterval with the nrmax-th largest error estimate.
    a1 = alist(maxerr)
    b1 = 0.5_8 * (alist(maxerr) + blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    ! call qgk(f, a1, b1, area1, qgkrule(keyf), error1, resabs, defab1)
    ! call qgk(f, a2, b2, area2, qgkrule(keyf), error2, resabs, defab2)
    call qgk(f, a1, b1, area1, rule_, error1, resabs, defab1)
    call qgk(f, a2, b2, area2, rule_, error2, resabs, defab2)

    !
    !  Improve previous approximations to integral and error
    !  and test for accuracy.
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)
    IF (defab1 == error1 .or. defab2 == error2) go to 15
    IF (abs(rlist(maxerr) - area12) > 1.e-5_8 * abs(area12)&
      & .or. (erro12 < 0.99_8 * errmax)) go to 10
    if (extrap) then
      iroff2 = iroff2 + 1
    else
      iroff1 = iroff1 + 1
    end if
10  continue
    IF (last > 10 .and. erro12 > errmax) iroff3 = iroff3 + 1
15  continue
    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max(epsabs, epsrel * abs(area))
    !
    !  Test for roundoff error and eventually set error flag.
    IF (iroff1 + iroff2 >= 10 .or. iroff3 >= 20) ier = 2
    IF (iroff2 >= 5) ierro = 3
    !
    !  Set error flag in the case that the number of subintervals
    !  equals limit.
    IF (last == limit) ier = 1
    !
    !  Set error flag in the case of bad integrand behavior
    !  at a point of the integration range.
    !
    IF (max(abs(a1), abs(b2)) <= (1._8 + 1000._8 * epsilon(a1)) * &
        (abs(a2) + 1000._8 * tiny(a2))) ier = 4
    !
    !  Append the newly-created intervals to the list.
    if (error2 <= error1) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with nrmax-th largest error estimate (to be bisected next).
    call qsort(limit, last, maxerr, errmax, elist, iord, nrmax)
    if (errsum <= errbnd) go to 115
    if (ier /= 0) then
      exit
    end if
    IF (last == 2) go to 80
    IF (noext) go to 90
    erlarg = erlarg - erlast
    IF (abs(b1 - a1) > small) erlarg = erlarg + erro12
    !
    !  Test whether the interval to be bisected next is the
    !  smallest interval.
    if (.not. extrap) then
      if (abs(blist(maxerr) - alist(maxerr)) > small) go to 90
      extrap = .true.
      nrmax = 2
    end if
    !40  continue
    !
    !  The smallest interval has the largest error.
    !  Before bisecting decrease the sum of the errors over the
    !  larger intervals (erlarg) and perform extrapolation.
    !
    if (ierro /= 3 .and. erlarg > ertest) then
      id = nrmax
      jupbnd = last
      IF (last > (2 + limit / 2)) jupbnd = limit + 3 - last
      do k = id, jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        IF (abs(blist(maxerr) - alist(maxerr)) > small) go to 90
        nrmax = nrmax + 1
      end do
    end if
    !
    !  Perform extrapolation.
    numrl2 = numrl2 + 1
    rlist2(numrl2) = area
    call qextr(numrl2, rlist2, reseps, abseps, res3la, nres)
    ktmin = ktmin + 1
    IF (ktmin > 5 .and. abserr < 1.e-3_8 * errsum) ier = 5

    if (abseps < abserr) then
      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max(epsabs, epsrel * abs(reseps))
      IF (abserr <= ertest) exit
    end if
    !
    !  Prepare bisection of the smallest interval.
    IF (numrl2 == 1) noext = .true.
    IF (ier == 5) exit

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 0.5_8
    erlarg = errsum
    go to 90
80  continue
    small = abs(b - a) * 0.375_8
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area
90  continue
  end do
!
!  Set final result and error estimate.
!
  IF (abserr == huge(abserr)) go to 115
  IF (ier + ierro == 0) go to 110
  IF (ierro == 3) abserr = abserr + correc
  IF (ier == 0) ier = 3
  IF (result /= 0._8 .and. area /= 0._8) go to 105
  IF (abserr > errsum) go to 115
  IF (area == 0._8) go to 130
  go to 110
105 continue
  IF (abserr / abs(result) > errsum / abs(area)) go to 115
!
!  Test on divergence.
110 continue
  IF (ksgn == (-1) .and. max(abs(result), abs(area)) <= defabs * 0.01_8) go to 130
  IF (0.01_8 > abs(result / area) .or. abs(result / area) > 100._8 &
    & .or. (errsum > abs(area))) ier = 6

  go to 130
!
!  Compute global integral sum.
115 continue
  result = sum(rlist(1:last))
  abserr = errsum
130 continue
  IF (2 < ier) ier = ier - 1

140 continue
  neval = 2 * Nrule * last - Nrule
  end PRNM(qagse)

  ! ========================================================================
!> QAGP computes a definite integral.
!!
!! The routine calculates an approximation RESULT to a definite integral
!!   I = integral of F over (A,B),
!! hopefully satisfying
!!   || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!!
!! Interior break points of the integration interval,
!! where local difficulties of the integrand may occur, such as
!! singularities or discontinuities, are provided by the user.
  PRNM(qagpe) (f, a, b, npts2, points, epsabs, epsrel, result, abserr, neval, ier, &
    & limit, alist, blist, rlist, elist, iord, last)
!*****************************************************************************80
!  Parameters:
!
!    Input, external real(dp) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(dp) :: f
!      real(dp) :: x
!    which evaluates the integrand function.
!
!    Input, real(dp) :: A, B, the limits of integration.
!
!    Input, integer NPTS2, the number of user-supplied break points within
!    the integration range, plus 2.  NPTS2 must be at least 2.
!
!    Input/output, real(dp) :: POINTS(NPTS2), contains the user provided interior
!    breakpoints in entries 1 through NPTS2-2.  If these points are not
!    in ascending order on input, they will be sorted.
!
!    Input, real(dp) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!
!    Output, real(dp) :: ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!    Output, integer IER, return flag.
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the data value
!                             of limit in qagp(and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement
!                             it is advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a local
!                             difficulty can be determined (i.e.
!                             singularity, discontinuity within the
!                             interval), it should be supplied to the
!                             routine as an element of the vector
!                             points. if necessary, an appropriate
!                             special-purpose integrator must be used,
!                             which is designed for handling the type
!                             of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved, and that
!                             the returned result is the best which
!                             can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier > 0.
!                         = 6 the input is invalid because
!                             npts2 < 2 or
!                             break points are specified outside
!                             the integration range or
!                             epsabs < 0 and epsrel < 0,
!                             or limit < npts2.
!                             result, abserr, neval are set to zero.
!
!  Local parameters:
!
!            the dimension of rlist2 is determined by the value of
!            limexp in QEXTR (rlist2 should be of dimension
!            (limexp+2) at least).
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2
!                       containing the part of the epsilon table which
!                       is still needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements in rlist2. if an appropriate
!                       approximation to the compounded integral has
!                       obtained, it is put in rlist2(numrl2) after
!                       numrl2 has been increased by one.
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine
!                       is attempting to perform extrapolation. i.e.
!                       before subdividing the smallest interval we
!                       try to decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation is
!                       no longer allowed (true-value)
!
  implicit none
  integer :: limit
  real(dp) :: a
  real(dp) :: abseps
  real(dp) :: abserr
  real(dp) :: alist(limit)
  NUMFOR_KINDR :: area, area1, area12, area2
  real(dp) :: a1, a2
  real(dp) :: b
  real(dp) :: blist(limit)
  real(dp) :: b1, b2
  real(dp) :: correc
  real(dp) :: defabs, defab1, defab2, dres
  real(dp) :: elist(limit)
  real(dp) :: epsabs
  real(dp) :: epsrel
  real(dp) :: erlarg, erlast, errbnd, errmax, error1, erro12, error2, errsum, ertest
  logical extrap, noext
  NUMFOR_KINDF :: f
  integer i, id, ier, ierro
  integer ind1, ind2
  integer iord(limit)
  integer iroff1, iroff2, iroff3
  integer j, jlow, jupbnd
  integer k, ksgn, ktmin
  integer last, levcur, levmax
  integer level(limit)
  integer maxerr
  integer ndin(40)
  integer neval
  integer nint, npts, npts2, nres, nrmax, numrl2
  real(dp)  :: points(npts2 - 2)
  real(dp)  :: pts(npts2)
  real(dp)  :: resa, resabs
  NUMFOR_KINDR :: result
  NUMFOR_KINDR :: reseps
  NUMFOR_KINDR :: res3la(3)
  NUMFOR_KINDR :: rlist(limit)
  NUMFOR_KINDR :: rlist2(52)
  real(dp) :: sign, temp
!
!  Test on validity of parameters.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0_8
  abserr = 0.0_8
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0_8
  elist(1) = 0.0_8
  iord(1) = 0
  level(1) = 0
  npts = npts2 - 2
  IF ((npts2 < 2) .or. limit <= npts .or. (epsabs < 0.0_8 .and. epsrel < 0.0_8)) return

!
!  If any break points are provided, sort them into an
!  ascending sequence.
  if (b < a) then
    sign = -1.0e+00
  else
    sign = +1.0E+00
  end if
  pts(1) = min(a, b)
  pts(2:npts + 1) = points(1:npts)
  pts(npts + 2) = max(a, b)
  nint = npts + 1
  a1 = pts(1)
  if (npts /= 0) then
    do i = 1, nint
      do j = i + 1, nint + 1
        if (pts(j) < pts(i)) then
          temp = pts(i)
          pts(i) = pts(j)
          pts(j) = temp
        end if
      end do
    end do
    IF (pts(1) /= min(a, b) .or. pts(nint + 1) /= max(a, b)) return
  end if
  ier = 0

!
!  Compute first integral and error approximations.
  resabs = 0.0_8
  do i = 1, nint
    b1 = pts(i + 1)
    call qk21(f, a1, b1, area1, error1, defabs, resa)
    abserr = abserr + error1
    result = result + area1
    ndin(i) = 0
    if (error1 == resa .and. error1 /= 0.0_8) then
      ndin(i) = 1
    end if
    resabs = resabs + defabs
    level(i) = 0
    elist(i) = error1
    alist(i) = a1
    blist(i) = b1
    rlist(i) = area1
    iord(i) = i
    a1 = b1
  end do
  errsum = 0.0_8
  do i = 1, nint
    if (ndin(i) == 1) then
      elist(i) = abserr
    end if
    errsum = errsum + elist(i)
  end do
!
!  Test on accuracy.
!
  last = nint
  neval = 21 * nint
  dres = abs(result)
  errbnd = max(epsabs, epsrel * dres)
  IF (abserr <= 100 * epsilon(resabs) * resabs .and. abserr > errbnd) ier = 2

  if (nint /= 1) then
    do i = 1, npts
      jlow = i + 1
      ind1 = iord(i)
      do j = jlow, nint
        ind2 = iord(j)
        if (elist(ind1) <= elist(ind2)) then
          ind1 = ind2
          k = j
        end if
      end do
      if (ind1 /= iord(i)) then
        iord(k) = iord(i)
        iord(i) = ind1
      end if
    end do
    IF (limit < npts2) ier = 1
  end if
  IF (ier /= 0 .or. abserr <= errbnd) return

!
!  Initialization
!
  rlist2(1) = result
  maxerr = iord(1)
  errmax = elist(maxerr)
  area = result
  nrmax = 1
  nres = 0
  numrl2 = 1
  ktmin = 0
  extrap = .False.
  noext = .False.
  erlarg = errsum
  ertest = errbnd
  levmax = 1
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ierro = 0
  abserr = huge(abserr)
  if (dres >= (1.0_8 - 0.5_8 * epsilon(resabs)) * resabs) then
    ksgn = 1
  else
    ksgn = -1
  end if
  do last = npts2, limit
    !
    !  Bisect the subinterval with the NRMAX-th largest error estimate.
    !
    levcur = level(maxerr) + 1
    a1 = alist(maxerr)
    b1 = 0.5_8 * (alist(maxerr) + blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk21(f, a1, b1, area1, error1, resa, defab1)
    call qk21(f, a2, b2, area2, error2, resa, defab2)
    !
    !  Improve previous approximations to integral and error
    !  and test for accuracy.
    !
    neval = neval + 42
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)
    if (defab1 /= error1 .and. defab2 /= error2) then
      if (abs(rlist(maxerr) - area12) <= 1.0e-05 * abs(area12) .and. erro12 >= 0.99_8 * errmax) then
        if (extrap) then
          iroff2 = iroff2 + 1
        else
          iroff1 = iroff1 + 1
        end if
      end if
      IF (last > 10 .and. erro12 > errmax) iroff3 = iroff3 + 1

    end if
    level(maxerr) = levcur
    level(last) = levcur
    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max(epsabs, epsrel * abs(area))
    !
    !  Test for roundoff error and eventually set error flag.
    IF (10 <= iroff1 + iroff2 .or. 20 <= iroff3) ier = 2
    IF (5 <= iroff2) ierro = 3
    !
    !  Set error flag in the case that the number of subintervals
    !  equals limit.
    if (last == limit) ier = 1
    !
    !  Set error flag in the case of bad integrand behavior
    !  at a point of the integration range
    IF (max(abs(a1), abs(b2)) <= (1.0e+00 + 1.0e+03 * epsilon(a1)) * &
        (abs(a2) + 1.0e+03 * tiny(a2))) ier = 4
    !
    !  Append the newly-created intervals to the list.
    if (error2 <= error1) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with nrmax-th largest error estimate (to be bisected next).
    call qsort(limit, last, maxerr, errmax, elist, iord, nrmax)
    IF (errsum <= errbnd) go to 190
    IF (ier /= 0) exit
    IF (noext) cycle
    erlarg = erlarg - erlast
    IF (levcur + 1 <= levmax) erlarg = erlarg + erro12
    !
    !  Test whether the interval to be bisected next is the
    !  smallest interval.
    if (.not. extrap) then
      IF (level(maxerr) + 1 <= levmax) cycle
      extrap = .true.
      nrmax = 2
    end if
    !
    !  The smallest interval has the largest error.
    !  Before bisecting decrease the sum of the errors over the
    !  larger intervals (erlarg) and perform extrapolation.
    if (ierro /= 3 .and. erlarg > ertest) then
      id = nrmax
      jupbnd = last
      IF (last > (2 + limit / 2)) jupbnd = limit + 3 - last
      do k = id, jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        IF (level(maxerr) + 1 <= levmax) go to 160
        nrmax = nrmax + 1
      end do
    end if
    !
    !  Perform extrapolation.
    numrl2 = numrl2 + 1
    rlist2(numrl2) = area
    IF (numrl2 <= 2) go to 155
    call qextr(numrl2, rlist2, reseps, abseps, res3la, nres)
    ktmin = ktmin + 1
    IF (5 < ktmin .and. abserr < 1.0e-03 * errsum) ier = 5
    if (abseps < abserr) then
      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max(epsabs, epsrel * abs(reseps))
      IF (abserr < ertest) exit
    end if
    !
    !  Prepare bisection of the smallest interval.
    IF (numrl2 == 1) noext = .true.
    IF (5 <= ier) exit

155 continue
    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    levmax = levmax + 1
    erlarg = errsum
160 continue
  end do
!
!  Set the final result.
  IF (abserr == huge(abserr)) go to 190
  IF ((ier + ierro) == 0) go to 180
  IF (ierro == 3) abserr = abserr + correc
  IF (ier == 0) ier = 3
  IF (result /= 0.0_8 .and. area /= 0.0_8) go to 175
  IF (errsum < abserr) go to 190
  IF (area == 0.0_8) go to 210

  go to 180
175 continue
  IF (abserr / abs(result) > errsum / abs(area)) go to 190

!
!  Test on divergence.
180 continue
  IF (ksgn == (-1) .and. max(abs(result), abs(area)) <= resabs * 1.0e-02) go to 210
  IF (0.01_8 > abs(result / area) .or. abs(result / area) > 100._8 .or. errsum > abs(area)) ier = 6

  go to 210
!
!  Compute global integral sum.
190 continue
  result = sum(rlist(1:last))
  abserr = errsum
210 continue
  IF (2 < ier) ier = ier - 1
  result = result * sign

  end PRNM(qagpe)

  ! ========================================================================
!! QAWSE estimates integrals with algebraico-logarithmic endpoint singularities.
!!
!! This routine calculates an approximation RESULT to an integral
!!   I = integral of F(X) * W(X) over (a,b),
!! where W(X) shows a singular behavior at the endpoints, hopefully
!! satisfying:
!!   | I - RESULT | <= max ( epsabs, epsrel * |I| ).
  PRNM(qawse) (f, a, b, alfa, beta, flgw, epsabs, epsrel, limit, &
               result, abserr, neval, ier, alist, blist, rlist, elist, iord, last)
!*****************************************************************************80
!  Parameters:
!    Input, external real(dp) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(dp) :: f
!      real(dp) :: x
!    which evaluates the integrand function.
!
!    Input, real(dp) :: A, B, the limits of integration.
!
!    Input, real(dp) :: ALFA, BETA, parameters used in the weight function.
!    ALFA and BETA should be greater than -1.
!
!    Input, integer FLGW, indicates which weight function is used:
!    = 1  (x-a)**alfa*(b-x)**beta
!    = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
!    = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
!    = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
!
!    Input, real(dp) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer LIMIT, an upper bound on the number of subintervals
!    in the partition of (A,B), LIMIT >= 2.  If LIMIT < 2, the routine
!     will end with IER = 6.
!
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!
!    Output, real(dp) :: ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for the integral and error
!                             are less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit. however, if this yields no
!                             improvement it is advised to analyze the
!                             integrand, in order to determine the
!                             integration difficulties which prevent
!                             the requested tolerance from being
!                             achieved. in case of a jump discontinuity
!                             or a local singularity of algebraico-
!                             logarithmic type at one or more interior
!                             points of the integration range, one
!                             should proceed by splitting up the
!                             interval at these points and calling the
!                             integrator on the subranges.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behavior occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             b <= a or alfa <= (-1) or beta <= (-1) or
!                             flgw < 1 or flgw > 4, or
!                             epsabs < 0 and epsrel < 0,
!                             or limit < 2.
!                             result, abserr, neval, rlist(1), elist(1),
!                             iord(1) and last are set to zero.
!                             alist(1) and blist(1) are set to a and b
!                             respectively.
!
!    Workspace, real(dp) :: ALIST(LIMIT), BLIST(LIMIT), contains in entries 1
!    through LAST the left and right ends of the partition subintervals.
!
!    Workspace, real(dp) :: RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.
!
!    Workspace, real(dp) :: ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.
!
!            iord   - integer
!                     vector of dimension at least limit, the first k
!                     elements of which are pointers to the error
!                     estimates over the subintervals, so that
!                     elist(iord(1)), ..., elist(iord(k)) with k = last
!                     if last <= (limit/2+2), and k = limit+1-last
!                     otherwise, form a decreasing sequence.
!
!    Output, integer LAST, the number of subintervals actually produced in
!    the subdivision process.
!
!  Local parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!
  implicit none
  NUMFOR_KINDF :: f
  NUMFOR_KINDR, intent(OUT) :: result

  integer limit
  real(dp) :: a
  real(dp) :: abserr
  real(dp) :: alfa
  real(dp) :: alist(limit)
  NUMFOR_KINDR :: area
  NUMFOR_KINDR :: area1
  NUMFOR_KINDR :: area12
  NUMFOR_KINDR :: area2
  real(dp) :: a1
  real(dp) :: a2
  real(dp) :: b
  real(dp) :: beta
  real(dp) :: blist(limit)
  real(dp) :: b1
  real(dp) :: b2
  real(dp) :: centre
  real(dp) :: elist(limit)
  real(dp) :: epsabs
  real(dp) :: epsrel
  real(dp) :: errbnd
  real(dp) :: errmax
  real(dp) :: error1
  real(dp) :: erro12
  real(dp) :: error2
  real(dp) :: errsum
  integer ier
  integer flgw
  integer iord(limit)
  integer iroff1
  integer iroff2
  integer last
  integer maxerr
  integer nev
  integer neval
  integer nrmax
  real(dp) :: resas1
  real(dp) :: resas2
  real(dp) :: rg(25)
  real(dp) :: rh(25)
  real(dp) :: ri(25)
  real(dp) :: rj(25)
  NUMFOR_KINDR :: rlist(limit)
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0

  rlist(1) = 0.0_8
  elist(1) = 0.0_8
  iord(1) = 0
  result = 0.0_8
  abserr = 0.0_8
  if (b <= a .or. (epsabs < 0.0_8 .and. epsrel < 0.0_8) .or. &
      (alfa <= -1.0_8) .or. (beta <= -1.0_8) .or. &
      (flgw < 1) .or. (flgw > 4) .or. limit < 2) then
    ier = 6
    return
  end if
!
!  Compute the modified Chebyshev moments.
!
  call qmomo(alfa, beta, ri, rj, rg, rh, flgw)
!
!  Integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
!
  centre = 0.5_8 * (b + a)
  call qc25s(f, a, centre, [a, b, alfa, beta], ri, rj, rg, rh, area1, &
             error1, resas1, flgw, nev)
  neval = nev
  call qc25s(f, centre, b, [a, b, alfa, beta], ri, rj, rg, rh, area2, &
             error2, resas2, flgw, nev)
  last = 2
  neval = neval + nev
  result = area1 + area2
  abserr = error1 + error2
!
!  Test on accuracy.
!
  errbnd = max(epsabs, epsrel * abs(result))
!
!  Initialization.
!
  if (error2 <= error1) then
    alist(1) = a
    alist(2) = centre
    blist(1) = centre
    blist(2) = b
    rlist(1) = area1
    rlist(2) = area2
    elist(1) = error1
    elist(2) = error2
  else
    alist(1) = centre
    alist(2) = a
    blist(1) = b
    blist(2) = centre
    rlist(1) = area2
    rlist(2) = area1
    elist(1) = error2
    elist(2) = error1
  end if
  iord(1) = 1
  iord(2) = 2
  if (limit == 2) then
    ier = 1
    return
  end if
  IF (abserr <= errbnd) return

  errmax = elist(1)
  maxerr = 1
  nrmax = 1
  area = result
  errsum = abserr
  iroff1 = 0
  iroff2 = 0
  do last = 3, limit
    !
    !  Bisect the subinterval with largest error estimate.
    !
    a1 = alist(maxerr)
    b1 = 0.5_8 * (alist(maxerr) + blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)
    call qc25s(f, a1, b1, [a, b, alfa, beta], ri, rj, rg, rh, area1, &
               error1, resas1, flgw, nev)
    neval = neval + nev
    call qc25s(f, a2, b2, [a, b, alfa, beta], ri, rj, rg, rh, area2, &
               error2, resas2, flgw, nev)
    neval = neval + nev
    !
    !  Improve previous approximations integral and error and
    !  test for accuracy.
    !
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)
    !
    !  Test for roundoff error.
    if (a /= a1 .and. b /= b2) then
      if (resas1 /= error1 .and. resas2 /= error2) then
        IF (abs(rlist(maxerr) - area12) < 1.0e-05 * abs(area12) .and. erro12 >= 9.9e-01 * errmax) &
          &        iroff1 = iroff1 + 1
        IF (last > 10 .and. erro12 > errmax) iroff2 = iroff2 + 1

      end if
    end if
    rlist(maxerr) = area1
    rlist(last) = area2
    !
    !  Test on accuracy.
    errbnd = max(epsabs, epsrel * abs(area))
    if (errsum > errbnd) then
      !  Set error flag in the case that the number of interval
      !  bisections exceeds limit.
      IF (last == limit) ier = 1
      !  Set error flag in the case of roundoff error.
      IF (iroff1 >= 6 .or. iroff2 >= 20) ier = 2
      !  Set error flag in the case of bad integrand behavior
      !  at interior points of integration range.
      IF (max(abs(a1), abs(b2)) <= (1.0e+00 + 1.0e+03 * epsilon(a1)) * &
          (abs(a2) + 1000._dp * tiny(a2))) ier = 3
    end if
    !
    !  Append the newly-created intervals to the list.
    !
    if (error2 <= error1) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with largest error estimate (to be bisected next).
    !
    call qsort(limit, last, maxerr, errmax, elist, iord, nrmax)
    if (ier /= 0 .or. errsum <= errbnd) then
      exit
    end if
  end do
!
!  Compute final result.
!
  result = sum(rlist(1:last))
  abserr = errsum

  end PRNM(qawse)
  ! ========================================================================
!> QAWCE computes a Cauchy principal value.
!!
!! The routine calculates an approximation RESULT to a Cauchy principal
!! value
!!   I = integral of F*W over (A,B),
!! with
!!   W(X) = 1 / ( X - C ),
!! with C distinct from A and B, hopefully satisfying
!!   | I - RESULT | <= max ( EPSABS, EPSREL * |I| ).
  PRNM(qawce) (f, a, b, c, epsabs, epsrel, limit, result, abserr, neval, &
               ier, alist, blist, rlist, elist, iord, last)
!*****************************************************************************80
!  Parameters:
!
!    Input, external real(dp) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(dp) :: f
!      real(dp) :: x
!    which evaluates the integrand function.
!
!    Input, real(dp) :: A, B, the limits of integration.
!
!    Input, real(dp) :: C, a parameter in the weight function, which cannot be
!    equal to A or B.
!
!    Input, real(dp) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer LIMIT, the upper bound on the number of subintervals that
!    will be used in the partition of [A,B].  LIMIT is typically 500.
!
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!
!    Output, real(dp) :: ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the value of
!                             limit. however, if this yields no
!                             improvement it is advised to analyze the
!                             integrand, in order to determine the
!                             integration difficulties.  if the position
!                             of a local difficulty can be determined
!                             (e.g. singularity, discontinuity within
!                             the interval) one will probably gain
!                             from splitting up the interval at this
!                             point and calling appropriate integrators
!                             on the subranges.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behavior occurs
!                             at some interior points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             c = a or c = b or
!                             epsabs < 0 and epsrel < 0,
!                             or limit < 1.
!                             result, abserr, neval, rlist(1), elist(1),
!                             iord(1) and last are set to zero.
!                             alist(1) and blist(1) are set to a and b
!                             respectively.
!
!    Workspace, real(dp) :: ALIST(LIMIT), BLIST(LIMIT), contains in entries 1
!    through LAST the left and right ends of the partition subintervals.
!
!    Workspace, real(dp) :: RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.
!
!    Workspace, real(dp) :: ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.
!
!            iord    - integer
!                      vector of dimension at least limit, the first k
!                      elements of which are pointers to the error
!                      estimates over the subintervals, so that
!                      elist(iord(1)), ...,  elist(iord(k)) with
!                      k = last if last <= (limit/2+2), and
!                      k = limit+1-last otherwise, form a decreasing
!                      sequence.
!
!            last    - integer
!                      number of subintervals actually produced in
!                      the subdivision process
!
!  Local parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!
  implicit none
  NUMFOR_KINDF :: f          !< function to integrate
  NUMFOR_KINDR, intent(OUT) :: result

  real(dp), intent(IN) :: a !<
  real(dp), intent(IN) :: b !<
  real(dp), intent(IN) :: c !<
  real(dp), intent(IN) :: epsabs !<
  real(dp), intent(IN) :: epsrel !<
  integer, intent(IN) :: limit !<
  real(dp), intent(OUT) :: abserr !<
  integer, intent(OUT) :: neval !<
  integer, intent(OUT) :: ier !<
  real(dp), dimension(limit), intent(OUT) :: alist !<
  real(dp), dimension(limit), intent(OUT) :: blist !<
  NUMFOR_KINDR, dimension(limit), intent(OUT) :: rlist !<
  real(dp), dimension(limit), intent(OUT) :: elist !<
  integer, dimension(limit), intent(OUT) :: iord !<
  integer, intent(OUT) :: last !<
!! Examples:
!!
  real(dp) :: aa
  NUMFOR_KINDR :: area
  NUMFOR_KINDR :: area1
  NUMFOR_KINDR :: area12
  NUMFOR_KINDR :: area2
  real(dp) :: a1
  real(dp) :: a2
  real(dp) :: bb
  real(dp) :: b1
  real(dp) :: b2
  real(dp) :: errbnd
  real(dp) :: errmax
  real(dp) :: error1
  real(dp) :: error2
  real(dp) :: erro12
  real(dp) :: errsum
  integer iroff1
  integer iroff2
  integer krule
  integer maxerr
  integer nev
  integer nrmax
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0_8
  elist(1) = 0.0_8
  iord(1) = 0
  result = 0.0_8
  abserr = 0.0_8
  if (c == a) then
    ier = 6
    return
  else if (c == b) then
    ier = 6
    return
  else if (epsabs < 0.0_8 .and. epsrel < 0.0_8) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  if (a <= b) then
    aa = a
    bb = b
  else
    aa = b
    bb = a
  end if
  krule = 1
  call qc25c(f, aa, bb, c, result, abserr, krule, neval)
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  alist(1) = a
  blist(1) = b
!
!  Test on accuracy.
!
  errbnd = max(epsabs, epsrel * abs(result))
  if (limit == 1) then
    ier = 1
    go to 70
  end if
  IF (abserr < min(0.01_8 * abs(result), errbnd)) go to 70
!
!  Initialization
!
  alist(1) = aa
  blist(1) = bb
  rlist(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  iroff1 = 0
  iroff2 = 0
  do last = 2, limit
    !
    !  Bisect the subinterval with nrmax-th largest error estimate.
    !
    a1 = alist(maxerr)
    b1 = 0.5_8 * (alist(maxerr) + blist(maxerr))
    b2 = blist(maxerr)
    IF (c <= b1 .and. a1 < c) b1 = 0.5_8 * (c + b2)
    IF (b1 < c .and. c < b2) b1 = 0.5_8 * (a1 + c)

    a2 = b1
    krule = 2
    call qc25c(f, a1, b1, c, area1, error1, krule, nev)
    neval = neval + nev
    call qc25c(f, a2, b2, c, area2, error2, krule, nev)
    neval = neval + nev
    !
    !  Improve previous approximations to integral and error
    !  and test for accuracy.
    !
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)
    IF (abs(rlist(maxerr) - area12) < 1.0e-05 * abs(area12) .and. erro12 >= 9.9e-01 * errmax .and. krule == 0) &
      iroff1 = iroff1 + 1
    IF (last > 10 .and. erro12 > errmax .and. krule == 0) iroff2 = iroff2 + 1

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max(epsabs, epsrel * abs(area))
    if (errsum > errbnd) then
      !  Test for roundoff error and eventually set error flag.
      IF (iroff1 >= 6 .and. iroff2 > 20) ier = 2
      !  Set error flag in the case that number of interval
      !  bisections exceeds limit.
      IF (last == limit) ier = 1
      !  Set error flag in the case of bad integrand behavior at
      !  a point of the integration range.
      IF (max(abs(a1), abs(b2)) <= (1.0e+00 + 1.0e+03 * epsilon(a1)) * (abs(a2) + 1.0e+03 * tiny(a2))) ier = 3

    end if
    !
    !  Append the newly-created intervals to the list.
    !
    if (error2 <= error1) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with NRMAX-th largest error estimate (to be bisected next).
    !
    call qsort(limit, last, maxerr, errmax, elist, iord, nrmax)
    IF (ier /= 0 .or. errsum <= errbnd) exit

  end do
!
!  Compute final result.
  result = sum(rlist(1:last))
  abserr = errsum
70 continue
  IF (aa == b) result = -result
  end PRNM(qawce)

  ! ========================================================================
!> qawfe calculates integrate(f(x) cos(w x), a, inf) or integrate(f(x) sin(w x), a, inf)
!!
!! The routine calculates an approximation RESULT to a definite integral
!!   I = integral of F*COS(OMEGA*X) or F*SIN(OMEGA*X) over (A,+Infinity),
!! hopefully satisfying
!!   || I - RESULT || <= EPSABS.
  PRNM(qawfe) (f, a, omega, flgw, epsabs, limlst, limit, maxp1, result, abserr, neval, ier, &
    & rslst, erlst, ierlst, lst, alist, blist, rlist, elist, iord, nnlog, chebmo)
!*****************************************************************************80
!  Parameters:
!
!    Input, external real(dp) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(dp) :: f
!      real(dp) :: x
!    which evaluates the integrand function.
!
!    Input, real(dp) :: A, the lower limit of integration.
!
!    Input, real(dp) :: OMEGA, the parameter in the weight function.
!
!    Input, integer FLGW, indicates which weight function is used
!    = 1      w(x) = cos(omega*x)
!    = 2      w(x) = sin(omega*x)
!
!    Input, real(dp) :: EPSABS, the absolute accuracy requested.
!
!    Input, integer LIMLST, an upper bound on the number of cycles.
!    LIMLST must be at least 1.  In fact, if LIMLST < 3, the routine
!    will end with IER= 6.
!
!    Input, integer LIMIT, an upper bound on the number of subintervals
!    allowed in the partition of each cycle, limit >= 1.
!
!            maxp1  - integer
!                     gives an upper bound on the number of
!                     Chebyshev moments which can be stored, i.e.
!                     for the intervals of lengths abs(b-a)*2**(-l),
!                     l=0,1, ..., maxp1-2, maxp1 >= 1
!
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!
!    Output, real(dp) :: ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!            ier    - ier = 0 normal and reliable termination of
!                             the routine. it is assumed that the
!                             requested accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error
!                             are less reliable. it is assumed that
!                             the requested accuracy has not been
!                             achieved.
!                    if omega /= 0
!                     ier = 6 the input is invalid because
!                             (flgw /= 1 and flgw /= 2) or
!                              epsabs <= 0 or limlst < 3.
!                              result, abserr, neval, lst are set
!                              to zero.
!                         = 7 bad integrand behavior occurs within
!                             one or more of the cycles. location
!                             and type of the difficulty involved
!                             can be determined from the vector ierlst.
!                             here lst is the number of cycles actually
!                             needed (see below).
!                             ierlst(k) = 1 the maximum number of
!                                           subdivisions (= limit)
!                                           has been achieved on the
!                                           k th cycle.
!                                       = 2 occurence of roundoff
!                                           error is detected and
!                                           prevents the tolerance
!                                           imposed on the k th cycle
!                                           from being acheived.
!                                       = 3 extremely bad integrand
!                                           behavior occurs at some
!                                           points of the k th cycle.
!                                       = 4 the integration procedure
!                                           over the k th cycle does
!                                           not converge (to within the
!                                           required accuracy) due to
!                                           roundoff in the
!                                           extrapolation procedure
!                                           invoked on this cycle. it
!                                           is assumed that the result
!                                           on this interval is the
!                                           best which can be obtained.
!                                       = 5 the integral over the k th
!                                           cycle is probably divergent
!                                           or slowly convergent. it
!                                           must be noted that
!                                           divergence can occur with
!                                           any other value of
!                                           ierlst(k).
!                         = 8 maximum number of  cycles  allowed
!                             has been achieved, i.e. of subintervals
!                             (a+(k-1)c,a+kc) where
!                             c = (2*int(abs(omega))+1)*pi/abs(omega),
!                             for k = 1, 2, ..., lst.
!                             one can allow more cycles by increasing
!                             the value of limlst (and taking the
!                             according dimension adjustments into
!                             account).
!                             examine the array iwork which contains
!                             the error flags over the cycles, in order
!                             to eventual look for local integration
!                             difficulties.
!                             if the position of a local difficulty can
!                             be determined (e.g. singularity,
!                             discontinuity within the interval)
!                             one will probably gain from splitting
!                             up the interval at this point and
!                             calling appopriate integrators on the
!                             subranges.
!                         = 9 the extrapolation table constructed for
!                             convergence acceleration of the series
!                             formed by the integral contributions
!                             over the cycles, does not converge to
!                             within the required accuracy.
!                             as in the case of ier = 8, it is advised
!                             to examine the array iwork which contains
!                             the error flags on the cycles.
!                    if omega = 0 and flgw = 1,
!                    the integral is calculated by means of qagi
!                    and ier = ierlst(1) (with meaning as described
!                    for ierlst(k), k = 1).
!
!            rslst  - real
!                     vector of dimension at least limlst
!                     rslst(k) contains the integral contribution
!                     over the interval (a+(k-1)c,a+kc) where
!                     c = (2*int(abs(omega))+1)*pi/abs(omega),
!                     k = 1, 2, ..., lst.
!                     note that, if omega = 0, rslst(1) contains
!                     the value of the integral over (a,infinity).
!
!            erlst  - real
!                     vector of dimension at least limlst
!                     erlst(k) contains the error estimate
!                     corresponding with rslst(k).
!
!            ierlst - integer
!                     vector of dimension at least limlst
!                     ierlst(k) contains the error flag corresponding
!                     with rslst(k). for the meaning of the local error
!                     flags see description of output parameter ier.
!
!            lst    - integer
!                     number of subintervals needed for the integration
!                     if omega = 0 then lst is set to 1.
!
!            alist, blist, rlist, elist - real
!                     vector of dimension at least limit,
!
!            iord, nnlog - integer
!                     vector of dimension at least limit, providing
!                     space for the quantities needed in the
!                     subdivision process of each cycle
!
!            chebmo - real
!                     array of dimension at least (maxp1,25),
!                     providing space for the Chebyshev moments
!                     needed within the cycles
!
!  Local parameters:
!
!           c1, c2    - end points of subinterval (of length
!                       cycle)
!           cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
!           psum      - vector of dimension at least (limexp+2)
!                       (see routine qextr)
!                       psum contains the part of the epsilon table
!                       which is still needed for further computations.
!                       each element of psum is a partial sum of
!                       the series which should sum to the value of
!                       the integral.
!           errsum    - sum of error estimates over the
!                       subintervals, calculated cumulatively
!           epsa      - absolute tolerance requested over current
!                       subinterval
  implicit none
  NUMFOR_KINDF :: f !<
  real(dp), intent(IN) :: a !<
  real(dp), intent(IN) :: omega !<
  integer, intent(IN) :: flgw !<
  real(dp), intent(IN) :: epsabs !< The absolute accuracy requested
  integer, intent(IN) :: limlst !< An upper bound on the number of cycles (>=3).
  integer, intent(IN) :: limit !< An upper bound on the number of subintervals
!! allowed in the partition of each cycle, limit >= 1.
  integer, intent(IN) :: maxp1 !< Upper bound on the number of Chebyshev moments which can be stored
  NUMFOR_KINDR, intent(OUT) :: result !< Result of integration
  real(dp), intent(OUT) :: abserr !< Estimate of || I - result ||
  integer, intent(OUT) :: neval !< Number of function evaluations
  integer, intent(OUT) :: ier !< Output code
  NUMFOR_KINDR, dimension(limlst), intent(OUT) :: rslst !< rslst(k) contains the integral contribution over the interval
!! (a+(k-1)c,a+kc) where c = (2*int(abs(omega))+1)*pi/abs(omega),
  real(dp), dimension(limlst), intent(OUT) :: erlst !< erlst(k) contains the error estimate corresponding to rslst(k).
  integer, dimension(limlst), intent(OUT) :: ierlst !< ierlst(k) contains the error flag corresponding with rslst(k)
  integer, intent(OUT) :: lst !< Number of subintervals needed for the integration
  real(dp), dimension(limit), intent(OUT) :: alist !< Vector for internal use
  real(dp), dimension(limit), intent(OUT) :: blist !< Vector for internal use
  NUMFOR_KINDR, dimension(limit), intent(OUT) :: rlist !< Vector for internal use
  real(dp), dimension(limit), intent(OUT) :: elist !< Vector for internal use
  integer, dimension(limit), intent(OUT) :: iord !< Used for subdivision
  integer, dimension(limit), intent(OUT) :: nnlog !< Used for subdivision
  real(dp), dimension(maxp1, 25), intent(OUT) :: chebmo !< Space to storage Chebyshev moments
!! Examples:
!!

  real(dp) :: abseps
  real(dp) :: correc
  real(dp) :: n_cycle
  real(dp) :: c1
  real(dp) :: c2
  real(dp) :: dl
  real(dp) :: drl
  real(dp) :: ep
  real(dp) :: eps
  real(dp) :: epsa
  real(dp) :: errsum
  real(dp) :: fact
  integer ktmin
  integer l
  integer ll
  integer momcom
  integer nev
  integer nres
  integer numrl2
  real(dp), parameter :: p = 0.9_8
  real(dp), parameter :: pi = M_PI
  real(dp) :: p1
  NUMFOR_KINDR :: psum(52)
  NUMFOR_KINDR :: reseps
  NUMFOR_KINDR :: res3la(3)
!
!  The dimension of  psum  is determined by the value of
!  limexp in QEXTR (psum must be
!  of dimension (limexp+2) at least).
!
!  Test on validity of parameters.
!
  result = 0.0_8
  abserr = 0.0_8
  neval = 0
  lst = 0
  ier = 0
  if ((flgw /= 1 .and. flgw /= 2) .or. &
      epsabs <= 0.0_8 .or. limlst < 3) then
    ier = 6
    return
  end if
  if (omega == 0.0_8) then
    if (flgw == 1) then
      call qags(f, 0._8, nf_inf, result, epsabs, 0.0_8, 'qk21', abserr, neval, ier)
      ! else
      !   result = 0.0_8
      !   abserr = 0.0_8
      !   neval = 0
      !   ier = 0
    end if
    rslst(1) = result
    erlst(1) = abserr
    ierlst(1) = ier
    lst = 1
    return
  end if
!
!  Initializations.
!
  l = int(abs(omega))
  dl = 2 * l + 1
  n_cycle = dl * pi / abs(omega)
  ier = 0
  ktmin = 0
  neval = 0
  numrl2 = 0
  nres = 0
  c1 = a
  c2 = n_cycle + a
  p1 = 1.0_8 - p
  eps = epsabs
  if (epsabs > tiny(epsabs) / p1) then
    eps = epsabs * p1
  end if
  ep = eps
  fact = 1.0_8
  correc = 0.0_8
  abserr = 0.0_8
  errsum = 0.0_8
  do lst = 1, limlst
    !
    !  Integrate over current subinterval.
    !
    !   dla = lst
    epsa = eps * fact
    call qfour(f, c1, c2, omega, flgw, epsa, 0.0d+00, limit, lst, maxp1, rslst(lst), erlst(lst), nev, &
      & ierlst(lst), alist, blist, rlist, elist, iord, nnlog, momcom, chebmo)
    neval = neval + nev
    fact = fact * p
    errsum = errsum + erlst(lst)
    drl = 50._8 * abs(rslst(lst))
    !
    !  Test on accuracy with partial sum.
    !
    if ((errsum + drl) <= epsabs .and. lst >= 6) then
      go to 80
    end if
    correc = max(correc, erlst(lst))
    if (ierlst(lst) /= 0) then
      eps = max(ep, correc * p1)
      ier = 7
    end if
    if (ier == 7 .and. (errsum + drl) <= correc * 10._8 .and. lst > 5) go to 80
    numrl2 = numrl2 + 1
    if (lst <= 1) then
      psum(1) = rslst(1)
      go to 40
    end if
    psum(numrl2) = psum(ll) + rslst(lst)
    IF (lst == 2) go to 40

    !
    !  Test on maximum number of subintervals
    !
    if (lst == limlst) then
      ier = 8
    end if
    !
    !  Perform new extrapolation
    !
    call qextr(numrl2, psum, reseps, abseps, res3la, nres)
    !
    !  Test whether extrapolated result is influenced by roundoff
    !
    ktmin = ktmin + 1
    if (ktmin >= 15 .and. abserr <= 1.0e-03 * (errsum + drl)) then
      ier = 9
    end if
    if (abseps <= abserr .or. lst == 3) then
      abserr = abseps
      result = reseps
      ktmin = 0
      !
      !  If IER is not 0, check whether direct result (partial sum) or
      !  extrapolated result yields the best integral approximation
      !
      IF ((abserr + 10 * correc) <= epsabs) exit
      IF (abserr <= epsabs .and. 10 * correc >= epsabs) exit

    end if
    IF (ier /= 0 .and. ier /= 7) exit

40  continue
    ll = numrl2
    c1 = c2
    c2 = c2 + n_cycle
  end do
!
!  Set final result and error estimate.
!
  abserr = abserr + 10 * correc
  IF (ier == 0) return

  IF (abs(result) /= 0.0_8 .and. abs(psum(numrl2)) /= 0.0_8) go to 70
  IF (abserr > errsum) go to 80
  IF (abs(psum(numrl2)) == 0.0_8) return

70 continue
  if (abserr / abs(result) <= (errsum + drl) / abs(psum(numrl2))) then
    IF (ier >= 1 .and. ier /= 7) abserr = abserr + drl
    return
  end if
80 continue
  result = psum(numrl2)
  abserr = errsum + drl
  end PRNM(qawfe)
