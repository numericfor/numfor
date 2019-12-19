!> qng  estimates an integral using non-adaptive integration.
!!
!! The routine calculates an approximation INTVAL to a definite integral
!!   I = integrate(f, x, a, b)
!! hopefully satisfying
!!   || I - INTVAL || <= max ( EPSABS, EPSREL * ||I|| ).
!!
!! The routine is a simple non-adaptive automatic integrator, based on a sequence
!! of rules with increasing degree of algebraic precision (Patterson, 1968).
#ifdef NF_FARG
  PRNM(qng) (ff, a, b, args, IntVal, epsabs, epsrel, abserr, neval, ier)
  implicit none
  NUMFOR_KINDF :: ff         !< function to integrate
  real(dp), dimension(:), intent(IN) :: args
  NF_FARG :: f
#else
  PRNM(qng) (f, a, b, IntVal, epsabs, epsrel, abserr, neval, ier)
  implicit none
  NUMFOR_KINDF :: f         !< function to integrate
#endif

  NUMFOR_KINDR, intent(OUT) :: IntVal
  real(dp), intent(IN) :: a
  real(dp), intent(IN) :: b
  real(dp), optional, intent(IN) :: epsabs
  real(dp), optional, intent(IN) :: epsrel
  real(dp), optional, intent(OUT) :: abserr
  integer, optional, intent(OUT) :: neval
  integer, optional, intent(OUT) :: ier
!*****************************************************************************80
!    Output, integer
!           ier    - ier = 0 normal and reliable termination of the
!                            routine. it is assumed that the requested
!                            accuracy has been achieved.
!                    ier > 0 abnormal termination of the routine. it is
!                            assumed that the requested accuracy has
!                            not been achieved.
!                    ier = 1 the maximum number of steps has been
!                            executed. the integral is probably too
!                            difficult to be calculated by qng.
!                        = 6 the input is invalid, because
!                            epsabs < 0 and epsrel < 0,
!                            IntVal, abserr and neval are set to zero.
!
!  Local Parameters:
!
!           centr  - mid point of the integration interval
!           hlgth  - half-length of the integration interval
!           fcentr - function value at mid point
!           absc   - abscissa
!           fval   - function value
!           savfun - array of function values which have already
!                    been computed
!           res10  - 10-point Gauss IntVal
!           res21  - 21-point Kronrod IntVal
!           res43  - 43-point IntVal
!           res87  - 87-point IntVal
!           resabs - approximation to the integral of abs(f)
!           resasc - approximation to the integral of abs(f-i/(b-a))
!
  real(dp) :: absc
  real(dp) :: centr
  real(dp) :: dhlgth
  NUMFOR_KINDR :: fcentr
  NUMFOR_KINDR :: fval
  NUMFOR_KINDR :: fval1
  NUMFOR_KINDR :: fval2
  NUMFOR_KINDR :: fv1(5)
  NUMFOR_KINDR :: fv2(5)
  NUMFOR_KINDR :: fv3(5)
  NUMFOR_KINDR :: fv4(5)
  real(dp) :: hlgth
  integer ipx
  integer k
  integer l

  NUMFOR_KINDR :: res10
  NUMFOR_KINDR :: res21
  NUMFOR_KINDR :: res43
  NUMFOR_KINDR :: res87
  real(dp) :: resabs
  real(dp) :: resasc
  NUMFOR_KINDR :: reskh
  NUMFOR_KINDR :: savfun(21)

  real(dp):: epsabs_
  real(dp):: epsrel_
  real(dp) :: abserr_

!    The following data contain the abscissae and weights of the integration rules used.
!
!    x1      abscissae common to the 10-, 21-, 43- and 87-point
!            rule
!    x2      abscissae common to the 21-, 43- and 87-point rule
!    x3      abscissae common to the 43- and 87-point rule
!    x4      abscissae of the 87-point rule
!    w10     weights of the 10-point formula
!    w21a    weights of the 21-point formula for abscissae x1
!    w21b    weights of the 21-point formula for abscissae x2
!    w43a    weights of the 43-point formula for absissae x1, x3
!    w43b    weights of the 43-point formula for abscissae x3
!    w87a    weights of the 87-point formula for abscissae x1,
!            x2 and x3
!    w87b    weights of the 87-point formula for abscissae x4
!
  real(8), dimension(5), parameter :: x1 = [0.9739065285171717_8, 0.8650633666889845_8, &
    & 0.6794095682990244_8, 0.4333953941292472_8, 0.1488743389816312_8]

  real(8), dimension(5), parameter :: x2 = [0.9956571630258081_8, 0.9301574913557082_8, &
    & 0.7808177265864169_8, 0.5627571346686047_8, 0.2943928627014602_8]

  real(8), dimension(*), parameter :: x3 = [0.9993333609019321_8, .9874334029080889_8,&
    & 0.9548079348142663_8, 0.9001486957483283_8, 0.8251983149831142_8, .7321483889893050_8,&
    & .6228479705377252_8, .4994795740710565_8, .3649016613465808_8, .2222549197766013_8&
    &, 0.07465061746138332_8]

  real(8), dimension(22), parameter :: x4 = [0.9999029772627292_8, 0.9979898959866787_8,&
    & 0.9921754978606872_8, .9813581635727128_8, .9650576238583846_8, .9431676131336706_8,&
    & 0.9158064146855072_8, .8832216577713165_8, .8457107484624157_8, .8035576580352310_8,&
    & 0.7570057306854956_8, .7062732097873218_8, .6515894665011779_8, .5932233740579611_8,&
    & 0.5314936059708319_8, .4667636230420228_8, .3994248478592188_8, .3298748771061883_8,&
    & 0.2585035592021616_8, .1856953965683467_8, .1118422131799075_8, .03735212339461987_8]

  real(dp), dimension(5), parameter :: w10 = [6.667134430868814e-02_8, 1.494513491505806e-01_8,&
    & 2.190863625159820e-01_8, 2.692667193099964e-01_8, 2.955242247147529e-01_8]

  real(dp), dimension(5), parameter :: w21a = [3.255816230796473e-02_8, 7.503967481091995e-02_8,&
    & 1.093871588022976e-01_8, 1.347092173114733e-01_8, 1.477391049013385e-01_8]

  real(dp), dimension(6), parameter :: w21b = [1.169463886737187e-02_8, 5.475589657435200e-02_8,&
    & 9.312545458369761e-02_8, 1.234919762620659e-01_8, 1.427759385770601e-01_8, 1.494455540029169e-01_8]

  real(dp), dimension(10), parameter :: w43a = [1.629673428966656e-02, 3.752287612086950e-02,&
    & 5.469490205825544e-02, 6.735541460947809e-02, 7.387019963239395e-02, 5.768556059769796e-03,&
    & 2.737189059324884e-02, 4.656082691042883e-02, 6.174499520144256e-02, 7.138726726869340e-02]

  real(dp), dimension(12), parameter :: w43b = [1.844477640212414e-03_8, 1.079868958589165e-02_8,&
    & 2.189536386779543e-02_8, 3.259746397534569e-02_8, 4.216313793519181e-02_8, 5.074193960018458e-02_8,&
    & 5.837939554261925e-02_8, 6.474640495144589e-02_8, 6.956619791235648e-02_8, 7.282444147183321e-02_8,&
    & 7.450775101417512e-02_8, 7.472214751740301e-02_8]

  real(dp), dimension(21), parameter :: w87a = [8.148377384149173e-03_8, 1.876143820156282e-02_8,&
    & 2.734745105005229e-02_8, 3.367770731163793e-02_8, 3.693509982042791e-02_8, 2.884872430211531e-03_8,&
    & 1.368594602271270e-02_8, 2.328041350288831e-02_8, 3.087249761171336e-02_8, 3.569363363941877e-02_8,&
    & 9.152833452022414e-04_8, 5.399280219300471e-03_8, 1.094767960111893e-02_8, 1.629873169678734e-02_8,&
    & 2.108156888920384e-02_8, 2.537096976925383e-02_8, 2.918969775647575e-02_8, 3.237320246720279e-02_8,&
    & 3.478309895036514e-02_8, 3.641222073135179e-02_8, 3.725387550304771e-02_8]

  real(dp), dimension(23), parameter :: w87b = [2.741455637620724e-04_8, 1.807124155057943e-03_8,&
    & 4.096869282759165e-03_8, 6.758290051847379e-03_8, 9.549957672201647e-03_8, 1.232944765224485e-02_8,&
    & 1.501044734638895e-02_8, 1.754896798624319e-02_8, 1.993803778644089e-02_8, 2.219493596101229e-02_8,&
    & 2.433914712600081e-02_8, 2.637450541483921e-02_8, 2.828691078877120e-02_8, 3.005258112809270e-02_8,&
    & 3.164675137143993e-02_8, 3.305041341997850e-02_8, 3.425509970422606e-02_8, 3.526241266015668e-02_8,&
    & 3.607698962288870e-02_8, 3.669860449845609e-02_8, 3.712054926983258e-02_8, 3.733422875193504e-02_8,&
    & 3.736107376267902e-02_8]

#ifdef NF_FARG
  f%f => ff
  allocate (f%args(size(args)))
  f%args = args
#endif
!
! Optional input arguments
  epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
  epsrel_ = DEFAULT_EPSREL; IF (Present(epsrel)) epsrel_ = epsrel
!
!  Test on validity of parameters.
!
  IntVal = 0.0_8
  abserr_ = 0.0_8
  if (epsabs_ < 0.0_8 .and. epsrel_ < 0.0_8) then
    IF (Present(ier)) ier = 6
    IF (Present(neval)) neval = 0
    return
  end if

  hlgth = 0.5_8 * (b - a)
  dhlgth = abs(hlgth)
  centr = 0.5_8 * (b + a)
  fcentr = NUMFOR_EVAL_F(f, centr)
  IF (Present(neval)) neval = 21
  IF (Present(ier)) ier = 1
!
!  Compute the integral using the 10- and 21-point formula.
!
  do l = 1, 3
    if (l == 1) then
      res10 = 0.0_8
      res21 = w21b(6) * fcentr
      resabs = w21b(6) * abs(fcentr)
      do k = 1, 5
        absc = hlgth * x1(k)
        fval1 = NUMFOR_EVAL_F(f, centr + absc)
        fval2 = NUMFOR_EVAL_F(f, centr - absc)
        fval = fval1 + fval2
        res10 = res10 + w10(k) * fval
        res21 = res21 + w21a(k) * fval
        resabs = resabs + w21a(k) * (abs(fval1) + abs(fval2))
        savfun(k) = fval
        fv1(k) = fval1
        fv2(k) = fval2
      end do
      ipx = 5
      do k = 1, 5
        ipx = ipx + 1
        absc = hlgth * x2(k)
        fval1 = NUMFOR_EVAL_F(f, centr + absc)
        fval2 = NUMFOR_EVAL_F(f, centr - absc)
        fval = fval1 + fval2
        res21 = res21 + w21b(k) * fval
        resabs = resabs + w21b(k) * (abs(fval1) + abs(fval2))
        savfun(ipx) = fval
        fv3(k) = fval1
        fv4(k) = fval2
      end do
      !
      !  Test for convergence.
      !
      IntVal = res21 * hlgth
      resabs = resabs * dhlgth
      reskh = 0.5_8 * res21
      resasc = w21b(6) * abs(fcentr - reskh)
      do k = 1, 5
        resasc = resasc + w21a(k) * (abs(fv1(k) - reskh) + abs(fv2(k) - reskh)) &
                 + w21b(k) * (abs(fv3(k) - reskh) + abs(fv4(k) - reskh))
      end do
      abserr_ = abs((res21 - res10) * hlgth)
      resasc = resasc * dhlgth
    else if (l == 2) then      !  Compute the integral using the 43-point formula.
      res43 = w43b(12) * fcentr
      IF (Present(neval)) neval = 43
      res43 = res43 + sum(savfun(1:10) * w43a(1:10))
      do k = 1, 11
        ipx = ipx + 1
        absc = hlgth * x3(k)
        fval = NUMFOR_EVAL_F(f, absc + centr) + NUMFOR_EVAL_F(f, centr - absc)
        res43 = res43 + fval * w43b(k)
        savfun(ipx) = fval
      end do
      !
      !  Test for convergence.
      IntVal = res43 * hlgth
      abserr_ = abs((res43 - res21) * hlgth)
      !
    else if (l == 3) then      !  Compute the integral using the 87-point formula.
      res87 = w87b(23) * fcentr
      IF (Present(neval)) neval = 87
      res87 = res87 + sum(savfun(1:21) * w87a(1:21))
      do k = 1, 22
        absc = hlgth * x4(k)
        res87 = res87 + w87b(k) * (NUMFOR_EVAL_F(f, absc + centr) + NUMFOR_EVAL_F(f, centr - absc))
      end do
      IntVal = res87 * hlgth
      abserr_ = abs((res87 - res43) * hlgth)
    end if

    IF (resasc /= 0._8 .and. abserr_ /= 0._8) abserr_ = resasc * min(1._8, (200 * abserr_ / resasc)**1.5_8)
    IF (resabs > tiny(resabs) / (50 * epsilon(resabs))) abserr_ = max((50 * epsilon(resabs)) * resabs, abserr_)
    if (abserr_ <= max(epsabs_, epsrel_ * abs(IntVal))) then
      IF (Present(ier)) ier = 0
      IF (Present(abserr)) abserr = abserr_
      return
    end if
  end do
  IF (Present(abserr)) abserr = abserr_

  end PRNM(qng)

! ========================================================================
!> QAG approximates an integral over a finite interval.
!!
!! The routine calculates an approximation RESULT to a definite integral
!!   I = integrate(f(x,...), x, a, b)
!! hopefully satisfying
!!   || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!!
!! QAG is a simple globally adaptive integrator using the strategy of
!! Aind (Piessens, 1973).  It is possible to choose between 6 pairs of
!! Gauss-Kronrod quadrature formulae for the rule evaluation component.
!! The pairs of high degree of precision are suitable for handling
!! integration difficulties due to a strongly oscillating integrand.
!!
#ifdef NF_FARG
  PRNM(qag) (ff, a, b, args, IntVal, epsabs, epsrel, gkrule, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: ff         !< function to integrate
  real(dp), dimension(:), intent(IN) :: args
  NF_FARG :: f
#else
  PRNM(qag) (f, a, b, IntVal, epsabs, epsrel, gkrule, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: f         !< function to integrate
#endif
  NUMFOR_KINDR, intent(OUT) :: IntVal !< Estimation value for the integral
  real(dp), intent(IN) :: a       !< lower limit
  real(dp), intent(IN) :: b       !< upper limit
  real(dp), optional, intent(IN) :: epsabs  !< Desired absolute error
  real(dp), optional, intent(IN) :: epsrel  !< Desired relative error
  character(len=4), optional, intent(IN) :: gkrule      !< rule to use
  real(dp), optional, intent(OUT) :: abserr !< Estimation of error
  integer, optional, intent(OUT) :: neval   !< Number of function evaluations
  integer, optional, intent(OUT) :: ier     !< Info (error) code
  NUMFOR_QP_INFO, optional, intent(INOUT) :: info !< Workspace array and info
  !
!*****************************************************************************80
!    Input, integer KEY, chooses the order of the local integration rule:
!    1,  7 Gauss points, 15 Gauss-Kronrod points,
!    2, 10 Gauss points, 21 Gauss-Kronrod points,
!    3, 15 Gauss points, 31 Gauss-Kronrod points,
!    4, 20 Gauss points, 41 Gauss-Kronrod points,
!    5, 25 Gauss points, 51 Gauss-Kronrod points,
!    6, 30 Gauss points, 61 Gauss-Kronrod points.
!
!    Output, real(8) :: INTVAL, the estimated value of the integral.
!
!    Output, real(8) :: ABSERR, an estimate of || I - INTVAL ||.
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
!  Local parameters:
!
!    LIMIT is the maximum number of subintervals allowed in
!    the subdivision process of QAGE.

  NUMFOR_QP_INFO :: inf_
  real(dp) :: epsabs_
  real(dp) :: epsrel_
  ! integer:: key_
  character(len=4) :: rule
  real(dp) :: abserr_
  integer :: neval_
  integer :: ier_
#ifdef NF_FARG
  f%f => ff
  allocate (f%args(size(args)))
  f%args = args
#endif

  epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
  epsrel_ = DEFAULT_EPSREL; IF (Present(epsrel)) epsrel_ = epsrel
  rule = 'qk15'; IF (Present(gkrule)) rule = gkrule
  ! key_ = 1; IF (Present(key)) key_ = max(min(key, 6), 1)
  if (Present(info) .and. info%size /= 0) then
    inf_ = info
  else
    inf_ = 500
  end if
  call qag(f, a, b, epsabs_, epsrel_, rule, IntVal, abserr_, neval_, ier_, &
    & inf_%size, inf_%alist, inf_%blist, inf_%rlist, inf_%elist, inf_%iord, inf_%last)
  IF (Present(neval)) neval = neval_
  IF (Present(ier)) ier = ier_
  IF (Present(abserr)) abserr = abserr_
  if (Present(info)) then
    info = inf_
    info%msg = "Routine qag:  "//qp_err_msg(ier_)
  end if
  end PRNM(qag)

! ! ========================================================================
! !> qags estimates the integral of a function.
! !!
! !! The routine calculates an approximation RESULT to a definite integral
! !!   I = integrate(f(x, ...), x, a, b)
! !! hopefully satisfying
! !!   || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
! #ifdef NF_FARG
!   PRNM(qags) (ff, a, b, args, IntVal, epsabs, epsrel, rule, abserr, neval, ier, info)
!   implicit none
!   NUMFOR_KINDF :: ff         !< function to integrate
!   real(dp), dimension(:), intent(IN) :: args
!   NF_FARG :: f
! #else
!   PRNM(qags) (f, a, b, IntVal, epsabs, epsrel, rule, abserr, neval, ier, info)
!   implicit none
!   NUMFOR_KINDF :: f         !< function to integrate
! #endif
!   character(len=4), optional, intent(IN) :: rule      !< rule to use
!   real(dp), intent(IN) :: a       !< lower limit
!   real(dp), intent(IN) :: b       !< upper limit
!   real(dp), optional, intent(IN) :: epsabs  !< Desired absolute error
!   real(dp), optional, intent(IN) :: epsrel  !< Desired relative error
!   NUMFOR_KINDR, intent(OUT) :: IntVal !< Estimation of integral
!   real(dp), optional, intent(OUT) :: abserr !< Estimation of error
!   integer, optional, intent(OUT) :: neval   !< Number of function evaluations
!   integer, optional, intent(OUT) :: ier     !< Info (error) code
!   NUMFOR_QP_INFO, optional, intent(INOUT) :: info !< Workspace array and info
! !! Examples:
! !!

!   NUMFOR_QP_INFO :: inf_
!   real(dp) :: epsabs_
!   real(dp) :: epsrel_
!   character(len=4):: rule_
!   real(dp) :: abserr_
!   integer :: neval_
!   integer :: ier_

! #ifdef NF_FARG
!   f%f => ff
!   allocate (f%args(size(args)))
!   f%args = args
! #endif

!   epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
!   epsrel_ = DEFAULT_EPSREL; IF (Present(epsrel)) epsrel_ = epsrel

!   rule_ = 'qk21'; IF (Present(rule)) rule_ = rule

!   if (Present(info)) then
!     inf_ = info
!   else
!     inf_ = 500
!   end if

!   call qags(f, a, b, epsabs_, epsrel_, rule_, IntVal, abserr_, neval_, ier_, &
!       & inf_%size, inf_%alist, inf_%blist, inf_%rlist, inf_%elist, inf_%iord, inf_%last)
!   IF (Present(neval)) neval = neval_
!   IF (Present(ier)) ier = ier_
!   IF (Present(abserr)) abserr = abserr_
!   if (Present(info)) then
!     info = inf_
!     info%msg = "Routine qags: "//qp_err_msg(ier_)
!   end if
!   end PRNM(qags)

! ========================================================================
#ifdef NF_FARG
  PRNM(qags) (ff, a, b, args, IntVal, epsabs, epsrel, gkrule, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: ff         !< function to integrate
  real(dp), dimension(:), intent(IN) :: args
  NF_FARG :: f
#else
  PRNM(qags) (f, a, b, IntVal, epsabs, epsrel, gkrule, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: f         !< function to integrate
#endif
  real(dp), intent(IN) :: a       !< lower limit
  real(dp), intent(IN) :: b       !< upper limit
  NUMFOR_KINDR, intent(OUT) :: IntVal !< Estimation of integral

  real(dp), optional, intent(IN) :: epsabs  !< Desired absolute error
  real(dp), optional, intent(IN) :: epsrel  !< Desired relative error
  character(len=4), optional, intent(IN) :: gkrule      !< rule to use
  real(dp), optional, intent(OUT) :: abserr !< Estimation of error
  integer, optional, intent(OUT) :: neval   !< Number of function evaluations
  integer, optional, intent(OUT) :: ier     !< Error code
  NUMFOR_QP_INFO, optional, intent(INOUT) :: info !< extra parameters
!! Examples:
!!
  NUMFOR_QP_INFO :: inf_
  real(dp) :: epsabs_
  real(dp) :: epsrel_
  character(len=4):: rule
  real(dp) :: abserr_
  integer :: neval_
  integer :: ier_
#ifdef NF_FARG
  f%f => ff
  allocate (f%args(size(args)))
  f%args = args
#endif

  epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
  epsrel_ = DEFAULT_EPSREL; IF (Present(epsrel)) epsrel_ = epsrel
  rule = 'qk15'; IF (Present(gkrule)) rule = gkrule
  if (Present(info) .and. info%size /= 0) then
    inf_ = info
  else
    inf_ = 500
  end if
  if (is_inf(a) /= 0 .and. is_inf(b) /= 0) then ! x ∈ (-∞,+∞)
    call qags(fb_transform, 0._8, 1._8, epsabs_, epsrel_, rule, IntVal, abserr_, neval_, ier_, &
      & inf_%size, inf_%alist, inf_%blist, inf_%rlist, inf_%elist, inf_%iord, inf_%last)
  else if (is_inf(a) == -1) then ! x ∈ (-∞,b)
    call qags(fib_transform, 0._8, 1._8, epsabs_, epsrel_, rule, IntVal, abserr_, neval_, ier_, &
      & inf_%size, inf_%alist, inf_%blist, inf_%rlist, inf_%elist, inf_%iord, inf_%last)
  else if (is_inf(b) == 1) then ! x ∈ (a,+∞)
    call qags(fai_transform, 0._8, 1._8, epsabs_, epsrel_, rule, IntVal, abserr_, neval_, ier_, &
      & inf_%size, inf_%alist, inf_%blist, inf_%rlist, inf_%elist, inf_%iord, inf_%last)
  else                          ! a and b finite
    IF (.not. Present(gkrule)) rule = 'qk21'
    call qags(f, a, b, epsabs_, epsrel_, rule, IntVal, abserr_, neval_, ier_, &
        & inf_%size, inf_%alist, inf_%blist, inf_%rlist, inf_%elist, inf_%iord, inf_%last)
  end if

  IF (Present(neval)) neval = neval_
  IF (Present(ier)) ier = ier_
  IF (Present(abserr)) abserr = abserr_
  if (Present(info)) then
    info = inf_
    info%msg = "Routine qags:  "//qp_err_msg(ier_)
  end if

contains

  ! integrate(f(x),-Inf,Inf) = integrate((f((1-t)/t) + f(-(1-t)/t))/t^2,0,1)
  function fb_transform(t) result(y)
    implicit none
    NUMFOR_KINDR :: y
    real(dp), intent(IN) :: t
    real(dp) :: x
    x = (1 - t) / t
    y = NUMFOR_EVAL_F(f, x) + NUMFOR_EVAL_F(f, -x)     ! f is taken from parent routine
    y = (y / (t * t))
  end function fb_transform

  ! integrate(f(x),a,Inf) = integrate(f(a+(1-t)/t)/t^2,0,1)
  function fai_transform(t) result(y)
    implicit none
    NUMFOR_KINDR :: y !<
    real(dp), intent(IN) :: t
    y = NUMFOR_EVAL_F(f, a + (1 - t) / t)    ! f, a are taken from parent routine
    y = (y / (t * t))
  end function fai_transform
  ! integrate(f(x),-Inf,b) = integrate(f(b-(1-t)/t)/t^2,0,1)
  function fib_transform(t) result(y)
    implicit none
    NUMFOR_KINDR :: y !<
    real(dp), intent(IN) :: t
    y = NUMFOR_EVAL_F(f, b - (1 - t) / t)    ! f, b are taken from parent routine
    y = (y / (t * t))
  end function fib_transform

  end PRNM(qags)

! ========================================================================
!! QAGP computes a definite integral with problematic points.
!!
!! The routine calculates an approximation INTVAL to a definite integral
!!   I = integrate(f, x, a, b)
!! hopefully satisfying
!!   || I - INTVAL || <= max ( EPSABS, EPSREL * ||I|| ).
!!
!! Interior break points of the integration interval, where local difficulties of the integrand may
!! occur, such as singularities or discontinuities, are provided by the user.
!!
#ifdef NF_FARG
  PRNM(qagp) (ff, a, b, points, args, IntVal, epsabs, epsrel, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: ff         !< function to integrate
  real(dp), dimension(:), intent(IN) :: args
  NF_FARG :: f
#else
  PRNM(qagp) (f, a, b, points, IntVal, epsabs, epsrel, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: f         !< function to integrate
#endif
  NUMFOR_KINDR, intent(OUT) :: IntVal !< Estimation value for the integral
  real(dp), intent(IN) :: a !<
  real(dp), intent(IN) :: b !<
  real(dp), dimension(:), intent(IN)  :: points

  real(dp), optional, intent(IN) :: epsabs !<
  real(dp), optional, intent(IN) :: epsrel !<
  real(dp), optional, intent(OUT) :: abserr !<
  integer, optional, intent(OUT) :: neval !<
  integer, optional, intent(OUT) :: ier !<
  NUMFOR_QP_INFO, optional, intent(INOUT) :: info
  NUMFOR_QP_INFO :: inf_
  real(dp) :: epsabs_
  real(dp) :: epsrel_
  real(dp) :: abserr_
  integer :: neval_
  integer :: ier_

#ifdef NF_FARG
  f%f => ff
  allocate (f%args(size(args)))
  f%args = args
#endif

  epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
  epsrel_ = DEFAULT_EPSREL; IF (Present(epsrel)) epsrel_ = epsrel
  if (Present(info) .and. info%size /= 0) then
    inf_ = info
  else
    inf_ = 500
  end if

  call qagp(f, a, b, size(points) + 2, points, epsabs_, epsrel_, IntVal, abserr_, neval_, ier_, &
        & inf_%size, inf_%alist, inf_%blist, inf_%rlist, inf_%elist, inf_%iord, inf_%last)

  IF (Present(neval)) neval = neval_
  IF (Present(ier)) ier = ier_
  IF (Present(abserr)) abserr = abserr_
  if (Present(info)) then
    info = inf_
    info%msg = "Routine qagp: "//qp_err_msg(ier_)
  end if

  end PRNM(qagp)

! ========================================================================
!> qawo computes the integrals of oscillatory integrands.
!!
!! The routine calculates an approximation INTVAL to a given
!! definite integral
!!   I = Integral ( A <= X <= B ) F(X) * cos ( OMEGA * X ) dx
!! or
!!   I = Integral ( A <= X <= B ) F(X) * sin ( OMEGA * X ) dx
!! hopefully satisfying following claim for accuracy
!!   | I - INTVAL | <= max ( epsabs, epsrel * |I| ).
#ifdef NF_FARG
  PRNM(qawo) (ff, a, b, omega, flgw, args, IntVal, epsabs, epsrel, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: ff         !< function to integrate
  real(dp), dimension(:), intent(IN) :: args
  NF_FARG :: f
#else
  PRNM(qawo) (f, a, b, omega, flgw, IntVal, epsabs, epsrel, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: f         !< function to integrate
#endif
  integer, parameter :: limit = 500
  integer, parameter :: maxp1 = 21
  NUMFOR_KINDR, INTENT(OUT) :: IntVal
  real(dp), intent(IN) :: a !< Lower limit of integration
  real(dp), intent(IN) :: b !< Upper limit of integration
  real(dp), intent(IN) :: omega !< Frequency of oscillation
  integer, intent(IN) :: flgw   !< Flag indicating if oscillation is cosine (flgw=1) or sine (flgw=2)
  real(dp), optional, intent(IN) :: epsabs  !< Desired absolute error
  real(dp), optional, intent(IN) :: epsrel  !< Desired relative error
  real(dp), optional, intent(OUT) :: abserr !< Estimation of error
  integer, optional, intent(OUT) :: neval   !< Number of function evaluations
  integer, optional, intent(OUT) :: ier     !< Output (error) code
  NUMFOR_QP_INFO, optional, intent(INOUT) :: info

  NUMFOR_QP_INFO :: inf_
  real(dp) :: epsabs_
  real(dp) :: epsrel_
  real(dp) :: abserr_
  integer :: neval_
  integer :: ier_
!*****************************************************************************80
!  Parameters:
!
!    Input, external real(8) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(8) :: f
!      real(8) :: x
!    which evaluates the integrand function.
!
!    Input, real(8) :: A, B, the limits of integration.
!
!    Input, real(8) :: OMEGA, the parameter in the weight function.
!
!    Input, integer FLGW, specifies the weight function:
!    1, W(X) = cos ( OMEGA * X )
!    2, W(X) = sin ( OMEGA * X )
!
!    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real(8) :: INTVAL, the estimated value of the integral.
!
!    Output, real(8) :: ABSERR, an estimate of || I - INTVAL ||.
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the
!                             requested accuracy has been achieved.
!                   - ier > 0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                     ier = 1 maximum number of subdivisions allowed
!                             (= leniw/2) has been achieved. one can
!                             allow more subdivisions by increasing the
!                             value of leniw (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement it
!                             is advised to analyze the integrand in
!                             order to determine the integration
!                             difficulties. if the position of a local
!                             difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the
!                             subranges. if possible, an appropriate
!                             special-purpose integrator should
!                             be used which is designed for handling
!                             the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some interior points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved due to
!                             roundoff in the extrapolation table,
!                             and that the returned IntVal is the best
!                             which can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             IntVal, abserr, neval are set to zero.
!
!  Local parameters:
!
!    limit is the maximum number of subintervals allowed in the
!    subdivision process of QFOUR. take care that limit >= 1.
!
!    maxp1 gives an upper bound on the number of Chebyshev moments
!    which can be stored, i.e. for the intervals of lengths
!    abs(b-a)*2**(-l), l = 0, 1, ... , maxp1-2. take care that
!    maxp1 >= 1.

  real(8) :: chebmo(maxp1, 25)
  integer momcom
  integer nnlog(limit)

#ifdef NF_FARG
  f%f => ff
  allocate (f%args(size(args)))
  f%args = args
#endif

  epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
  epsrel_ = DEFAULT_EPSREL; IF (Present(epsrel)) epsrel_ = epsrel
  if (Present(info) .and. info%size /= 0) then
    inf_ = info
  else
    inf_ = limit
  end if

  call qfour(f, a, b, omega, flgw, epsabs_, epsrel_, inf_%size, 1, maxp1, IntVal, abserr_, &
    & neval_, ier_, inf_%alist, inf_%blist, inf_%rlist, inf_%elist, inf_%iord, &
    & nnlog, momcom, chebmo)

  IF (Present(neval)) neval = neval_
  IF (Present(ier)) ier = ier_
  IF (Present(abserr)) abserr = abserr_
  if (Present(info)) then
    info = inf_
    info%msg = "Routine qawo: "//qp_err_msg(ier_)
  end if

  end PRNM(qawo)

! ========================================================================
!> QAWS estimates integrals with algebraico-logarithmic endpoint singularities.
!!
!! This routine calculates an approximation INTVAL to a given
!! definite integral
!!   I = integrate(f(x,...)* w(x), x, a, b)
!! where w shows a singular behavior at the end points (see parameter flgw).
!! Hopefully satisfying following claim for accuracy
!!   abs(i-IntVal) <= max(epsabs,epsrel*abs(i)).
#ifdef NF_FARG
  PRNM(qaws) (ff, a, b, alfa, beta, flgw, args, IntVal, epsabs, epsrel, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: ff         !< function to integrate
  real(dp), dimension(:), intent(IN) :: args
  NF_FARG :: f
#else
  PRNM(qaws) (f, a, b, alfa, beta, flgw, IntVal, epsabs, epsrel, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: f         !< function to integrate
#endif

  NUMFOR_KINDR, INTENT(OUT) :: IntVal
  real(dp), intent(IN) :: a !< Lower limit of integration
  real(dp), intent(IN) :: b !< Upper limit of integration
  real(dp), intent(IN) :: alfa, beta !< exponents of weight function
  integer, intent(IN) :: flgw   !< Flag indicating if oscillation is cosine (flgw=1) or sine (flgw=2)
  real(dp), optional, intent(IN) :: epsabs  !< Desired absolute error
  real(dp), optional, intent(IN) :: epsrel  !< Desired relative error
  real(dp), optional, intent(OUT) :: abserr !< Estimation of error
  integer, optional, intent(OUT) :: neval   !< Number of function evaluations
  integer, optional, intent(OUT) :: ier     !< Output (error) code
  NUMFOR_QP_INFO, optional, intent(INOUT) :: info ! Structure holding working arrays
!*****************************************************************************80
!  Parameters:
!
!    Input, external real(8) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(8) :: f
!      real(8) :: x
!    which evaluates the integrand function.
!
!    Input, real(8) :: A, B, the limits of integration.
!
!    Input, real(8) :: ALFA, BETA, parameters used in the weight function.
!    ALFA and BETA should be greater than -1.
!
!    Input, integer FLGW, indicates which weight function is to be used
!    = 1  (x-a)**alfa*(b-x)**beta
!    = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
!    = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
!    = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
!
!    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real(8) :: INTVAL, the estimated value of the integral.
!
!    Output, real(8) :: ABSERR, an estimate of || I - INTVAL ||.
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
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the data value
!                             of limit in qaws (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement it
!                             is advised to analyze the integrand, in
!                             order to determine the integration
!                             difficulties which prevent the requested
!                             tolerance from being achieved. in case of
!                             a jump discontinuity or a local
!                             singularity of algebraico-logarithmic type
!                             at one or more interior points of the
!                             integration range, one should proceed by
!                             splitting up the interval at these points
!                             and calling the integrator on the
!                             subranges.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behavior occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             b <= a or alfa <= (-1) or beta <= (-1) or
!                             flgw < 1 or flgw > 4 or
!                             epsabs < 0 and epsrel < 0,
!                             IntVal, abserr, neval are set to zero.
!
!  Local parameters:
!
!    LIMIT is the maximum number of subintervals allowed in the
!    subdivision process of qawse. take care that limit >= 2.
!

  integer, parameter :: limit = 500
  NUMFOR_QP_INFO :: inf_
  real(dp) :: epsabs_
  real(dp) :: epsrel_
  real(dp) :: abserr_
  integer :: neval_
  integer :: ier_

#ifdef NF_FARG
  f%f => ff
  allocate (f%args(size(args)))
  f%args = args
#endif

  epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
  epsrel_ = DEFAULT_EPSREL; IF (Present(epsrel)) epsrel_ = epsrel

  if (Present(info) .and. info%size /= 0) then
    inf_ = info
  else
    inf_ = limit
  end if

  call qaws(f, a, b, alfa, beta, flgw, epsabs_, epsrel_, inf_%size, IntVal, &
    & abserr_, neval_, ier_,&
    & inf_%alist, inf_%blist, inf_%rlist, inf_%elist, inf_%iord, inf_%last)

  IF (Present(neval)) neval = neval_
  IF (Present(ier)) ier = ier_
  IF (Present(abserr)) abserr = abserr_
  if (Present(info)) then
    info = inf_
    info%msg = "Routine qaws: "//qp_err_msg(ier_)
  end if
  end PRNM(qaws)

! ========================================================================
!> QAWC computes a Cauchy principal value.
!!
!! The routine calculates an approximation INTVAL to a Cauchy principal
!! value
!!   I = integral of F*W over (A,B),
!! with
!!   W(X) = 1 / (X-C),
!! with C distinct from A and B, hopefully satisfying
!!   || I - INTVAL || <= max ( EPSABS, EPSREL * ||I|| ).
#ifdef NF_FARG
  PRNM(qawc) (ff, a, b, c, args, IntVal, epsabs, epsrel, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: ff         !< function to integrate
  real(dp), dimension(:), intent(IN) :: args
  NF_FARG :: f
#else
  PRNM(qawc) (f, a, b, c, IntVal, epsabs, epsrel, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: f         !< function to integrate
#endif
  NUMFOR_KINDR, INTENT(OUT) :: IntVal
  real(dp), intent(IN) :: a !< Lower limit of integration
  real(dp), intent(IN) :: b !< Upper limit of integration
  real(dp), intent(IN) :: c !< Location of the pole
  real(dp), optional, intent(IN) :: epsabs  !< Desired absolute error
  real(dp), optional, intent(IN) :: epsrel  !< Desired relative error
  real(dp), optional, intent(OUT) :: abserr !< Estimation of error
  integer, optional, intent(OUT) :: neval   !< Number of function evaluations
  integer, optional, intent(OUT) :: ier     !< Output (error) code
  NUMFOR_QP_INFO, optional, intent(INOUT) :: info ! Structure holding working arrays
!*****************************************************************************80
!  Parameters:
!
!    Input, external real(8) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(8) :: f
!      real(8) :: x
!    which evaluates the integrand function.
!
!    Input, real(8) :: A, B, the limits of integration.
!
!    Input, real(8) :: C, a parameter in the weight function, which must
!    not be equal to A or B.
!
!    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real(8) :: INTVAL, the estimated value of the integral.
!
!    Output, real(8) :: ABSERR, an estimate of || I - INTVAL ||.
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
!                             divisions by increasing the data value of
!                             limit in qawc (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement it
!                             is advised to analyze the integrand in
!                             order to determine the integration
!                             difficulties. if the position of a local
!                             difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval one will probably gain from
!                             splitting up the interval at this point
!                             and calling appropriate integrators on the
!                             subranges.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behavior occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             c = a or c = b or
!                             epsabs < 0 and epsrel < 0,
!                             IntVal, abserr, neval are set to zero.
!
!  Local parameters:
!
!    LIMIT is the maximum number of subintervals allowed in the
!    subdivision process of qawce. take care that limit >= 1.
!

  NUMFOR_QP_INFO :: inf_
  real(dp) :: epsabs_
  real(dp) :: epsrel_
  real(dp) :: abserr_
  integer :: neval_
  integer :: ier_

#ifdef NF_FARG
  f%f => ff
  allocate (f%args(size(args)))
  f%args = args
#endif

  epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
  epsrel_ = DEFAULT_EPSREL; IF (Present(epsrel)) epsrel_ = epsrel

  if (Present(info) .and. info%size /= 0) then
    inf_ = info
  else
    inf_ = 500
  end if

  call qawc(f, a, b, c, epsabs_, epsrel_, inf_%size, IntVal, abserr_, neval_, ier_,&
    & inf_%alist, inf_%blist, inf_%rlist, inf_%elist, inf_%iord, inf_%last)

  IF (Present(neval)) neval = neval_
  IF (Present(ier)) ier = ier_
  IF (Present(abserr)) abserr = abserr_
  if (Present(info)) then
    info = inf_
    info%msg = "Routine qawc: "//qp_err_msg(ier_)
  end if
  end PRNM(qawc)

! ========================================================================
!> qawf computes Fourier integrals over the interval [ A, +Infinity ).
!!
!! The routine calculates an approximation INTVAL to a definite integral
!!
!!   I = integral of F*COS(OMEGA*X)
!! or
!!   I = integral of F*SIN(OMEGA*X)
!!
!! over the interval [A,+Infinity), hopefully satisfying
!!
!!   || I - INTVAL || <= EPSABS.
!!
!! If OMEGA = 0 and FLGW = 1, the integral is calculated by means
!! of QAGS, and IER has the meaning as described in the comments of QAGS.
!!
#ifdef NF_FARG
  PRNM(qawf) (ff, a, omega, flgw, args, IntVal, epsabs, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: ff         !< function to integrate
  real(dp), dimension(:), intent(IN) :: args
  NF_FARG :: f
#else
  PRNM(qawf) (f, a, omega, flgw, IntVal, epsabs, abserr, neval, ier, info)
  implicit none
  NUMFOR_KINDF :: f         !< function to integrate
#endif
  NUMFOR_KINDR, INTENT(OUT) :: IntVal
  real(dp), intent(IN) :: a !< Lower limit of integration
  real(dp), intent(IN) :: omega !< Frequency of oscillation
  integer, intent(IN) :: flgw   !< Flag indicating if oscillation is cosine (flgw=1) or sine (flgw=2)
  real(dp), optional, intent(IN) :: epsabs  !< Desired absolute error
  real(dp), optional, intent(OUT) :: abserr !< Estimation of error
  integer, optional, intent(OUT) :: neval   !< Number of function evaluations
  integer, optional, intent(OUT) :: ier     !< Output (error) code
  NUMFOR_QP_INFO, optional, intent(INOUT) :: info ! Structure holding working arrays
!*****************************************************************************80
!  Parameters:
!
!    Input, external real(8) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(8) :: f
!      real(8) :: x
!    which evaluates the integrand function.
!
!    Input, real(8) :: A, the lower limit of integration.
!
!    Input, real(8) :: OMEGA, the parameter in the weight function.
!
!    Input, integer FLGW, indicates which weight functions is used
!    = 1, w(x) = cos(omega*x)
!    = 2, w(x) = sin(omega*x)
!
!    Input, real(8) :: EPSABS, the absolute accuracy requested.
!
!    Output, real(8) :: INTVAL, the estimated value of the integral.
!
!    Output, real(8) :: ABSERR, an estimate of || I - INTVAL ||.
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the
!                             requested accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                    if omega /= 0
!                     ier = 6 the input is invalid because
!                             (flgw /= 1 and flgw /= 2) or
!                              epsabs <= 0
!                              IntVal, abserr, neval, lst are set to
!                              zero.
!                         = 7 abnormal termination of the computation
!                             of one or more subintegrals
!                         = 8 maximum number of cycles allowed
!                             has been achieved, i.e. of subintervals
!                             (a+(k-1)c,a+kc) where
!                             c = (2*int(abs(omega))+1)*pi/abs(omega),
!                             for k = 1, 2, ...
!                         = 9 the extrapolation table constructed for
!                             convergence acceleration of the series
!                             formed by the integral contributions
!                             over the cycles, does not converge to
!                             within the requested accuracy.
!
!  Local parameters:
!
!    Integer LIMLST, gives an upper bound on the number of cycles, LIMLST >= 3.
!    if limlst < 3, the routine will end with ier = 6.
!
!    Integer MAXP1, an upper bound on the number of Chebyshev moments which
!    can be stored, i.e. for the intervals of lengths abs(b-a)*2**(-l),
!    l = 0,1, ..., maxp1-2, maxp1 >= 1.  if maxp1 < 1, the routine will end
!    with ier = 6.

  NUMFOR_QP_INFO :: inf_
  real(dp) :: epsabs_

  real(dp) :: abserr_
  integer :: neval_
  integer :: ier_
  integer, parameter :: limit = 500
  integer, parameter :: limlst = 50
  integer, parameter :: maxp1 = 21
  real(dp) :: chebmo(maxp1, 25)

  real(dp) :: erlst(limlst)
  integer ierlst(limlst)
  integer lst
  integer nnlog(limit)
  NUMFOR_KINDR, dimension(limlst) :: rslst

#ifdef NF_FARG
  f%f => ff
  allocate (f%args(size(args)))
  f%args = args
#endif

! Check input parameters
  epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
  if (Present(info) .and. info%size /= 0) then
    inf_ = info
  else
    inf_ = limit
  end if
  if (limlst < 3 .or. maxp1 < 1) then
    ier_ = 6; neval_ = 0; IntVal = 0.0_8; abserr_ = 0.0_8
  else
    call qawf(f, a, omega, flgw, epsabs_, limlst, inf_%size, maxp1, IntVal, abserr_, neval_, ier_,&
      & rslst, erlst, ierlst, lst, inf_%alist, inf_%blist, inf_%rlist, inf_%elist, inf_%iord, nnlog, chebmo)
  end if

  IF (Present(neval)) neval = neval_
  IF (Present(ier)) ier = ier_
  IF (Present(abserr)) abserr = abserr_
  if (Present(info)) then
    info = inf_
    info%msg = "Routine qawf: "//qp_err_msg(ier_)
  end if

  end PRNM(qawf)

! ! ========================================================================
! !> quads is a wrapper for simple use of a subset of quadpack integrators
! !! with no weight function
! !!
! !! Examples:
! !!
! #ifdef NF_FARG
!   PRNM(quad) (ff, a, b, args, IntVal, abserr, epsabs, epsrel, maxsub, points, neval, ier, info)
!   implicit none
!   NUMFOR_KINDF :: ff         !< function to integrate
!   real(dp), dimension(:), intent(IN) :: args
!   NF_FARG :: f
! #else
!   PRNM(quad) (f, a, b, IntVal, abserr, epsabs, epsrel, maxsub, points, neval, ier, info)
!   implicit none
!   NUMFOR_KINDF :: f         !< function to integrate
! #endif

!   NUMFOR_KINDR, intent(OUT) :: IntVal !<
!   real(dp), intent(IN) :: a !<
!   real(dp), intent(IN) :: b !<
!   real(dp), intent(OUT) :: abserr !<
!   real(dp), optional, intent(IN) :: epsabs !<
!   real(dp), optional, intent(IN) :: epsrel !<
!   integer, optional, intent(IN) :: maxsub !<
!   real(dp), dimension(:), optional, intent(IN) :: points !<
!   integer, optional, intent(OUT) :: ier !<
!   integer, optional, intent(OUT) :: neval !<
!   NUMFOR_QP_INFO, optional, intent(INOUT) :: info ! Structure holding working arrays

!   NUMFOR_QP_INFO :: inf_
!   real(dp) :: epsabs_, epsrel_
!   integer :: maxsub_, ier_, neval_
!   character(len=4) :: gkrule = "qk21"

! #ifdef NF_FARG
!   f%f => ff
!   allocate (f%args(size(args)))
!   f%args = args
! #endif

!   epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
!   epsrel_ = DEFAULT_EPSREL; IF (Present(epsrel)) epsrel_ = epsrel
!   maxsub_ = 500; IF (Present(maxsub)) maxsub_ = maxsub ! "limit" in the original routines
!   if (Present(info)) then
!     inf_ = info
!     IF (inf_%size == 0) inf_ = maxsub_
!   else
!     inf_ = maxsub_
!   end if

! ! If there are any break points, we use them
!   if (Present(points)) then
!     ! Incompatible parameters
!     IF (is_inf(a) /= 0 .or. is_inf(b) /= 0) &
!       & call print_msg("Infinite integration limits cannot be used with break points.", errcode=-1)
!     call qagp(f, a, b, points, IntVal, epsabs_, epsrel_, abserr, neval_, ier_, inf_)
!   else                          ! infinite integration domain
!     call qags(f, a, b, IntVal, epsabs_, epsrel_, gkrule, abserr, neval_, ier_, inf_)
!   end if

!   IF (Present(ier)) ier = ier_
!   IF (Present(neval)) neval = neval_
!   IF (Present(info)) info = inf_

!   end PRNM(quad)
