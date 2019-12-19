!*****************************************************************************80
#ifdef NF_FARG
PRNM(qgk) (ff, a, b, args, IntVal, rule, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qgk) (f, a, b, IntVal, rule, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
real(dp), intent(IN) :: a !<
real(dp), intent(IN) :: b !<
NUMFOR_KINDR, intent(OUT) :: IntVal !<
character(len=4), optional, intent(IN) :: rule      !< rule to use
! integer, optional, intent(IN) :: rule
real(dp), intent(OUT) :: abserr !<
real(dp), intent(OUT) :: resabs !<
real(dp), intent(OUT) :: resasc !<

character(len=4) :: rule_
#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif

rule_ = 'qk21'; IF (Present(rule)) rule_ = rule

select case (rule_)
case ('qk15')
  call qk15(f, a, b, IntVal, abserr, resabs, resasc)
case ('qk31')
  call qk31(f, a, b, IntVal, abserr, resabs, resasc)
case ('qk41')
  call qk41(f, a, b, IntVal, abserr, resabs, resasc)
case ('qk51')
  call qk51(f, a, b, IntVal, abserr, resabs, resasc)
case ('qk61')
  call qk61(f, a, b, IntVal, abserr, resabs, resasc)
case default
  call qk21(f, a, b, IntVal, abserr, resabs, resasc)
end select

end PRNM(qgk)

!! QK15 carries out a 15 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
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
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 15-point Kronrod rule (RESK)
!    obtained by optimal addition of abscissae to the 7-point Gauss rule
!    (RESG).
!
!    Output, real(dp) :: ABSERR, an estimate of | I - RESULT |.
!
!    Output, real(dp) :: RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real(dp) :: RESASC, approximation to the integral | F-I/(B-A) |
!    over [A,B].
!
!  Local Parameters:
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 15-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 7-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point Gauss rule
!
!           wgk    - weights of the 15-point Kronrod rule
!
!           wg     - weights of the 7-point Gauss rule
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 7-point Gauss formula
!           resk   - result of the 15-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
#ifdef NF_FARG
PRNM(qk15) (ff, a, b, args, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qk15) (f, a, b, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
real(dp), intent(IN) :: a !<
real(dp), intent(IN) :: b !<
NUMFOR_KINDR, intent(OUT) :: result !<
real(dp), intent(OUT) :: abserr !<
real(dp), intent(OUT) :: resabs !<
real(dp), intent(OUT) :: resasc !<

real(dp) :: absc
real(dp) :: centr
real(dp) :: dhlgth
real(dp), parameter :: mindiff = 50 * epsilon(resabs)
NUMFOR_KINDR :: fc
NUMFOR_KINDR :: fsum
NUMFOR_KINDR :: fval1
NUMFOR_KINDR :: fval2
NUMFOR_KINDR :: fv1(7)
NUMFOR_KINDR :: fv2(7)
real(dp) :: hlgth
integer j
integer jtw
integer jtwm1
NUMFOR_KINDR :: resg
NUMFOR_KINDR :: resk
NUMFOR_KINDR :: reskh

real(dp), dimension(8), parameter :: xgk = [9.914553711208126e-01_8, 9.491079123427585e-01_8, &
  8.648644233597691e-01_8, 7.415311855993944e-01_8, 5.860872354676911e-01_8, 4.058451513773972e-01_8, &
  & 2.077849550078985e-01_8, 0.0_8]
real(dp), dimension(8), parameter :: wgk = [2.293532201052922e-02_8, 6.309209262997855e-02_8, &
  1.047900103222502e-01_8, 1.406532597155259e-01_8, 1.690047266392679e-01_8, 1.903505780647854e-01_8, &
  & 2.044329400752989e-01_8, 2.094821410847278e-01_8]
real(dp), dimension(4), parameter :: wg = [1.294849661688697e-01_8, 2.797053914892767e-01_8, &
  & 3.818300505051189e-01_8, 4.179591836734694e-01_8]
#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif

!
centr = 0.5_8 * (a + b)
hlgth = 0.5_8 * (b - a)
dhlgth = abs(hlgth)
!
!  Compute the 15-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!
fc = NUMFOR_EVAL_F(f, centr)
resg = fc * wg(4)
resk = fc * wgk(8)
resabs = abs(resk)

do j = 1, 3
  jtw = j * 2
  absc = hlgth * xgk(jtw)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtw) = fval1
  fv2(jtw) = fval2
  fsum = fval1 + fval2
  resg = resg + wg(j) * fsum
  resk = resk + wgk(jtw) * fsum
  resabs = resabs + wgk(jtw) * (abs(fval1) + abs(fval2))
end do

do j = 1, 4
  jtwm1 = j * 2 - 1
  absc = hlgth * xgk(jtwm1)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtwm1) = fval1
  fv2(jtwm1) = fval2
  fsum = fval1 + fval2
  resk = resk + wgk(jtwm1) * fsum
  resabs = resabs + wgk(jtwm1) * (abs(fval1) + abs(fval2))
end do

reskh = resk * 0.5_8
resasc = wgk(8) * abs(fc - reskh)

resasc = resasc + sum(wgk(:7) * (abs(fv1(:7) - reskh) + abs(fv2(:7) - reskh)))

result = resk * hlgth
resabs = resabs * dhlgth
resasc = resasc * dhlgth
abserr = abs((resk - resg) * hlgth)

IF (resasc /= 0.0_8 .and. abserr /= 0.0_8) abserr = resasc * min(1.0_8, (200 * abserr / resasc)**1.5_8)
IF (resabs * mindiff > tiny(resabs)) abserr = max(mindiff * resabs, abserr)

end PRNM(qk15)

!> qk15w applies a 15 point Gauss-Kronrod rule for a weighted integrand.
!!
!! Integral is computed by applying the 15-point Kronrod rule obtained
!! by optimal addition of abscissae to the 7-point Gauss rule
!!
!!  Author:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner
!!
!!  Reference:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
!!    QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
!!
#ifdef NF_FARG
PRNM(qk15w) (ff, w, p, kp, a, b, args, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qk15w) (f, w, p, kp, a, b, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
procedure(qfwght) :: w !< weight function w(x). Needs to be declared external in the calling program
real(dp), dimension(:), intent(IN) :: p !< Parameters for the weight function
integer, intent(IN) :: kp               !< Key for the weight function
real(dp), intent(IN) :: a               !< Lower limit of integration
real(dp), intent(IN) :: b               !< Upper limit of integration
NUMFOR_KINDR, intent(OUT) :: result !< Integral value of f*w over (a,b)
!
real(dp), intent(OUT) :: abserr !< Estimation of error
real(dp), intent(OUT) :: resabs !< integral of abs(f*w) over (a,b)
real(dp), intent(OUT) :: resasc !< integrate(|f-I/(b-a)|, x, a, b)
!! Examples:
!!
!*****************************************************************************80
!  Parameters:
!
!    Input, external real(dp) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(dp) :: f
!      real(dp) :: x
!    which evaluates the integrand function.
!
!              w      - real
!                       function subprogram defining the integrand
!                       weight function w(x). the actual name for w
!                       needs to be declared e x t e r n a l in the
!                       calling program.
!
!    ?, real(dp) :: P1, P2, P3, P4, parameters in the weight function
!
!    Input, integer KP, key for indicating the type of weight function
!
!    Input, real(dp) :: A, B, the limits of integration.
!
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained by
!    optimal addition of abscissae to the 7-point Gauss rule (RESG).
!
!    Output, real(dp) :: ABSERR, an estimate of | I - RESULT |.
!
!    Output, real(dp) :: RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real(dp) :: RESASC, approximation to the integral | F-I/(B-A) |
!    over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc*  - abscissa
!           fval*  - function value
!           resg   - result of the 7-point Gauss formula
!           resk   - result of the 15-point Kronrod formula
!           reskh  - approximation to the mean value of f*w over (a,b),
!                    i.e. to i/(b-a)
!
real(dp), parameter :: mindiff = 50 * epsilon(resabs)
real(dp) :: absc
real(dp) :: absc1
real(dp) :: absc2
real(dp) :: centr
real(dp) :: dhlgth
NUMFOR_KINDR :: fc
NUMFOR_KINDR :: fsum
NUMFOR_KINDR :: fval1
NUMFOR_KINDR :: fval2
NUMFOR_KINDR :: fv1(7)
NUMFOR_KINDR :: fv2(7)
real(dp) :: hlgth
integer j
integer jtw
integer jtwm1
NUMFOR_KINDR :: resg
NUMFOR_KINDR :: resk
NUMFOR_KINDR :: reskh
!
!  the abscissae and weights are given for the interval (-1,1).
!  because of symmetry only the positive abscissae and their
!  corresponding weights are given.
!
!           xgk    - abscissae of the 15-point Gauss-Kronrod rule
!                    xgk(2), xgk(4), ... abscissae of the 7-point Gauss
!                    rule
!                    xgk(1), xgk(3), ... abscissae which are optimally
!                    added to the 7-point Gauss rule
!
!           wgk    - weights of the 15-point Gauss-Kronrod rule
!
!           wg     - weights of the 7-point Gauss rule
!
real(dp), parameter :: wg(4) = [0.1294849661688697_8, 0.2797053914892767_8, 0.3818300505051889_8, &
  & 0.4179591836734694_8]
real(dp), parameter :: wgk(8) = [2.293532201052922e-02_8, 6.309209262997855e-02_8, .1047900103222502_8,&
  & 1.406532597155259e-01_8, 1.690047266392679e-01_8, 1.903505780647854e-01_8, 2.044329400752989e-01_8, &
  & 2.094821410847278e-01_8]
real(dp), parameter :: xgk(8) = [0.9914553711208126_8, 0.9491079123427585_8, .8648644233597691_8, &
  & 0.7415311855993944_8, .5860872354676911_8, .4058451513773972_8, .2077849550789850_8, 0.00_8]

#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif

centr = 0.5_8 * (a + b)
hlgth = 0.5_8 * (b - a)
dhlgth = abs(hlgth)
!
!  Compute the 15-point Kronrod approximation to the integral,
!  and estimate the error.
!
fc = NUMFOR_EVAL_F(f, centr)
fc = fc * w(centr, p, kp)
resg = wg(4) * fc
resk = wgk(8) * fc
resabs = abs(resk)
do j = 1, 3
  jtw = j * 2
  absc = hlgth * xgk(jtw)
  absc1 = centr - absc
  absc2 = centr + absc
  fval1 = NUMFOR_EVAL_F(f, absc1) * w(absc1, p, kp)
  fval2 = NUMFOR_EVAL_F(f, absc2) * w(absc2, p, kp)
  fv1(jtw) = fval1
  fv2(jtw) = fval2
  fsum = fval1 + fval2
  resg = resg + wg(j) * fsum
  resk = resk + wgk(jtw) * fsum
  resabs = resabs + wgk(jtw) * (abs(fval1) + abs(fval2))
end do

do j = 1, 4
  jtwm1 = j * 2 - 1
  absc = hlgth * xgk(jtwm1)
  absc1 = centr - absc
  absc2 = centr + absc
  fval1 = NUMFOR_EVAL_F(f, absc1) * w(absc1, p, kp)
  fval2 = NUMFOR_EVAL_F(f, absc2) * w(absc2, p, kp)
  fv1(jtwm1) = fval1
  fv2(jtwm1) = fval2
  fsum = fval1 + fval2
  resk = resk + wgk(jtwm1) * fsum
  resabs = resabs + wgk(jtwm1) * (abs(fval1) + abs(fval2))
end do

reskh = resk * 0.5_8
resasc = wgk(8) * abs(fc - reskh) + sum(wgk(:7) * (abs(fv1(:7) - reskh) + abs(fv2(:7) - reskh)))

result = resk * hlgth
resabs = resabs * dhlgth
resasc = resasc * dhlgth
abserr = abs((resk - resg) * hlgth)

IF (resasc /= 0._8 .and. abserr /= 0._8) abserr = resasc * min(1.0_8, (200 * abserr / resasc)**1.5_8)
IF (resabs > tiny(resabs) / mindiff) abserr = max(mindiff * resabs, abserr)

end PRNM(qk15w)

!*****************************************************************************80
!
!! QK15I applies a 15 point Gauss-Kronrod quadrature on an infinite interval.
!
!  Discussion:
!
!    The original infinite integration range is mapped onto the interval
!    (0,1) and (a,b) is a part of (0,1).  The routine then computes:
!
!    i = integral of transformed integrand over (a,b),
!    j = integral of abs(transformed integrand) over (a,b).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real(dp) :: F, the name of the function routine, of the form
!      function f ( x )
!      real(dp) :: f
!      real(dp) :: x
!    which evaluates the integrand function.
!
!    Input, real(dp) :: BOUN, the finite bound of the original integration range,
!    or zero if INF is 2.
!
!    Input, integer INF, indicates the type of the interval.
!    -1: the original interval is (-infinity,BOUN),
!    +1, the original interval is (BOUN,+infinity),
!    +2, the original interval is (-infinity,+infinity) and
!    the integral is computed as the sum of two integrals, one
!    over (-infinity,0) and one over (0,+infinity).
!
!    Input, real(dp) :: A, B, the limits of integration, over a subrange of [0,1].
!
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained
!    by optimal addition of abscissae to the 7-point Gauss rule (RESG).
!
!    Output, real(dp) :: ABSERR, an estimate of | I - RESULT |.
!
!    Output, real(dp) :: RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real(dp) :: RESASC, approximation to the integral of the
!    transformated integrand | F-I/(B-A) | over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc*  - abscissa
!           tabsc* - transformed abscissa
!           fval*  - function value
!           resg   - result of the 7-point Gauss formula
!           resk   - result of the 15-point Kronrod formula
!           reskh  - approximation to the mean value of the transformed
!                    integrand over (a,b), i.e. to i/(b-a)
!
#ifdef NF_FARG
PRNM(qk15i) (ff, boun, inf, a, b, args, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qk15i) (f, boun, inf, a, b, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif

real(dp) :: a
real(dp) :: absc
real(dp) :: absc1
real(dp) :: absc2
real(dp) :: abserr
real(dp) :: b
real(dp) :: boun
real(dp) :: centr
real(dp) :: dinf
NUMFOR_KINDR :: fc
NUMFOR_KINDR :: fsum
NUMFOR_KINDR :: fval1
NUMFOR_KINDR :: fval2
NUMFOR_KINDR :: fv1(7)
NUMFOR_KINDR :: fv2(7)
real(dp) :: hlgth
integer inf
integer j
real(dp) :: resabs
real(dp) :: resasc
real(dp), parameter :: mindiff = 50 * epsilon(resabs)
NUMFOR_KINDR :: resg
NUMFOR_KINDR :: resk
NUMFOR_KINDR :: reskh
NUMFOR_KINDR, intent(OUT) :: result
real(dp) :: tabsc1
real(dp) :: tabsc2
real(dp) :: wg(8)
real(dp) :: wgk(8)
real(dp) :: xgk(8)
!
!  the abscissae and weights are supplied for the interval
!  (-1,1).  because of symmetry only the positive abscissae and
!  their corresponding weights are given.
!
!           xgk    - abscissae of the 15-point Kronrod rule
!                    xgk(2), xgk(4), ... abscissae of the 7-point Gauss
!                    rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point Gauss rule
!
!           wgk    - weights of the 15-point Kronrod rule
!
!           wg     - weights of the 7-point Gauss rule, corresponding
!                    to the abscissae xgk(2), xgk(4), ...
!                    wg(1), wg(3), ... are set to zero.
!
data xgk(1), xgk(2), xgk(3), xgk(4), xgk(5), xgk(6), xgk(7), xgk(8)/ &
  9.914553711208126e-01, 9.491079123427585e-01, &
  8.648644233597691e-01, 7.415311855993944e-01, &
  5.860872354676911e-01, 4.058451513773972e-01, &
  2.077849550078985e-01, 0.0000000000000000_8/

data wgk(1), wgk(2), wgk(3), wgk(4), wgk(5), wgk(6), wgk(7), wgk(8)/ &
  2.293532201052922e-02, 6.309209262997855e-02, &
  1.047900103222502e-01, 1.406532597155259e-01, &
  1.690047266392679e-01, 1.903505780647854e-01, &
  2.044329400752989e-01, 2.094821410847278e-01/

data wg(1), wg(2), wg(3), wg(4), wg(5), wg(6), wg(7), wg(8)/ &
  0.0000000000000000_8, 1.294849661688697e-01, &
  0.0000000000000000_8, 2.797053914892767e-01, &
  0.0000000000000000_8, 3.818300505051189e-01, &
  0.0000000000000000_8, 4.179591836734694e-01/

#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif

dinf = min(1, inf)

centr = 0.5_8 * (a + b)
hlgth = 0.5_8 * (b - a)
tabsc1 = boun + dinf * (1.0_8 - centr) / centr
fval1 = NUMFOR_EVAL_F(f, tabsc1)
if (inf == 2) fval1 = fval1 + NUMFOR_EVAL_F(f, -tabsc1)
fc = (fval1 / centr) / centr
!
!  Compute the 15-point Kronrod approximation to the integral,
!  and estimate the error.
!
resg = wg(8) * fc
resk = wgk(8) * fc
resabs = abs(resk)

do j = 1, 7

  absc = hlgth * xgk(j)
  absc1 = centr - absc
  absc2 = centr + absc
  tabsc1 = boun + dinf * (1.0_8 - absc1) / absc1
  tabsc2 = boun + dinf * (1.0_8 - absc2) / absc2
  fval1 = NUMFOR_EVAL_F(f, tabsc1)
  fval2 = NUMFOR_EVAL_F(f, tabsc2)

  if (inf == 2) then
    fval1 = fval1 + NUMFOR_EVAL_F(f, -tabsc1)
    fval2 = fval2 + NUMFOR_EVAL_F(f, -tabsc2)
  end if

  fval1 = (fval1 / absc1) / absc1
  fval2 = (fval2 / absc2) / absc2
  fv1(j) = fval1
  fv2(j) = fval2
  fsum = fval1 + fval2
  resg = resg + wg(j) * fsum
  resk = resk + wgk(j) * fsum
  resabs = resabs + wgk(j) * (abs(fval1) + abs(fval2))
end do

reskh = resk * 0.5_8
resasc = wgk(8) * abs(fc - reskh) + sum(wgk(:7) * (abs(fv1(:7) - reskh) + abs(fv2(:7) - reskh)))

result = resk * hlgth
resasc = resasc * hlgth
resabs = resabs * hlgth
abserr = abs((resk - resg) * hlgth)

IF (resasc /= 0._8 .and. abserr /= 0._8) abserr = resasc * min(1.0_8, (200 * abserr / resasc)**1.5_8)
IF (resabs > tiny(resabs) / mindiff) abserr = max(mindiff * resabs, abserr)

end PRNM(qk15i)

!> qk21
#ifdef NF_FARG
PRNM(qk21) (ff, a, b, args, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qk21) (f, a, b, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
!*****************************************************************************80
!
!! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
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
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 21-point Kronrod rule (resk)
!    obtained by optimal addition of abscissae to the 10-point Gauss
!    rule (resg).
!
!    Output, real(dp) :: ABSERR, an estimate of | I - RESULT |.
!
!    Output, real(dp) :: RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real(dp) :: RESASC, approximation to the integral | F-I/(B-A) |
!    over [A,B].
!
real(dp), intent(IN) :: a !<
real(dp), intent(IN) :: b !<
NUMFOR_KINDR, intent(OUT) :: result !<
real(dp), intent(OUT) :: abserr !<
real(dp), intent(OUT) :: resabs !<
real(dp), intent(OUT) :: resasc !<
real(dp) :: absc
real(dp) :: centr
real(dp) :: dhlgth
NUMFOR_KINDR :: fc
NUMFOR_KINDR :: fsum
NUMFOR_KINDR :: fval1
NUMFOR_KINDR :: fval2
NUMFOR_KINDR :: fv1(10)
NUMFOR_KINDR :: fv2(10)
real(dp) :: hlgth
integer j
integer jtw
integer jtwm1
NUMFOR_KINDR :: resg
NUMFOR_KINDR :: resk
NUMFOR_KINDR :: reskh
real(dp) :: wg(5)
real(dp) :: wgk(11)
real(dp) :: xgk(11)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 21-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 10-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 10-point Gauss rule
!
!           wgk    - weights of the 21-point Kronrod rule
!
!           wg     - weights of the 10-point Gauss rule
!
data xgk(1), xgk(2), xgk(3), xgk(4), xgk(5), xgk(6), xgk(7), xgk(8), &
  xgk(9), xgk(10), xgk(11)/ &
  9.956571630258081e-01, 9.739065285171717e-01, &
  9.301574913557082e-01, 8.650633666889845e-01, &
  7.808177265864169e-01, 6.794095682990244e-01, &
  5.627571346686047e-01, 4.333953941292472e-01, &
  2.943928627014602e-01, 1.488743389816312e-01, &
  0.000000000000000_8/
!
data wgk(1), wgk(2), wgk(3), wgk(4), wgk(5), wgk(6), wgk(7), wgk(8), &
  wgk(9), wgk(10), wgk(11)/ &
  1.169463886737187e-02, 3.255816230796473e-02, &
  5.475589657435200e-02, 7.503967481091995e-02, &
  9.312545458369761e-02, 1.093871588022976e-01, &
  1.234919762620659e-01, 1.347092173114733e-01, &
  1.427759385770601e-01, 1.477391049013385e-01, &
  1.494455540029169e-01/
!
data wg(1), wg(2), wg(3), wg(4), wg(5)/ &
  6.667134430868814e-02, 1.494513491505806e-01, &
  2.190863625159820e-01, 2.692667193099964e-01, &
  2.955242247147529e-01/

#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 10-point Gauss formula
!           resk   - result of the 21-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
centr = 0.5_8 * (a + b)
hlgth = 0.5_8 * (b - a)
dhlgth = abs(hlgth)
!
!  Compute the 21-point Kronrod approximation to the
!  integral, and estimate the absolute error.
!
resg = 0.0_8
fc = NUMFOR_EVAL_F(f, centr)
resk = wgk(11) * fc
resabs = abs(resk)

do j = 1, 5
  jtw = 2 * j
  absc = hlgth * xgk(jtw)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtw) = fval1
  fv2(jtw) = fval2
  fsum = fval1 + fval2
  resg = resg + wg(j) * fsum
  resk = resk + wgk(jtw) * fsum
  resabs = resabs + wgk(jtw) * (abs(fval1) + abs(fval2))
end do

do j = 1, 5
  jtwm1 = 2 * j - 1
  absc = hlgth * xgk(jtwm1)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtwm1) = fval1
  fv2(jtwm1) = fval2
  fsum = fval1 + fval2
  resk = resk + wgk(jtwm1) * fsum
  resabs = resabs + wgk(jtwm1) * (abs(fval1) + abs(fval2))
end do

reskh = resk * 0.5_8
resasc = wgk(11) * abs(fc - reskh)

do j = 1, 10
  resasc = resasc + wgk(j) * (abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
end do

result = resk * hlgth
resabs = resabs * dhlgth
resasc = resasc * dhlgth
abserr = abs((resk - resg) * hlgth)

if (resasc /= 0.0_8 .and. abserr /= 0.0_8) then
  abserr = resasc * min(1.0_8, (2.0e+02 * abserr / resasc)**1.5_8)
end if

if (resabs > tiny(resabs) / (50._8 * epsilon(resabs))) then
  abserr = max((epsilon(resabs) * 50._8) * resabs, abserr)
end if

end PRNM(qk21)

#ifdef NF_FARG
PRNM(qk31) (ff, a, b, args, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qk31) (f, a, b, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
!*****************************************************************************80
!
!! QK31 carries out a 31 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
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
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!                       result is computed by applying the 31-point
!                       Gauss-Kronrod rule (resk), obtained by optimal
!                       addition of abscissae to the 15-point Gauss
!                       rule (resg).
!
!    Output, real(dp) :: ABSERR, an estimate of | I - RESULT |.
!
!    Output, real(dp) :: RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real(dp) :: RESASC, approximation to the integral | F-I/(B-A) |
!    over [A,B].
!
real(dp), intent(IN) :: a !<
real(dp), intent(IN) :: b !<
NUMFOR_KINDR, intent(OUT) :: result !<
real(dp), intent(OUT) :: abserr !<
real(dp), intent(OUT) :: resabs !<
real(dp), intent(OUT) :: resasc !<

real(dp) :: absc
real(dp) :: centr
real(dp) :: dhlgth
NUMFOR_KINDR :: fc
NUMFOR_KINDR :: fsum
NUMFOR_KINDR :: fval1
NUMFOR_KINDR :: fval2
NUMFOR_KINDR :: fv1(15)
NUMFOR_KINDR :: fv2(15)
real(dp) :: hlgth
integer j
integer jtw
integer jtwm1
NUMFOR_KINDR :: resg
NUMFOR_KINDR :: resk
NUMFOR_KINDR :: reskh
real(dp) :: wg(8)
real(dp) :: wgk(16)
real(dp) :: xgk(16)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 31-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 15-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 15-point Gauss rule
!
!           wgk    - weights of the 31-point Kronrod rule
!
!           wg     - weights of the 15-point Gauss rule
!
data xgk(1), xgk(2), xgk(3), xgk(4), xgk(5), xgk(6), xgk(7), xgk(8), &
  xgk(9), xgk(10), xgk(11), xgk(12), xgk(13), xgk(14), xgk(15), xgk(16)/ &
  9.980022986933971e-01, 9.879925180204854e-01, &
  9.677390756791391e-01, 9.372733924007059e-01, &
  8.972645323440819e-01, 8.482065834104272e-01, &
  7.904185014424659e-01, 7.244177313601700e-01, &
  6.509967412974170e-01, 5.709721726085388e-01, &
  4.850818636402397e-01, 3.941513470775634e-01, &
  2.991800071531688e-01, 2.011940939974345e-01, &
  1.011420669187175e-01, 0.0_8/
data wgk(1), wgk(2), wgk(3), wgk(4), wgk(5), wgk(6), wgk(7), wgk(8), &
  wgk(9), wgk(10), wgk(11), wgk(12), wgk(13), wgk(14), wgk(15), wgk(16)/ &
  5.377479872923349e-03, 1.500794732931612e-02, &
  2.546084732671532e-02, 3.534636079137585e-02, &
  4.458975132476488e-02, 5.348152469092809e-02, &
  6.200956780067064e-02, 6.985412131872826e-02, &
  7.684968075772038e-02, 8.308050282313302e-02, &
  8.856444305621177e-02, 9.312659817082532e-02, &
  9.664272698362368e-02, 9.917359872179196e-02, &
  1.007698455238756e-01, 1.013300070147915e-01/
data wg(1), wg(2), wg(3), wg(4), wg(5), wg(6), wg(7), wg(8)/ &
  3.075324199611727e-02, 7.036604748810812e-02, &
  1.071592204671719e-01, 1.395706779261543e-01, &
  1.662692058169939e-01, 1.861610000155622e-01, &
  1.984314853271116e-01, 2.025782419255613e-01/
#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif
!
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 15-point Gauss formula
!           resk   - result of the 31-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
centr = 0.5_8 * (a + b)
hlgth = 0.5_8 * (b - a)
dhlgth = abs(hlgth)
!
!  Compute the 31-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!
fc = NUMFOR_EVAL_F(f, centr)
resg = wg(8) * fc
resk = wgk(16) * fc
resabs = abs(resk)

do j = 1, 7
  jtw = j * 2
  absc = hlgth * xgk(jtw)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtw) = fval1
  fv2(jtw) = fval2
  fsum = fval1 + fval2
  resg = resg + wg(j) * fsum
  resk = resk + wgk(jtw) * fsum
  resabs = resabs + wgk(jtw) * (abs(fval1) + abs(fval2))
end do

do j = 1, 8
  jtwm1 = j * 2 - 1
  absc = hlgth * xgk(jtwm1)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtwm1) = fval1
  fv2(jtwm1) = fval2
  fsum = fval1 + fval2
  resk = resk + wgk(jtwm1) * fsum
  resabs = resabs + wgk(jtwm1) * (abs(fval1) + abs(fval2))
end do

reskh = resk * 0.5_8
resasc = wgk(16) * abs(fc - reskh)

do j = 1, 15
  resasc = resasc + wgk(j) * (abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
end do

result = resk * hlgth
resabs = resabs * dhlgth
resasc = resasc * dhlgth
abserr = abs((resk - resg) * hlgth)

IF (resasc /= 0.0_8 .and. abserr /= 0.0_8) &
  abserr = resasc * min(1.0_8, (2.0e+02 * abserr / resasc)**1.5_8)

IF (resabs > tiny(resabs) / (50._8 * epsilon(resabs)))&
  & abserr = max((epsilon(resabs) * 50._8) * resabs, abserr)

end PRNM(qk31)

#ifdef NF_FARG
PRNM(qk41) (ff, a, b, args, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qk41) (f, a, b, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
!*****************************************************************************80
!
!! QK41 carries out a 41 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
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
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!                       result is computed by applying the 41-point
!                       Gauss-Kronrod rule (resk) obtained by optimal
!                       addition of abscissae to the 20-point Gauss
!                       rule (resg).
!
!    Output, real(dp) :: ABSERR, an estimate of | I - RESULT |.
!
!    Output, real(dp) :: RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real(dp) :: RESASC, approximation to the integral | F-I/(B-A) |
!    over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 20-point Gauss formula
!           resk   - result of the 41-point Kronrod formula
!           reskh  - approximation to mean value of f over (a,b), i.e.
!                    to i/(b-a)
!
real(dp), intent(IN) :: a !<
real(dp), intent(IN) :: b !<
NUMFOR_KINDR, intent(OUT) :: result !<
real(dp), intent(OUT) :: abserr !<
real(dp), intent(OUT) :: resabs !<
real(dp), intent(OUT) :: resasc !<

real(dp) :: absc
real(dp) :: centr
real(dp) :: dhlgth
NUMFOR_KINDR :: fc
NUMFOR_KINDR :: fsum
NUMFOR_KINDR :: fval1
NUMFOR_KINDR :: fval2
NUMFOR_KINDR :: fv1(20)
NUMFOR_KINDR :: fv2(20)
real(dp) :: hlgth
integer j
integer jtw
integer jtwm1
NUMFOR_KINDR :: resg
NUMFOR_KINDR :: resk
NUMFOR_KINDR :: reskh

real(dp) :: wg(10)
real(dp) :: wgk(21)
real(dp) :: xgk(21)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 41-point Gauss-Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 20-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 20-point Gauss rule
!
!           wgk    - weights of the 41-point Gauss-Kronrod rule
!
!           wg     - weights of the 20-point Gauss rule
!
data xgk(1), xgk(2), xgk(3), xgk(4), xgk(5), xgk(6), xgk(7), xgk(8), &
  xgk(9), xgk(10), xgk(11), xgk(12), xgk(13), xgk(14), xgk(15), xgk(16), &
  xgk(17), xgk(18), xgk(19), xgk(20), xgk(21)/ &
  9.988590315882777e-01, 9.931285991850949e-01, &
  9.815078774502503e-01, 9.639719272779138e-01, &
  9.408226338317548e-01, 9.122344282513259e-01, &
  8.782768112522820e-01, 8.391169718222188e-01, &
  7.950414288375512e-01, 7.463319064601508e-01, &
  6.932376563347514e-01, 6.360536807265150e-01, &
  5.751404468197103e-01, 5.108670019508271e-01, &
  4.435931752387251e-01, 3.737060887154196e-01, &
  3.016278681149130e-01, 2.277858511416451e-01, &
  1.526054652409227e-01, 7.652652113349733e-02, &
  0.0_8/
data wgk(1), wgk(2), wgk(3), wgk(4), wgk(5), wgk(6), wgk(7), wgk(8), &
  wgk(9), wgk(10), wgk(11), wgk(12), wgk(13), wgk(14), wgk(15), wgk(16), &
  wgk(17), wgk(18), wgk(19), wgk(20), wgk(21)/ &
  3.073583718520532e-03, 8.600269855642942e-03, &
  1.462616925697125e-02, 2.038837346126652e-02, &
  2.588213360495116e-02, 3.128730677703280e-02, &
  3.660016975820080e-02, 4.166887332797369e-02, &
  4.643482186749767e-02, 5.094457392372869e-02, &
  5.519510534828599e-02, 5.911140088063957e-02, &
  6.265323755478117e-02, 6.583459713361842e-02, &
  6.864867292852162e-02, 7.105442355344407e-02, &
  7.303069033278667e-02, 7.458287540049919e-02, &
  7.570449768455667e-02, 7.637786767208074e-02, &
  7.660071191799966e-02/
data wg(1), wg(2), wg(3), wg(4), wg(5), wg(6), wg(7), wg(8), wg(9), wg(10)/ &
  1.761400713915212e-02, 4.060142980038694e-02, &
  6.267204833410906e-02, 8.327674157670475e-02, &
  1.019301198172404e-01, 1.181945319615184e-01, &
  1.316886384491766e-01, 1.420961093183821e-01, &
  1.491729864726037e-01, 1.527533871307259e-01/
!
#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif
!
centr = 0.5_8 * (a + b)
hlgth = 0.5_8 * (b - a)
dhlgth = abs(hlgth)
!
!  Compute 41-point Gauss-Kronrod approximation to the
!  the integral, and estimate the absolute error.
!
resg = 0.0_8
fc = NUMFOR_EVAL_F(f, centr)
resk = wgk(21) * fc
resabs = abs(resk)

do j = 1, 10
  jtw = j * 2
  absc = hlgth * xgk(jtw)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtw) = fval1
  fv2(jtw) = fval2
  fsum = fval1 + fval2
  resg = resg + wg(j) * fsum
  resk = resk + wgk(jtw) * fsum
  resabs = resabs + wgk(jtw) * (abs(fval1) + abs(fval2))
end do

do j = 1, 10
  jtwm1 = j * 2 - 1
  absc = hlgth * xgk(jtwm1)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtwm1) = fval1
  fv2(jtwm1) = fval2
  fsum = fval1 + fval2
  resk = resk + wgk(jtwm1) * fsum
  resabs = resabs + wgk(jtwm1) * (abs(fval1) + abs(fval2))
end do

reskh = resk * 0.5_8
resasc = wgk(21) * abs(fc - reskh)

do j = 1, 20
  resasc = resasc + wgk(j) * (abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
end do

result = resk * hlgth
resabs = resabs * dhlgth
resasc = resasc * dhlgth
abserr = abs((resk - resg) * hlgth)

if (resasc /= 0.0_8 .and. abserr /= 0.0_8) &
  abserr = resasc * min(1.0_8, (2.0e+02 * abserr / resasc)**1.5_8)

if (resabs > tiny(resabs) / (50._8 * epsilon(resabs))) then
  abserr = max((epsilon(resabs) * 50._8) * resabs, abserr)
end if

return
end PRNM(qk41)

!*****************************************************************************80
!
!! QK51 carries out a 51 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
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
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!                       result is computed by applying the 51-point
!                       Kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 25-point Gauss rule (resg).
!
!    Output, real(dp) :: ABSERR, an estimate of | I - RESULT |.
!
!    Output, real(dp) :: RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real(dp) :: RESASC, approximation to the integral | F-I/(B-A) |
!    over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 25-point Gauss formula
!           resk   - result of the 51-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
#ifdef NF_FARG
PRNM(qk51) (ff, a, b, args, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qk51) (f, a, b, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
real(dp), intent(IN) :: a !<
real(dp), intent(IN) :: b !<
NUMFOR_KINDR, intent(OUT) :: result !<
real(dp), intent(OUT) :: abserr !<
real(dp), intent(OUT) :: resabs !<
real(dp), intent(OUT) :: resasc !<

real(dp) :: absc
real(dp) :: centr
real(dp) :: dhlgth

NUMFOR_KINDR :: fc
NUMFOR_KINDR :: fsum
NUMFOR_KINDR :: fval1
NUMFOR_KINDR :: fval2
NUMFOR_KINDR :: fv1(25)
NUMFOR_KINDR :: fv2(25)
real(dp) :: hlgth
integer j
integer jtw
integer jtwm1
NUMFOR_KINDR :: resg
NUMFOR_KINDR :: resk
NUMFOR_KINDR :: reskh
real(dp) :: wg(13)
real(dp) :: wgk(26)
real(dp) :: xgk(26)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 51-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 25-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 25-point Gauss rule
!
!           wgk    - weights of the 51-point Kronrod rule
!
!           wg     - weights of the 25-point Gauss rule
!
data xgk(1), xgk(2), xgk(3), xgk(4), xgk(5), xgk(6), xgk(7), xgk(8), &
  xgk(9), xgk(10), xgk(11), xgk(12), xgk(13), xgk(14)/ &
  9.992621049926098e-01, 9.955569697904981e-01, &
  9.880357945340772e-01, 9.766639214595175e-01, &
  9.616149864258425e-01, 9.429745712289743e-01, &
  9.207471152817016e-01, 8.949919978782754e-01, &
  8.658470652932756e-01, 8.334426287608340e-01, &
  7.978737979985001e-01, 7.592592630373576e-01, &
  7.177664068130844e-01, 6.735663684734684e-01/
data xgk(15), xgk(16), xgk(17), xgk(18), xgk(19), xgk(20), xgk(21), &
  xgk(22), xgk(23), xgk(24), xgk(25), xgk(26)/ &
  6.268100990103174e-01, 5.776629302412230e-01, &
  5.263252843347192e-01, 4.730027314457150e-01, &
  4.178853821930377e-01, 3.611723058093878e-01, &
  3.030895389311078e-01, 2.438668837209884e-01, &
  1.837189394210489e-01, 1.228646926107104e-01, &
  6.154448300568508e-02, 0.0_8/
data wgk(1), wgk(2), wgk(3), wgk(4), wgk(5), wgk(6), wgk(7), wgk(8), &
  wgk(9), wgk(10), wgk(11), wgk(12), wgk(13), wgk(14)/ &
  1.987383892330316e-03, 5.561932135356714e-03, &
  9.473973386174152e-03, 1.323622919557167e-02, &
  1.684781770912830e-02, 2.043537114588284e-02, &
  2.400994560695322e-02, 2.747531758785174e-02, &
  3.079230016738749e-02, 3.400213027432934e-02, &
  3.711627148341554e-02, 4.008382550403238e-02, &
  4.287284502017005e-02, 4.550291304992179e-02/
data wgk(15), wgk(16), wgk(17), wgk(18), wgk(19), wgk(20), wgk(21), &
  wgk(22), wgk(23), wgk(24), wgk(25), wgk(26)/ &
  4.798253713883671e-02, 5.027767908071567e-02, &
  5.236288580640748e-02, 5.425112988854549e-02, &
  5.595081122041232e-02, 5.743711636156783e-02, &
  5.868968002239421e-02, 5.972034032417406e-02, &
  6.053945537604586e-02, 6.112850971705305e-02, &
  6.147118987142532e-02, 6.158081806783294e-02/
data wg(1), wg(2), wg(3), wg(4), wg(5), wg(6), wg(7), wg(8), wg(9), wg(10), &
  wg(11), wg(12), wg(13)/ &
  1.139379850102629e-02, 2.635498661503214e-02, &
  4.093915670130631e-02, 5.490469597583519e-02, &
  6.803833381235692e-02, 8.014070033500102e-02, &
  9.102826198296365e-02, 1.005359490670506e-01, &
  1.085196244742637e-01, 1.148582591457116e-01, &
  1.194557635357848e-01, 1.222424429903100e-01, &
  1.231760537267155e-01/

#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif
!
centr = 0.5_8 * (a + b)
hlgth = 0.5_8 * (b - a)
dhlgth = abs(hlgth)
!
!  Compute the 51-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!
fc = NUMFOR_EVAL_F(f, centr)
resg = wg(13) * fc
resk = wgk(26) * fc
resabs = abs(resk)

do j = 1, 12
  jtw = j * 2
  absc = hlgth * xgk(jtw)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtw) = fval1
  fv2(jtw) = fval2
  fsum = fval1 + fval2
  resg = resg + wg(j) * fsum
  resk = resk + wgk(jtw) * fsum
  resabs = resabs + wgk(jtw) * (abs(fval1) + abs(fval2))
end do

do j = 1, 13
  jtwm1 = j * 2 - 1
  absc = hlgth * xgk(jtwm1)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtwm1) = fval1
  fv2(jtwm1) = fval2
  fsum = fval1 + fval2
  resk = resk + wgk(jtwm1) * fsum
  resabs = resabs + wgk(jtwm1) * (abs(fval1) + abs(fval2))
end do

reskh = resk * 0.5_8
resasc = wgk(26) * abs(fc - reskh)

do j = 1, 25
  resasc = resasc + wgk(j) * (abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
end do

result = resk * hlgth
resabs = resabs * dhlgth
resasc = resasc * dhlgth
abserr = abs((resk - resg) * hlgth)

IF (resasc /= 0.0_8 .and. abserr /= 0.0_8)&
  & abserr = resasc * min(1.0_8, (2.0e+02 * abserr / resasc)**1.5_8)
IF (resabs > tiny(resabs) / (50._8 * epsilon(resabs)))&
  & abserr = max((epsilon(resabs) * 50._8) * resabs, abserr)

end PRNM(qk51)

#ifdef NF_FARG
PRNM(qk61) (ff, a, b, args, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qk61) (f, a, b, result, abserr, resabs, resasc)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
!*****************************************************************************80
!
!! QK61 carries out a 61 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
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
!    Output, real(dp) :: RESULT, the estimated value of the integral.
!                    result is computed by applying the 61-point
!                    Kronrod rule (resk) obtained by optimal addition of
!                    abscissae to the 30-point Gauss rule (resg).
!
!    Output, real(dp) :: ABSERR, an estimate of | I - RESULT |.
!
!    Output, real(dp) :: RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real(dp) :: RESASC, approximation to the integral | F-I/(B-A) |
!    over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 30-point Gauss rule
!           resk   - result of the 61-point Kronrod rule
!           reskh  - approximation to the mean value of f
!                    over (a,b), i.e. to i/(b-a)
!
real(dp), intent(IN) :: a !<
real(dp), intent(IN) :: b !<
NUMFOR_KINDR, intent(OUT) :: result !<
real(dp), intent(OUT) :: abserr !<
real(dp), intent(OUT) :: resabs !<
real(dp), intent(OUT) :: resasc !<

real(dp) :: absc
real(dp) :: centr
real(dp) :: dhlgth
NUMFOR_KINDR :: fc
NUMFOR_KINDR :: fsum
NUMFOR_KINDR :: fval1
NUMFOR_KINDR :: fval2
NUMFOR_KINDR :: fv1(30)
NUMFOR_KINDR :: fv2(30)
real(dp) :: hlgth
integer j
integer jtw
integer jtwm1
NUMFOR_KINDR :: resg
NUMFOR_KINDR :: resk
NUMFOR_KINDR :: reskh
real(dp) :: wg(15)
real(dp) :: wgk(31)
real(dp) :: xgk(31)
!
!           the abscissae and weights are given for the
!           interval (-1,1). because of symmetry only the positive
!           abscissae and their corresponding weights are given.
!
!           xgk   - abscissae of the 61-point Kronrod rule
!                   xgk(2), xgk(4)  ... abscissae of the 30-point
!                   Gauss rule
!                   xgk(1), xgk(3)  ... optimally added abscissae
!                   to the 30-point Gauss rule
!
!           wgk   - weights of the 61-point Kronrod rule
!
!           wg    - weigths of the 30-point Gauss rule
!
data xgk(1), xgk(2), xgk(3), xgk(4), xgk(5), xgk(6), xgk(7), xgk(8), &
  xgk(9), xgk(10)/ &
  9.994844100504906e-01, 9.968934840746495e-01, &
  9.916309968704046e-01, 9.836681232797472e-01, &
  9.731163225011263e-01, 9.600218649683075e-01, &
  9.443744447485600e-01, 9.262000474292743e-01, &
  9.055733076999078e-01, 8.825605357920527e-01/
data xgk(11), xgk(12), xgk(13), xgk(14), xgk(15), xgk(16), xgk(17), &
  xgk(18), xgk(19), xgk(20)/ &
  8.572052335460611e-01, 8.295657623827684e-01, &
  7.997278358218391e-01, 7.677774321048262e-01, &
  7.337900624532268e-01, 6.978504947933158e-01, &
  6.600610641266270e-01, 6.205261829892429e-01, &
  5.793452358263617e-01, 5.366241481420199e-01/
data xgk(21), xgk(22), xgk(23), xgk(24), xgk(25), xgk(26), xgk(27), &
  xgk(28), xgk(29), xgk(30), xgk(31)/ &
  4.924804678617786e-01, 4.470337695380892e-01, &
  4.004012548303944e-01, 3.527047255308781e-01, &
  3.040732022736251e-01, 2.546369261678898e-01, &
  2.045251166823099e-01, 1.538699136085835e-01, &
  1.028069379667370e-01, 5.147184255531770e-02, &
  0.0_8/
data wgk(1), wgk(2), wgk(3), wgk(4), wgk(5), wgk(6), wgk(7), wgk(8), &
  wgk(9), wgk(10)/ &
  1.389013698677008e-03, 3.890461127099884e-03, &
  6.630703915931292e-03, 9.273279659517763e-03, &
  1.182301525349634e-02, 1.436972950704580e-02, &
  1.692088918905327e-02, 1.941414119394238e-02, &
  2.182803582160919e-02, 2.419116207808060e-02/
data wgk(11), wgk(12), wgk(13), wgk(14), wgk(15), wgk(16), wgk(17), &
  wgk(18), wgk(19), wgk(20)/ &
  2.650995488233310e-02, 2.875404876504129e-02, &
  3.090725756238776e-02, 3.298144705748373e-02, &
  3.497933802806002e-02, 3.688236465182123e-02, &
  3.867894562472759e-02, 4.037453895153596e-02, &
  4.196981021516425e-02, 4.345253970135607e-02/
data wgk(21), wgk(22), wgk(23), wgk(24), wgk(25), wgk(26), wgk(27), &
  wgk(28), wgk(29), wgk(30), wgk(31)/ &
  4.481480013316266e-02, 4.605923827100699e-02, &
  4.718554656929915e-02, 4.818586175708713e-02, &
  4.905543455502978e-02, 4.979568342707421e-02, &
  5.040592140278235e-02, 5.088179589874961e-02, &
  5.122154784925877e-02, 5.142612853745903e-02, &
  5.149472942945157e-02/
data wg(1), wg(2), wg(3), wg(4), wg(5), wg(6), wg(7), wg(8)/ &
  7.968192496166606e-03, 1.846646831109096e-02, &
  2.878470788332337e-02, 3.879919256962705e-02, &
  4.840267283059405e-02, 5.749315621761907e-02, &
  6.597422988218050e-02, 7.375597473770521e-02/
data wg(9), wg(10), wg(11), wg(12), wg(13), wg(14), wg(15)/ &
  8.075589522942022e-02, 8.689978720108298e-02, &
  9.212252223778613e-02, 9.636873717464426e-02, &
  9.959342058679527e-02, 1.017623897484055e-01, &
  1.028526528935588e-01/
#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif
centr = 0.5_8 * (b + a)
hlgth = 0.5_8 * (b - a)
dhlgth = abs(hlgth)
!
!  Compute the 61-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!
resg = 0.0_8
fc = NUMFOR_EVAL_F(f, centr)
resk = wgk(31) * fc
resabs = abs(resk)

do j = 1, 15
  jtw = j * 2
  absc = hlgth * xgk(jtw)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtw) = fval1
  fv2(jtw) = fval2
  fsum = fval1 + fval2
  resg = resg + wg(j) * fsum
  resk = resk + wgk(jtw) * fsum
  resabs = resabs + wgk(jtw) * (abs(fval1) + abs(fval2))
end do

do j = 1, 15
  jtwm1 = j * 2 - 1
  absc = hlgth * xgk(jtwm1)
  fval1 = NUMFOR_EVAL_F(f, centr - absc)
  fval2 = NUMFOR_EVAL_F(f, centr + absc)
  fv1(jtwm1) = fval1
  fv2(jtwm1) = fval2
  fsum = fval1 + fval2
  resk = resk + wgk(jtwm1) * fsum
  resabs = resabs + wgk(jtwm1) * (abs(fval1) + abs(fval2))
end do

reskh = resk * 0.5_8
resasc = wgk(31) * abs(fc - reskh)

do j = 1, 30
  resasc = resasc + wgk(j) * (abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
end do

result = resk * hlgth
resabs = resabs * dhlgth
resasc = resasc * dhlgth
abserr = abs((resk - resg) * hlgth)

IF (resasc /= 0.0_8 .and. abserr /= 0.0_8) &
  &  abserr = resasc * min(1.0_8, (2.0e+02 * abserr / resasc)**1.5_8)
IF (resabs > tiny(resabs) / (50._8 * epsilon(resabs))) &
  & abserr = max((epsilon(resabs) * 50._8) * resabs, abserr)

end PRNM(qk61)

