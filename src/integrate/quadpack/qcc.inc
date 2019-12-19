!> qc25o returns integration rules for integrands with a COS or SIN factor.
!!
!! This routine estimates the integral
!!      I = integral of f(x) * w(x) over (a,b)
!! where
!!      w(x) = cos(omega*x)
!! or
!!      w(x) = sin(omega*x),
!! and estimates
!!      J = integral ( A <= X <= B ) |F(X)| dx.
!!
!! For small values of `omega` or small intervals `(a,b)` the 15-point Gauss-Kronrod rule is used.  In all
!! other cases a generalized Clenshaw-Curtis method is used, that is, a truncated Chebyshev expansion
!! of the function F is computed on (a,b), so that the integrand can be written as a sum of terms of
!! the form W(X)*T(K,X), where T(K,X) is the Chebyshev polynomial of degree K.  The Chebyshev moments
!! are computed with use of a linear recurrence relation.
!!
!!  Author:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner
!!
!!  Reference:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
!!    QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
!!    (Section: 2.2.3.4)
!!
#ifdef NF_FARG
PRNM(qc25o) (ff, a, b, omega, flgw, nrmom, maxp1, ksave, args, result, abserr, neval, resabs, resasc, momcom, chebmo)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qc25o) (f, a, b, omega, flgw, nrmom, maxp1, ksave, result, abserr, neval, resabs, resasc, momcom, chebmo)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
real(dp), intent(IN) :: a !< Lower limit of integration
real(dp), intent(IN) :: b !< Upper limit of integration
real(dp), intent(IN) :: omega !< Frequency of oscillation
integer, intent(IN) :: flgw !< Flag. If flgw=1 => w(x) = cos(omega*x). flgw=2 => sin(omega*x)
integer, intent(INOUT) :: nrmom !< Must be zero at the first call
integer, intent(IN) :: maxp1   !< Upper bound on the number of Chebyshev moments
integer, intent(IN) :: ksave !< flag indicating if the momenta for an interval have been computed
NUMFOR_KINDR, intent(OUT) :: result !<
real(dp), intent(OUT) :: abserr !<
integer, intent(OUT) :: neval !<
real(dp), intent(OUT) :: resabs !<
real(dp), intent(OUT) :: resasc !<
integer, intent(OUT) :: momcom !<
real(dp), dimension(maxp1, 25), intent(INOUT) :: chebmo !<
!! Examples:
!!

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
!    Input, integer FLGW, indicates which weight function is to be used
!    = 1, w(x) = cos(omega*x)
!    = 2, w(x) = sin(omega*x)
!
!    ?, integer NRMOM, the length of interval (a,b) is equal to the length of
!    the original integration interval divided by 2**nrmom (we suppose that the
!    routine is used in an adaptive integration process, otherwise set nrmom = 0).
!    nrmom must be zero at the first call.
!
!           maxp1  - integer
!                    gives an upper bound on the number of Chebyshev
!                    moments which can be stored, i.e. for the intervals
!                    of lengths abs(bb-aa)*2**(-l), l = 0,1,2, ...,
!                    maxp1-2.
!
!           ksave  - integer
!                    key which is one when the moments for the
!                    current interval have been computed
!
!    Output, real(8) :: RESULT, the estimated value of the integral.
!
!           abserr - real
!                    estimate of the modulus of the absolute
!                    error, which should equal or exceed abs(i-result)
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!    Output, real(8) :: RESABS, approximation to the integral J.
!
!    Output, real(8) :: RESASC, approximation to the integral of abs(F-I/(B-A)).
!
!         on entry and return
!           momcom - integer
!                    for each interval length we need to compute
!                    the Chebyshev moments. momcom counts the number
!                    of intervals for which these moments have already
!                    been computed. if nrmom < momcom or ksave = 1,
!                    the Chebyshev moments for the interval (a,b)
!                    have already been computed and stored, otherwise
!                    we compute them and we increase momcom.
!
!           chebmo - real
!                    array of dimension at least (maxp1,25) containing
!                    the modified Chebyshev moments for the first momcom
!                    interval lengths
!
!  Local parameters:
!
!    maxp1 gives an upper bound
!           on the number of Chebyshev moments which can be
!           computed, i.e. for the interval (bb-aa), ...,
!           (bb-aa)/2**(maxp1-2).
!           should this number be altered, the first dimension of
!           chebmo needs to be adapted.
!
!    x contains the values cos(k*pi/24)
!           k = 1, ...,11, to be used for the Chebyshev expansion of f
!
!           centr  - mid point of the integration interval
!           hlgth  - half length of the integration interval
!           fval   - value of the function f at the points
!                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5
!                    k = 0, ...,24
!           cheb12 - coefficients of the Chebyshev series expansion
!                    of degree 12, for the function f, in the
!                    interval (a,b)
!           cheb24 - coefficients of the Chebyshev series expansion
!                    of degree 24, for the function f, in the
!                    interval (a,b)
!           resc12 - approximation to the integral of
!                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
!                    over (-1,+1), using the Chebyshev series
!                    expansion of degree 12
!           resc24 - approximation to the same integral, using the
!                    Chebyshev series expansion of degree 24
!           ress12 - the analogue of resc12 for the sine
!           ress24 - the analogue of resc24 for the sine
!
real(8) :: ac
real(8) :: an
real(8) :: an2
real(8) :: as
real(8) :: asap
real(8) :: ass
real(8) :: centr

NUMFOR_KINDR :: cheb12(13)
NUMFOR_KINDR :: cheb24(25)
real(8) :: conc
real(8) :: cons
real(8) :: cospar
real(8) :: d(28)
real(8) :: d1(28)
real(8) :: d2(28)
real(8) :: d3(28)
real(8) :: estc
real(8) :: ests

NUMFOR_KINDR :: fval(25)
real(8) :: hlgth
integer i
integer isym
integer j
integer k
integer m
integer, parameter :: nmac = 28
integer noeq1
integer noequ
real(8) :: parint
real(8) :: par2
real(8) :: par22
! real(8) :: omega_v(1)

NUMFOR_KINDR :: resc12
NUMFOR_KINDR :: resc24
NUMFOR_KINDR :: ress12
NUMFOR_KINDR :: ress24
real(8) :: sinpar
real(8) :: v(28)
real(8), parameter, dimension(11) :: x = [.9914448613738104_8, .9659258262890683_8, .9238795325112868_8, &
  & 0.8660254037844386_8, .7933533402912352_8, .7071067811865475_8, .6087614290087206_8, &
  & 0.5000000000000000_8, .3826834323650898_8, .2588190451025208_8, .1305261922200516_8]

#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif

centr = 0.5_8 * (b + a)
hlgth = 0.5_8 * (b - a)
parint = omega * hlgth
! omega_v(1) = omega
!
!  Compute the integral using the 15-point Gauss-Kronrod
!  formula if the value of the parameter in the integrand
!  is small or if the length of the integration interval
!  is less than (bb-aa)/2**(maxp1-2), where (aa,bb) is the
!  original integration interval.
!
if (abs(parint) <= 2.0_8) then
  call qk15w(f, qwgto, [omega], flgw, a, b, result, abserr, resabs, resasc)
  neval = 15
  return
end if
!
!  Compute the integral using the generalized clenshaw-curtis method.
!
conc = hlgth * cos(centr * omega)
cons = hlgth * sin(centr * omega)
resasc = huge(resasc)
neval = 25
!
!  Check whether the Chebyshev moments for this interval
!  have already been computed.
!
IF (nrmom < momcom .or. ksave == 1) go to 140
!
!  Compute a new set of Chebyshev moments.
!
m = momcom + 1
par2 = parint * parint
par22 = par2 + 2.0_8
sinpar = sin(parint)
cospar = cos(parint)
!
!  Compute the Chebyshev moments with respect to cosine.
!
v(1) = 2.0_8 * sinpar / parint
v(2) = (8.0_8 * cospar + (par2 + par2 - 8.0_8) * sinpar / parint) / par2
v(3) = (3.2e+01 * (par2 - 1.2e+01) * cospar + (2.0_8 * &
                                               ((par2 - 8.0e+01) * par2 + 1.92e+02) * sinpar) / &
        parint) / (par2 * par2)
ac = 8.0_8 * cospar
as = 24.0_8 * parint * sinpar

IF (abs(parint) > 24.0_8) go to 70
!
!  Compute the Chebyshev moments as the solutions of a boundary value
!  problem with one initial value (v(3)) and one end value computed
!  using an asymptotic formula.
!
noequ = nmac - 3
noeq1 = noequ - 1
an = 6.0_8

do k = 1, noeq1
  an2 = an * an
  d(k) = -2.0_8 * (an2 - 4.0_8) * (par22 - an2 - an2)
  d2(k) = (an - 1.0_8) * (an - 2.0_8) * par2
  d1(k) = (an + 3.0_8) * (an + 4.0_8) * par2
  v(k + 3) = as - (an2 - 4.0_8) * ac
  an = an + 2.0_8
end do

an2 = an * an
d(noequ) = -2.0_8 * (an2 - 4.0_8) * (par22 - an2 - an2)
v(noequ + 3) = as - (an2 - 4.0_8) * ac
v(4) = v(4) - 5.6e+01 * par2 * v(3)
ass = parint * sinpar
asap = (((((2.10e+02 * par2 - 1.0_8) * cospar - (1.05e+02 * par2 &
                                                 - 6.3e+01) * ass) / an2 - (1.0_8 - 1.5e+01 * par2) * cospar &
          + 1.5e+01 * ass) / an2 - cospar + 3.0_8 * ass) / an2 - cospar) / an2
v(noequ + 3) = v(noequ + 3) - 2.0_8 * asap * par2 * (an - 1.0_8) * &
               (an - 2.0_8)
!
!  Solve the tridiagonal system by means of Gaussian
!  elimination with partial pivoting.
!
d3(1:noequ) = 0.0_8

d2(noequ) = 0.0_8

do i = 1, noeq1

  if (abs(d1(i)) > abs(d(i))) then
    an = d1(i)
    d1(i) = d(i)
    d(i) = an
    an = d2(i)
    d2(i) = d(i + 1)
    d(i + 1) = an
    d3(i) = d2(i + 1)
    d2(i + 1) = 0.0_8
    an = v(i + 4)
    v(i + 4) = v(i + 3)
    v(i + 3) = an
  end if

  d(i + 1) = d(i + 1) - d2(i) * d1(i) / d(i)
  d2(i + 1) = d2(i + 1) - d3(i) * d1(i) / d(i)
  v(i + 4) = v(i + 4) - v(i + 3) * d1(i) / d(i)

end do

v(noequ + 3) = v(noequ + 3) / d(noequ)
v(noequ + 2) = (v(noequ + 2) - d2(noeq1) * v(noequ + 3)) / d(noeq1)

do i = 2, noeq1
  k = noequ - i
  v(k + 3) = (v(k + 3) - d3(k) * v(k + 5) - d2(k) * v(k + 4)) / d(k)
end do

go to 90
!
!  Compute the Chebyshev moments by means of forward recursion
!
70 continue

an = 4.0_8

do i = 4, 13
  an2 = an * an
  v(i) = ((an2 - 4.0_8) * (2.0_8 * (par22 - an2 - an2) * v(i - 1) - ac) &
          + as - par2 * (an + 1.0_8) * (an + 2.0_8) * v(i - 2)) / &
         (par2 * (an - 1.0_8) * (an - 2.0_8))
  an = an + 2.0_8
end do

90 continue

do j = 1, 13
  chebmo(m, 2 * j - 1) = v(j)
end do
!
!  Compute the Chebyshev moments with respect to sine.
!
v(1) = 2.0_8 * (sinpar - parint * cospar) / par2
v(2) = (1.8e+01 - 4.8e+01 / par2) * sinpar / par2 &
       + (-2.0_8 + 4.8e+01 / par2) * cospar / parint
ac = -2.4e+01 * parint * cospar
as = -8.0_8 * sinpar
chebmo(m, 2) = v(1)
chebmo(m, 4) = v(2)
if (abs(parint) <= 2.4e+01) then

  do k = 3, 12
    an = k
    chebmo(m, 2 * k) = -sinpar / (an * (2.0_8 * an - 2.0_8)) &
                       - 2.5e-01 * parint * (v(k + 1) / an - v(k) / (an - 1.0_8))
  end do
  !
  !  Compute the Chebyshev moments by means of forward recursion.
  !
else

  an = 3.0_8

  do i = 3, 12
    an2 = an * an
    v(i) = ((an2 - 4.0_8) * (2.0_8 * (par22 - an2 - an2) * v(i - 1) + as) &
            + ac - par2 * (an + 1.0_8) * (an + 2.0_8) * v(i - 2)) &
           / (par2 * (an - 1.0_8) * (an - 2.0_8))
    an = an + 2.0_8
    chebmo(m, 2 * i) = v(i)
  end do

end if

140 continue

IF (nrmom < momcom) m = nrmom + 1

IF (momcom < maxp1 - 1 .and. nrmom >= momcom) momcom = momcom + 1

!
!  Compute the coefficients of the Chebyshev expansions
!  of degrees 12 and 24 of the function F.
!
fval(1) = 0.5_8 * NUMFOR_EVAL_F(f, centr + hlgth)
fval(13) = NUMFOR_EVAL_F(f, centr)
fval(25) = 0.5_8 * NUMFOR_EVAL_F(f, centr - hlgth)

do i = 2, 12
  isym = 26 - i
  fval(i) = NUMFOR_EVAL_F(f, hlgth * x(i - 1) + centr)
  fval(isym) = NUMFOR_EVAL_F(f, centr - hlgth * x(i - 1))
end do

call qcheb(x, fval, cheb12, cheb24)
!
!  Compute the integral and error estimates.
!
resc12 = cheb12(13) * chebmo(m, 13)
ress12 = 0.0_8
estc = abs(cheb24(25) * chebmo(m, 25)) + &
  &    abs((cheb12(13) - cheb24(13)) * chebmo(m, 13))
ests = 0.0_8
k = 11

do j = 1, 6
  resc12 = resc12 + cheb12(k) * chebmo(m, k)
  ress12 = ress12 + cheb12(k + 1) * chebmo(m, k + 1)
  estc = estc + abs((cheb12(k) - cheb24(k)) * chebmo(m, k))
  ests = ests + abs((cheb12(k + 1) - cheb24(k + 1)) * chebmo(m, k + 1))
  k = k - 2
end do

resc24 = cheb24(25) * chebmo(m, 25)
ress24 = 0.0_8
resabs = abs(cheb24(25))
k = 23

do j = 1, 12

  resc24 = resc24 + cheb24(k) * chebmo(m, k)
  ress24 = ress24 + cheb24(k + 1) * chebmo(m, k + 1)
  resabs = resabs + abs(cheb24(k)) + abs(cheb24(k + 1))

  if (j <= 5) then
    estc = estc + abs(cheb24(k) * chebmo(m, k))
    ests = ests + abs(cheb24(k + 1) * chebmo(m, k + 1))
  end if

  k = k - 2

end do

resabs = resabs * abs(hlgth)
if (flgw == 1) then
  result = conc * resc24 - cons * ress24
  abserr = abs(conc * estc) + abs(cons * ests)
else
  result = conc * ress24 + cons * resc24
  abserr = abs(conc * ests) + abs(cons * estc)
end if
end PRNM(qc25o)

!> QC25S returns rules for algebraico-logarithmic end point singularities.
!!
!!  This routine computes the integral \f$ I = \int_{bl}^{br} f(X) w(X) dx \f$ with error estimate,
!!  where the weight function w(x) has a singular behavior of algebraico-logarithmic type at the points
!!  a and/or b.  The interval (bl,br) is a subinterval of (a,b).
!!
!! The estimated value of the integral is computed by using a generalized clenshaw-curtis method
!! if b1 = a or br = b. In all other cases the 15-point Kronrod rule is applied, obtained by optimal
!! addition of abscissae to the 7-point Gauss rule.
!!
!!  Author:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner
!!
!!  Reference:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
!!    QUADPACK, a Subroutine Package for Automatic Flgwation, Springer Verlag, 1983
!!    (Section 2.2.3.5)
!!
#ifdef NF_FARG
PRNM(qc25s) (ff, bl, br, params, ri, rj, rg, rh, args, result, abserr, resasc, flagw, neval)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qc25s) (f, bl, br, params, ri, rj, rg, rh, result, abserr, resasc, flagw, neval)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
real(dp), intent(IN) :: bl !<
real(dp), intent(IN) :: br !<
real(dp), dimension(4), intent(IN) :: params !< Vector with parameters for the weigth functions:
!! - params(1) -> a (bl <= a < b)
!! - params(2) -> b (b <= br)
!! - params(3) -> alfa
!! - params(4) -> beta
real(dp), dimension(25), intent(IN) :: ri !< Modified Chebyshev moments for generalized Clenshaw-Curtis method
real(dp), dimension(25), intent(IN) :: rj !<
real(dp), dimension(25), intent(IN) :: rg !<
real(dp), dimension(25), intent(IN) :: rh !<
NUMFOR_KINDR, intent(OUT) :: result !< I = integrate(f(x)*w(x), x, bl,br)
real(dp), intent(OUT) :: abserr !< Estimation of absolute error
real(dp), intent(OUT) :: resasc !< integrate(abs(F(x)*w(x)-I/(b-a)),x,bl,br)
integer, intent(IN) :: flagw !< Determines the weight function
!!     - 1: -> w(x) = (x-a)**alfa*(b-x)**beta
!!     - 2: -> w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
!!     - 3: -> w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
!!     - 4: -> w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
integer, intent(OUT) :: neval !< Number of function evaluations performed
!! Examples:
!!
!*****************************************************************************80
!  Local Parameters:
!
!           fval   - value of the function f at the points
!                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
!                    k = 0, ..., 24
!           cheb12 - coefficients of the Chebyshev series expansion
!                    of degree 12, for the function f, in the interval
!                    (bl,br)
!           cheb24 - coefficients of the Chebyshev series expansion
!                    of degree 24, for the function f, in the interval
!                    (bl,br)
!           res12  - approximation to the integral obtained from cheb12
!           res24  - approximation to the integral obtained from cheb24
!           qwgts  - external function subprogram defining the four
!                    possible weight functions
!           hlgth  - half-length of the interval (bl,br)
!           centr  - mid point of the interval (bl,br)
!
!           the vector x contains the values cos(k*pi/24)
!           k = 1, ..., 11, to be used for the computation of the
!           Chebyshev series expansion of f.
!
! real(8), external :: qwgts
real(8) :: centr
NUMFOR_KINDR :: cheb12(13)
NUMFOR_KINDR :: cheb24(25)
real(8) :: dc
real(8) :: factor
real(8) :: fix
NUMFOR_KINDR :: fval(25)
real(8) :: hlgth
integer i
integer isym
real(8) :: resabs             ! Not used
NUMFOR_KINDR :: res12
NUMFOR_KINDR :: res24
real(8) :: u
real(dp) :: a, b, alfa, beta
real(8), dimension(11), parameter :: x = [ &
  & 9.914448613738104e-01, 9.659258262890683e-01, 9.238795325112868e-01, &
  & 8.660254037844386e-01, 7.933533402912352e-01, 7.071067811865475e-01, &
  & 6.087614290087206e-01, 5.000000000000000e-01, 3.826834323650898e-01, &
  & 2.588190451025208e-01, 1.305261922200516e-01]
#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif
neval = 25
! associate (a=>params(1), b=>params(2), alfa=>params(3), beta=>params(4))
a = params(1); b = params(2)
alfa = params(3); beta = params(4)
if (bl == a .and. (alfa /= 0.0_8 .or. flagw == 2 .or. flagw == 4)) then
  go to 10
end if

if (br == b .and. (beta /= 0.0_8 .or. flagw == 3 .or. flagw == 4)) &
  go to 140
!
!  If a > bl and b < br, apply the 15-point Gauss-Kronrod scheme.
!
call qk15w(f, qwgts, params, flagw, bl, br, result, abserr, resabs, resasc)
neval = 15
return
!
!  This part of the program is executed only if a = bl.
!
!  Compute the Chebyshev series expansion of the function
!  f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta*f(0.5*(br-a)*x+0.5*(br+a))
!
10 continue

hlgth = 0.5_8 * (br - bl)
centr = 0.5_8 * (br + bl)
fix = b - centr
fval(1) = 0.5_8 * NUMFOR_EVAL_F(f, hlgth + centr) * (fix - hlgth)**beta
fval(13) = NUMFOR_EVAL_F(f, centr) * (fix**beta)
fval(25) = 0.5_8 * NUMFOR_EVAL_F(f, centr - hlgth) * (fix + hlgth)**beta

do i = 2, 12
  u = hlgth * x(i - 1)
  isym = 26 - i
  fval(i) = NUMFOR_EVAL_F(f, u + centr) * (fix - u)**beta
  fval(isym) = NUMFOR_EVAL_F(f, centr - u) * (fix + u)**beta
end do

factor = hlgth**(alfa + 1.0_8)
result = 0.0_8
abserr = 0.0_8
res12 = 0.0_8
res24 = 0.0_8

if (flagw > 2) go to 70

call qcheb(x, fval, cheb12, cheb24)
!
!  flagw = 1  (or 2)
!
do i = 1, 13
  res12 = res12 + cheb12(i) * ri(i)
  res24 = res24 + cheb24(i) * ri(i)
end do

do i = 14, 25
  res24 = res24 + cheb24(i) * ri(i)
end do

if (flagw == 1) go to 130
!
!  flagw = 2
!
dc = log(br - bl)
result = res24 * dc
abserr = abs((res24 - res12) * dc)
res12 = 0.0_8
res24 = 0.0_8

do i = 1, 13
  res12 = res12 + cheb12(i) * rg(i)
  res24 = res24 + cheb24(i) * rg(i)
end do

do i = 14, 25
  res24 = res24 + cheb24(i) * rg(i)
end do

go to 130
!
!  Compute the Chebyshev series expansion of the function
!  F4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
!
70 continue

fval(1) = fval(1) * log(fix - hlgth)
fval(13) = fval(13) * log(fix)
fval(25) = fval(25) * log(fix + hlgth)

do i = 2, 12
  u = hlgth * x(i - 1)
  isym = 26 - i
  fval(i) = fval(i) * log(fix - u)
  fval(isym) = fval(isym) * log(fix + u)
end do

call qcheb(x, fval, cheb12, cheb24)
!
!  flagw = 3  (or 4)
!
do i = 1, 13
  res12 = res12 + cheb12(i) * ri(i)
  res24 = res24 + cheb24(i) * ri(i)
end do

do i = 14, 25
  res24 = res24 + cheb24(i) * ri(i)
end do

if (flagw == 3) then
  go to 130
end if
!
!  flagw = 4
!
dc = log(br - bl)
result = res24 * dc
abserr = abs((res24 - res12) * dc)
res12 = 0.0_8
res24 = 0.0_8

do i = 1, 13
  res12 = res12 + cheb12(i) * rg(i)
  res24 = res24 + cheb24(i) * rg(i)
end do

do i = 14, 25
  res24 = res24 + cheb24(i) * rg(i)
end do

130 continue

result = (result + res24) * factor
abserr = (abserr + abs(res24 - res12)) * factor
return
!
!  This part of the program is executed only if b = br.
!
!  Compute the Chebyshev series expansion of the function
!  f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa*f(0.5*(b-bl)*x+0.5*(b+bl))
!
140 continue

hlgth = 0.5_8 * (br - bl)
centr = 0.5_8 * (br + bl)
fix = centr - a
fval(1) = 0.5_8 * NUMFOR_EVAL_F(f, hlgth + centr) * (fix + hlgth)**alfa
fval(13) = NUMFOR_EVAL_F(f, centr) * (fix**alfa)
fval(25) = 0.5_8 * NUMFOR_EVAL_F(f, centr - hlgth) * (fix - hlgth)**alfa

do i = 2, 12
  u = hlgth * x(i - 1)
  isym = 26 - i
  fval(i) = NUMFOR_EVAL_F(f, u + centr) * (fix + u)**alfa
  fval(isym) = NUMFOR_EVAL_F(f, centr - u) * (fix - u)**alfa
end do

factor = hlgth**(beta + 1.0_8)
result = 0.0_8
abserr = 0.0_8
res12 = 0.0_8
res24 = 0.0_8

if (flagw == 2 .or. flagw == 4) then
  go to 200
end if
!
!  flagw = 1  (or 3)
!
call qcheb(x, fval, cheb12, cheb24)

do i = 1, 13
  res12 = res12 + cheb12(i) * rj(i)
  res24 = res24 + cheb24(i) * rj(i)
end do

do i = 14, 25
  res24 = res24 + cheb24(i) * rj(i)
end do

if (flagw == 1) go to 260
!
!  flagw = 3
!
dc = log(br - bl)
result = res24 * dc
abserr = abs((res24 - res12) * dc)
res12 = 0.0_8
res24 = 0.0_8

do i = 1, 13
  res12 = res12 + cheb12(i) * rh(i)
  res24 = res24 + cheb24(i) * rh(i)
end do

do i = 14, 25
  res24 = res24 + cheb24(i) * rh(i)
end do

go to 260
!
!  Compute the Chebyshev series expansion of the function
!  f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
!
200 continue

fval(1) = fval(1) * log(hlgth + fix)
fval(13) = fval(13) * log(fix)
fval(25) = fval(25) * log(fix - hlgth)

do i = 2, 12
  u = hlgth * x(i - 1)
  isym = 26 - i
  fval(i) = fval(i) * log(u + fix)
  fval(isym) = fval(isym) * log(fix - u)
end do

call qcheb(x, fval, cheb12, cheb24)
!
!  flagw = 2  (or 4)
!
do i = 1, 13
  res12 = res12 + cheb12(i) * rj(i)
  res24 = res24 + cheb24(i) * rj(i)
end do

do i = 14, 25
  res24 = res24 + cheb24(i) * rj(i)
end do

if (flagw == 2) go to 260

dc = log(br - bl)
result = res24 * dc
abserr = abs((res24 - res12) * dc)
res12 = 0.0_8
res24 = 0.0_8
!
!  flagw = 4
!
do i = 1, 13
  res12 = res12 + cheb12(i) * rh(i)
  res24 = res24 + cheb24(i) * rh(i)
end do

do i = 14, 25
  res24 = res24 + cheb24(i) * rh(i)
end do

260 continue

result = (result + res24) * factor
abserr = (abserr + abs(res24 - res12)) * factor
! end associate

end PRNM(qc25s)

!> qc25c returns integration rules for Cauchy Principal Value integrals.
!!
!!    This routine estimates
!!      I = integral of F(X) * W(X) over (a,b)
!!    with error estimate, where
!!      w(x) = 1/(x-c)
!!
!! The result is computed by using a generalized Clenshaw-Curtis method if C lies within ten percent
!! of the integration interval.  In the other case the 15-point Kronrod rule obtained by optimal
!! addition of abscissae to the 7-point Gauss rule, is applied.
!!
!!  Author:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner
!!
!!  Reference:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
!!    QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983.
!!    (Section 2.2.3.6)
!!
#ifdef NF_FARG
PRNM(qc25c) (ff, a, b, c, args, result, abserr, krul, neval)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qc25c) (f, a, b, c, result, abserr, krul, neval)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
real(dp), intent(IN) :: a !<
real(dp), intent(IN) :: b !<
real(dp), intent(IN) :: c !< Position of the pole in the weight function
NUMFOR_KINDR, intent(OUT) :: result !<
real(dp), intent(OUT) :: abserr !<
integer, intent(INOUT) :: krul !<
integer, intent(OUT) :: neval !<
!! Examples:
!!

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
!    Input, real(8) :: C, the parameter in the weight function.
!
!    Output, real(8) :: RESULT, the estimated value of the integral.

!    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
!
!           krul   - integer
!                    key which is decreased by 1 if the 15-point
!                    Gauss-Kronrod scheme has been used
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!  Local parameters:
!
!           fval   - value of the function f at the points
!                    cos(k*pi/24),  k = 0, ..., 24
!           cheb12 - Chebyshev series expansion coefficients, for the
!                    function f, of degree 12
!           cheb24 - Chebyshev series expansion coefficients, for the
!                    function f, of degree 24
!           res12  - approximation to the integral corresponding to the
!                    use of cheb12
!           res24  - approximation to the integral corresponding to the
!                    use of cheb24
!           qwgtc  - external function subprogram defining the weight
!                    function
!           hlgth  - half-length of the interval
!           centr  - mid point of the interval
!
real(8) :: ak22
real(8) :: amom0, amom1, amom2
real(8) :: cc, centr
NUMFOR_KINDR :: cheb12(13)
NUMFOR_KINDR :: cheb24(25)
NUMFOR_KINDR :: fval(25)
real(8) :: hlgth
integer i
integer isym
integer k
integer :: kp
real(8) :: resabs, resasc
NUMFOR_KINDR :: res12, res24
real(8) :: u
real(8), parameter, dimension(11) :: x = [9.914448613738104e-01_8, 9.659258262890683e-01_8, &
  & 9.238795325112868e-01_8, 8.660254037844386e-01_8, 7.933533402912352e-01_8, 7.071067811865475e-01_8,&
  & 6.087614290087206e-01_8, 5.000000000000000e-01_8, 3.826834323650898e-01_8, 2.588190451025208e-01_8,&
  & 1.305261922200516e-01_8]

#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif
!
!  Check the position of C.
!
cc = (2.0_8 * c - b - a) / (b - a)
!
!  Apply the 15-point Gauss-Kronrod scheme.
!
if (abs(cc) >= 1.1_8) then
  krul = krul - 1
  call qk15w(f, qwgtc, [c], kp, a, b, result, abserr, resabs, resasc)
  neval = 15
  if (resasc == abserr) then
    krul = krul + 1
  end if
  return
end if
!
!  Use the generalized Clenshaw-Curtis method.
!
hlgth = 0.5_8 * (b - a)
centr = 0.5_8 * (b + a)
neval = 25
fval(1) = 0.5_8 * NUMFOR_EVAL_F(f, hlgth + centr)
fval(13) = NUMFOR_EVAL_F(f, centr)
fval(25) = 0.5_8 * NUMFOR_EVAL_F(f, centr - hlgth)

do i = 2, 12
  u = hlgth * x(i - 1)
  isym = 26 - i
  fval(i) = NUMFOR_EVAL_F(f, u + centr)
  fval(isym) = NUMFOR_EVAL_F(f, centr - u)
end do
!
!  Compute the Chebyshev series expansion.
!
call qcheb(x, fval, cheb12, cheb24)
!
!  The modified Chebyshev moments are computed by forward
!  recursion, using AMOM0 and AMOM1 as starting values.
!
amom0 = log(abs((1.0_8 - cc) / (1.0_8 + cc)))
amom1 = 2.0_8 + cc * amom0
res12 = cheb12(1) * amom0 + cheb12(2) * amom1
res24 = cheb24(1) * amom0 + cheb24(2) * amom1

do k = 3, 13
  amom2 = 2.0_8 * cc * amom1 - amom0
  ak22 = (k - 2) * (k - 2)
  if ((k / 2) * 2 == k) then
    amom2 = amom2 - 4.0_8 / (ak22 - 1.0_8)
  end if
  res12 = res12 + cheb12(k) * amom2
  res24 = res24 + cheb24(k) * amom2
  amom0 = amom1
  amom1 = amom2
end do

do k = 14, 25
  amom2 = 2.0_8 * cc * amom1 - amom0
  ak22 = (k - 2) * (k - 2)
  if ((k / 2) * 2 == k) then
    amom2 = amom2 - 4.0_8 / (ak22 - 1.0_8)
  end if
  res24 = res24 + cheb24(k) * amom2
  amom0 = amom1
  amom1 = amom2
end do

result = res24
abserr = abs(res24 - res12)

end PRNM(qc25c)

!> qfour estimates the integrals of oscillatory functions.
!!
!! This routine calculates an approximation RESULT to a definite integral
!!   I = integral of F(X) * COS(OMEGA*X)
!! or
!!   I = integral of F(X) * SIN(OMEGA*X)
!! over (A,B), hopefully satisfying:
!!   | I - RESULT | <= max ( epsabs, epsrel * |I| ) ).
!!
!! QFOUR is called by QAWO and QAWF.  It can also be called directly in a user-written program.  In the
!! latter case it is possible for the user to determine the first dimension of array CHEBMO(MAXP1,25).
!! See also parameter description of `maxp1`.  Additionally see parameter description of `icall` for
!! eventually re-using Chebyshev moments computed during former call on subinterval of equal length
!! abs(B-A).
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
PRNM(qfour) (ff, a, b, omega, flgw, epsabs, epsrel, limit, icall, maxp1, &
  & args, result, abserr, neval, ier, alist, blist, rlist, elist, iord, nnlog, momcom, chebmo)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(qfour) (f, a, b, omega, flgw, epsabs, epsrel, limit, icall, maxp1, &
  & result, abserr, neval, ier, alist, blist, rlist, elist, iord, nnlog, momcom, chebmo)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
real(dp), intent(IN) :: a !<
real(dp), intent(IN) :: b !<
real(dp), intent(IN) :: omega !<
integer, intent(IN) :: flgw !<
real(dp), intent(IN) :: epsabs !<
real(dp), intent(IN) :: epsrel !<
integer, intent(IN) :: limit !<
integer, intent(IN) :: icall !<
integer, intent(IN) :: maxp1 !<
NUMFOR_KINDR, intent(OUT) :: result !<
real(dp), intent(OUT) :: abserr !<
integer, intent(OUT) :: neval !<
integer, intent(OUT) :: ier !<
real(dp), dimension(limit), intent(OUT) :: alist !<
real(dp), dimension(limit), intent(OUT) :: blist !<
NUMFOR_KINDR, dimension(limit), intent(OUT) :: rlist !<
real(dp), dimension(limit), intent(OUT) :: elist !<
integer, dimension(limit), intent(OUT) :: iord !<
integer, dimension(limit), intent(OUT) :: nnlog !<
integer, intent(INOUT) :: momcom !<
real(dp), dimension(maxp1, 25), intent(INOUT) :: chebmo !<
!! Examples:
!!

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
!    Input, real(8) :: OMEGA, the multiplier of X in the weight function.
!
!    Input, integer FLGW, indicates the weight functions to be used.
!    = 1, w(x) = cos(omega*x)
!    = 2, w(x) = sin(omega*x)
!
!    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer LIMIT, the maximum number of subintervals of [A,B]
!    that can be generated.
!
!    icall  - integer
!                     if qfour is to be used only once, ICALL must
!                     be set to 1.  assume that during this call, the
!                     Chebyshev moments (for clenshaw-curtis integration
!                     of degree 24) have been computed for intervals of
!                     lenghts (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
!                     the Chebyshev moments already computed can be
!                     re-used in subsequent calls, if qfour must be
!                     called twice or more times on intervals of the
!                     same length abs(b-a). from the second call on, one
!                     has to put then ICALL > 1.
!                     if ICALL < 1, the routine will end with ier = 6.
!
!            maxp1  - integer
!                     gives an upper bound on the number of
!                     Chebyshev moments which can be stored, i.e.
!                     for the intervals of lenghts abs(b-a)*2**(-l),
!                     l=0,1, ..., maxp1-2, maxp1 >= 1.
!                     if maxp1 < 1, the routine will end with ier = 6.
!                     increasing (decreasing) the value of maxp1
!                     decreases (increases) the computational time but
!                     increases (decreases) the required memory space.
!
!    Output, real(8) :: RESULT, the estimated value of the integral.
!
!    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
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
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand, in order to
!                             determine the integration difficulties.
!                             if the position of a local difficulty can
!                             be determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used which is designed for
!                             handling the type of difficulty involved.
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
!                             tolerance cannot be achieved due to
!                             roundoff in the extrapolation table, and
!                             that the returned result is the best which
!                             can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier > 0.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             or (flgw /= 1 and flgw /= 2) or
!                             ICALL < 1 or maxp1 < 1.
!                             result, abserr, neval, last, rlist(1),
!                             elist(1), iord(1) and nnlog(1) are set to
!                             zero. alist(1) and blist(1) are set to a
!                             and b respectively.
!
!    Workspace, real(8) :: ALIST(LIMIT), BLIST(LIMIT), contains in entries 1
!    through LAST the left and right ends of the partition subintervals.
!
!    Workspace, real(8) :: RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.
!
!    Workspace, real(8) :: ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.
!
!            iord   - integer
!                     vector of dimension at least limit, the first k
!                     elements of which are pointers to the error
!                     estimates over the subintervals, such that
!                     elist(iord(1)), ..., elist(iord(k)), form
!                     a decreasing sequence, with k = last
!                     if last <= (limit/2+2), and
!                     k = limit+1-last otherwise.
!
!            nnlog  - integer
!                     vector of dimension at least limit, indicating the
!                     subdivision levels of the subintervals, i.e.
!                     iwork(i) = l means that the subinterval numbered
!                     i is of length abs(b-a)*2**(1-l)
!
!         on entry and return
!            momcom - integer
!                     indicating that the Chebyshev moments have been
!                     computed for intervals of lengths
!                     (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
!                     momcom < maxp1
!
!            chebmo - real
!                     array of dimension (maxp1,25) containing the
!                     Chebyshev moments
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
!                       been obtained it is put in rlist2(numrl2) after
!                       numrl2 has been increased by one
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation, i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true value)
!
real(8) :: abseps
NUMFOR_KINDR :: area
NUMFOR_KINDR :: area1
NUMFOR_KINDR :: area12
NUMFOR_KINDR :: area2
real(8) :: a1
real(8) :: a2
real(8) :: b1
real(8) :: b2
real(8) :: correc
real(8) :: defab1
real(8) :: defab2
real(8) :: defabs
real(8) :: domega
real(8) :: dres
real(8) :: erlarg
real(8) :: erlast
real(8) :: errbnd
real(8) :: errmax
real(8) :: error1
real(8) :: erro12
real(8) :: error2
real(8) :: errsum
real(8) :: ertest
logical extall
logical extrap
integer id
integer ierro
integer iroff1
integer iroff2
integer iroff3
integer jupbnd
integer k
integer ksgn
integer ktmin
integer last
integer maxerr
integer nev
logical noext
integer nres
integer nrmax
integer nrmom
integer numrl2
real(8) :: resabs
NUMFOR_KINDR :: reseps
NUMFOR_KINDR :: res3la(3)
NUMFOR_KINDR :: rlist2(52)
real(8) :: small
real(8) :: width
!
!  the dimension of rlist2 is determined by  the value of
!  limexp in QEXTR (rlist2 should be of dimension
!  (limexp+2) at least).
!
#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif
!
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
nnlog(1) = 0

if ((flgw /= 1 .and. flgw /= 2) .or. (epsabs <= 0.0_8 .and. epsrel <= 0.0_8)&
  & .or. icall < 1 .or. maxp1 < 1) then
  ier = 6
  return
end if
!
!  First approximation to the integral.
!
domega = abs(omega)
nrmom = 0

IF (icall <= 1) momcom = 0
call qc25o(f, a, b, domega, flgw, nrmom, maxp1, 0, result, abserr, &
           neval, defabs, resabs, momcom, chebmo)
!
!  Test on accuracy.
!
dres = abs(result)
errbnd = max(epsabs, epsrel * dres)
rlist(1) = result
elist(1) = abserr
iord(1) = 1
IF (abserr <= 1.0e+02 * epsilon(defabs) * defabs .and. abserr > errbnd) ier = 2
IF (limit == 1) ier = 1
IF (ier /= 0 .or. abserr <= errbnd) go to 200

!
!  Initializations
!
errmax = abserr
maxerr = 1
area = result
errsum = abserr
abserr = huge(abserr)
nrmax = 1
extrap = .false.
noext = .false.
ierro = 0
iroff1 = 0
iroff2 = 0
iroff3 = 0
ktmin = 0
small = abs(b - a) * 0.75_8
nres = 0
numrl2 = 0
extall = .False.

if (0.5_8 * abs(b - a) * domega <= 2.0_8) then
  numrl2 = 1
  extall = .True.
  rlist2(1) = result
end if

if (2.5e-01 * abs(b - a) * domega <= 2.0_8) then
  extall = .true.
end if

if (dres >= (1.0_8 - 5.0e+01 * epsilon(defabs)) * defabs) then
  ksgn = 1
else
  ksgn = -1
end if
!
!  main do-loop
!
do last = 2, limit
  !
  !  Bisect the subinterval with the nrmax-th largest error estimate.
  !
  nrmom = nnlog(maxerr) + 1
  a1 = alist(maxerr)
  b1 = 0.5_8 * (alist(maxerr) + blist(maxerr))
  a2 = b1
  b2 = blist(maxerr)
  erlast = errmax
  call qc25o(f, a1, b1, domega, flgw, nrmom, maxp1, 0, area1, &
             error1, nev, resabs, defab1, momcom, chebmo)
  neval = neval + nev

  call qc25o(f, a2, b2, domega, flgw, nrmom, maxp1, 1, area2, &
             error2, nev, resabs, defab2, momcom, chebmo)
  neval = neval + nev
  !
  !  Improve previous approximations to integral and error and
  !  test for accuracy.
  !
  area12 = area1 + area2
  erro12 = error1 + error2
  errsum = errsum + erro12 - errmax
  area = area + area12 - rlist(maxerr)
  if (defab1 == error1 .or. defab2 == error2) go to 25
  if (abs(rlist(maxerr) - area12) > 1.0e-05 * abs(area12) &
      .or. erro12 < 9.9e-01 * errmax) go to 20
  if (extrap) iroff2 = iroff2 + 1

  if (.not. extrap) then
    iroff1 = iroff1 + 1
  end if

20 continue

  IF (last > 10 .and. erro12 > errmax) iroff3 = iroff3 + 1
25 continue

  rlist(maxerr) = area1
  rlist(last) = area2
  nnlog(maxerr) = nrmom
  nnlog(last) = nrmom
  errbnd = max(epsabs, epsrel * abs(area))
  !
  !  Test for roundoff error and eventually set error flag
  !
  if (iroff1 + iroff2 >= 10 .or. iroff3 >= 20) ier = 2

  if (iroff2 >= 5) ierro = 3
  !
  !  Set error flag in the case that the number of subintervals
  !  equals limit.
  !
  if (last == limit) then
    ier = 1
  end if
  !
  !  Set error flag in the case of bad integrand behavior at
  !  a point of the integration range.
  !
  if (max(abs(a1), abs(b2)) <= (1.0_8 + 1.0e+03 * epsilon(a1)) &
      * (abs(a2) + 1.0e+03 * tiny(a2))) then
    ier = 4
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
  !  with nrmax-th largest error estimate (to be bisected next).
  !
  call qsort(limit, last, maxerr, errmax, elist, iord, nrmax)

  IF (errsum <= errbnd) go to 170

  IF (ier /= 0) exit
  IF (last == 2 .and. extall) go to 120
  IF (noext) cycle

  IF (.not. extall) go to 50
  erlarg = erlarg - erlast
  IF (abs(b1 - a1) > small) erlarg = erlarg + erro12
  IF (extrap) go to 70
  !
  !  Test whether the interval to be bisected next is the
  !  smallest interval.
  !
50 continue
  width = abs(blist(maxerr) - alist(maxerr))

  IF (width > small) cycle
  IF (extall) go to 60
  !
  !  Test whether we can start with the extrapolation procedure
  !  (we do this if we integrate over the next interval with
  !  use of a Gauss-Kronrod rule - see QC25O).
  !
  small = small * 0.5_8

  IF (2.5e-01 * width * domega > 2.0_8) cycle

  extall = .True.
  go to 130

60 continue

  extrap = .True.
  nrmax = 2

70 continue

  IF (ierro == 3 .or. erlarg <= ertest) go to 90
  !
  !  The smallest interval has the largest error.
  !  Before bisecting decrease the sum of the errors over the
  !  larger intervals (ERLARG) and perform extrapolation.
  !
  jupbnd = last
  IF (last > (limit / 2 + 2)) jupbnd = limit + 3 - last
  id = nrmax

  do k = id, jupbnd
    maxerr = iord(nrmax)
    errmax = elist(maxerr)
    if (abs(blist(maxerr) - alist(maxerr)) > small) go to 140
    nrmax = nrmax + 1
  end do

  !
  !  Perform extrapolation.
  !
90 continue

  numrl2 = numrl2 + 1
  rlist2(numrl2) = area

  IF (numrl2 < 3) go to 110
  call qextr(numrl2, rlist2, reseps, abseps, res3la, nres)
  ktmin = ktmin + 1

  IF (ktmin > 5 .and. abserr < 1.0e-03 * errsum) ier = 5

  IF (abseps >= abserr) go to 100

  ktmin = 0
  abserr = abseps
  result = reseps
  correc = erlarg
  ertest = max(epsabs, epsrel * abs(reseps))

  IF (abserr <= ertest) exit
  !
  !  Prepare bisection of the smallest interval.
  !
100 continue

  IF (numrl2 == 1) noext = .True.
  IF (ier == 5) exit

110 continue

  maxerr = iord(1)
  errmax = elist(maxerr)
  nrmax = 1
  extrap = .False.
  small = small * 0.5_8
  erlarg = errsum
  cycle

120 continue

  small = small * 0.5_8
  numrl2 = numrl2 + 1
  rlist2(numrl2) = area

130 continue

  ertest = errbnd
  erlarg = errsum
140 continue

end do
!
!  set the final result.
!
IF (abserr == huge(abserr) .or. nres == 0) go to 170

IF (ier + ierro == 0) go to 165
IF (ierro == 3) abserr = abserr + correc
IF (ier == 0) ier = 3
IF (result /= 0.0_8 .and. area /= 0.0_8) go to 160
IF (abserr > errsum) go to 170
IF (area == 0.0_8) go to 190
go to 165

160 continue

IF (abserr / abs(result) > errsum / abs(area)) go to 170
!
!  Test on divergence.
165 continue

if (ksgn == (-1) .and. max(abs(result), abs(area)) <= &
    defabs * 1.0e-2) go to 190

if (1.0e-2 > abs(result / area) .or. abs(result / area) > 100._8 &
    .or. errsum >= abs(area)) ier = 6

go to 190
!
!  Compute global integral sum.
170 continue

result = sum(rlist(1:last))
abserr = errsum

190 continue

IF (ier > 2) ier = ier - 1
200 continue

IF (flgw == 2 .and. omega < 0.0_8) result = -result
end PRNM(qfour)
