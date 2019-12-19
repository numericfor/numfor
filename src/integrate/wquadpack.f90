!> @file wquadpack.f90 provides wrapper routines for integration with QUADPACK
!! @date "2019-12-19 08:39:35"

!> Wrapper of (slightly modified) QUADPACK routines
!! Description: @ref docintegrate
!!
!! ## From QUADPACK original paper:
!!
!! ### Introduction
!!
!! QUADPACK is a Fortran subroutine package for the numerical computation of definite one-dimensional integrals.
!! It originated from a joint project of R. Piessens and E. de Doncker
!! (Appl.  Math. and Progr. Div.- K.U.Leuven, Belgium),
!! C. Ueberhuber (Inst.  Fuer Math.-Techn.U.Wien, Austria),
!! and D. Kahaner (Nation. Bureau of Standards- Washington D.C., U.S.A.).
!!
!! All routines calculate approximations \f$J\f$ to integrals
!! \f[ J \approx I := \int_{a}^{b} f(x, args) dx \f]
!! or integrals of weighted functions:
!! \f[ J \approx I := \int_{a}^{b} f(x) W(x) dx \f]
!! hopefully satisfying error conditions given by the user (epsabs, epsrel)
!!   \f[ || I - J || \le \max ( \mathrm{epsabs}, \mathrm{epsrel} \cdot ||I|| ). \f]
!!
!! ### Guidelines for the use of QUADPACK
!!
!!  Here it is not our purpose to investigate the question when
!!  automatic quadrature should be used. we shall rather attempt
!!  to help the user who already made the decision to use quadpack,
!!  with selecting an appropriate routine or a combination of
!!  several routines for handling his problem.
!!
!!  For both quadrature over finite and over infinite intervals,
!!  one of the first questions to be answered by the user is
!!  related to the amount of computer time he wants to spend,
!!  versus his -own- time which would be needed, for example, for
!!  manual subdivision of the interval or other analytic
!!  manipulations.
!!
!!  1. The user may not care about computer time, or not be
!!      willing to do any analysis of the problem. especially when
!!      only one or a few integrals must be calculated, this attitude
!!      can be perfectly reasonable. in this case it is clear that
!!      either the most sophisticated of the routines for finite
!!      intervals, qags, must be used, or its analogue for infinite
!!      intervals, qagi. these routines are able to cope with
!!      rather difficult, even with improper integrals.
!!      this way of proceeding may be expensive. but the integrator
!!      is supposed to give you an answer in return, with additional
!!      information in the case of a failure, through its error
!!      estimate and flag. yet it must be stressed that the programs
!!      cannot be totally reliable.
!!
!!  2. The user may want to examine the integrand function.
!!      if bad local difficulties occur, such as a discontinuity, a
!!      singularity, derivative singularity or high peak at one or
!!      more points within the interval, the first advice is to
!!      split up the interval at these points. the integrand must
!!      then be examinated over each of the subintervals separately,
!!      so that a suitable integrator can be selected for each of
!!      them. if this yields problems involving relative accuracies
!!      to be imposed on -finite- subintervals, one can make use of
!!      qagp, which must be provided with the positions of the local
!!      difficulties. however, if strong singularities are present
!!      and a high accuracy is requested, application of qags on the
!!      subintervals may yield a better result.
!!
!!    for quadrature over finite intervals we thus dispose of qags
!!    and
!!    - qng for well-behaved integrands,
!!    - qag for functions with an oscillating behavior of a non
!!      specific type,
!!    - qawo for functions, eventually singular, containing a
!!      factor cos(omega*x) or sin(omega*x) where omega is known,
!!    - qaws for integrands with algebraico-logarithmic end point
!!      singularities of known type,
!!    - qawc for cauchy principal values.
!!
!!    @remark
!!      On return, the info type `[d|c]_qp_extra` in the last argument
!!    adaptive integrators contain information about the interval
!!    subdivision process and hence about the integrand behavior:
!!    the end points of the subintervals, the local integral
!!    contributions and error estimates, and eventually other
!!    characteristics. for this reason, and because of its simple
!!    globally adaptive nature, the routine qag in particular is
!!    well-suited for integrand examination. difficult spots can
!!    be located by investigating the error estimates on the
!!    subintervals.
!!
!!    for infinite intervals we provide only one general-purpose
!!    routine, qagi. it is based on the qags algorithm applied
!!    after a transformation of the original interval into (0,1).
!!    yet it may eventuate that another type of transformation is
!!    more appropriate, or one might prefer to break up the
!!    original interval and use qagi only on the infinite part
!!    and so on. these kinds of actions suggest a combined use of
!!    different quadpack integrators. note that, when the only
!!    difficulty is an integrand singularity at the finite
!!    integration limit, it will in general not be necessary to
!!    break up the interval, as qagi deals with several types of
!!    singularity at the boundary point of the integration range.
!!    it also handles slowly convergent improper integrals, on
!!    the condition that the integrand does not oscillate over
!!    the entire infinite interval. if it does we would advise
!!    to sum succeeding positive and negative contributions to
!!    the integral -e.g. integrate between the zeros- with one
!!    or more of the finite-range integrators, and apply
!!    convergence acceleration eventually by means of quadpack
!!    subroutine qext which implements the epsilon algorithm.
!!    such quadrature problems include the fourier transform as
!!    a special case. yet for the latter we have an automatic
!!    integrator available, qawf.
!!
!! ## How to decide what routine to use.
!!
!! ### If your integration region is finite:
!!
!! - If you can factor the integrand as \f$ f(x)=W(x)g(x)\f$, where \f$g(x) \f$ is smooth on [A,B] and \f$ W(x)=\cos(\omega x)\f$ or \f$ W(x)= \sin(\omega x) \f$ then use [qawo()](@ref quadpack::qawo).
!! - Otherwise, if you can factor \f$f(x)=W(x) g(x)\f$ where \f$g(x)\f$ is smooth and
!! \f$ W(x)=(x-a)^{\alpha} (b-x)^{\beta}  (\log(x-a))^L  (\log(b-x))^K  \f$
!! with \f$K, L = 0\f$ or 1, and \f$ \alpha, \beta > -1 \f$, then use  [qaws()](@ref quadpack::qaws).
!! - Otherwise, if you can factor \f$ f(x)=W(x)g(x)\f$, where \f$g(x) \f$ is smooth and \f$ W(x)=1/(x-c)\f$ for some constant c, use [qawc()](@ref quadpack::qawc).
!! - Otherwise, if you do not care too much about possible inefficient use of computer time, and do not want to further analyze the problem, use [qags()](@ref quadpack::qags).
!! - Otherwise, if the integrand is smooth, use [qng()](@ref quadpack::qng) or [qag()](@ref quadpack::qag).
!! - Otherwise, if there are discontinuities or singularities of the integrand or of its derivative, and you know where they are, split the integration range at these points and analyze each subinterval. You can also use `qagp`, which is to be provided with the x-locations of the singularities or discontinuities.
!! - Otherwise, if the integrand has end point singularities, use `qags`.
!! - Otherwise, if the integrand has an oscillatory behavior of nonspecific type, and no singularities, use [qag()](@ref quadpack::qag) with rule="qk61".
!! - Otherwise, use  [qags()](@ref quadpack::qags).
!!
!! ### Routines for an infinite region:
!!
!! - If the integrand decays rapidly to zero, truncate the interval and use the finite interval decision tree.
!! - Otherwise, if the integrand oscillates over the entire infinite range, and the integral is a Fourier transform, use `qawf`.
!! - Or, if the integrand oscillates over the entire infinite range, but is not a Fourier transform, then sum the successive positive and negative contributions by integrating between the zeroes of the integrand. Apply convergence acceleration with `qextr`.
!! - Otherwise, if you are not constrained by computer time, and do not wish to analyze the problem further, use `qagi`.
!! - Otherwise, if the integrand has a non-smooth behavior in the range, and you know where it occurs, split off these regions and use the appropriate finite range routines to integrate over them. Then begin this tree again to handle the remainder of the region.
!! - Otherwise, truncation of the interval, or application of a suitable transformation for reducing the problem to a finite range may be possible. And you may also call `qagi`.
!!
!! ## @anchor integrable Integrable functions
!! All routines described here may integrate functions with one of the following signatures:
!!
!! * real function `f(x)` of a real argument `x`
!! * complex function `f(x)` of a real argument `x`
!! * real function `f(x,args)` of real arguments `x`, and array `args`
!! * complex function `f(x,args)` of real arguments `x`, and array `args`
!!
!! For functions that need extra arguments besides the integration variable `x`, the argument array must be passed to the integration routines.
!!
!! There is also the possibility of using
!! `type(nf_rfunction)` for real functions or
!! `type(nf_cfunction)` for complex functions
!!
module quadpack
  USE utils, only: dp, M_PI, str, str2i, is_inf, nf_minf, nf_inf, print_msg
  USE func_integ
  implicit none

#define DEFAULT_EPSABS 1.e-7_dp
#define DEFAULT_EPSREL 1.e-5_dp

#ifdef __GFORTRAN__
#define PASTE(a) a
#define JOIN(a,b) PASTE(a)b
#else
#define PASTE(a,b) a ## b
#define JOIN(a,b) PASTE(a,b)
#endif

  ! -------- Quadpack ----------------------------------------------------------
  !> qng  estimates an integral using non-adaptive integration.
  !!
  !! The routine calculates an approximation \f$J\f$ to a definite integral
  !! \f[ J \approx I =\int_{a}^{b} f(x, args) dx \f]
  !! hopefully satisfying
  !!   \f[ || I - J || \le \max ( epsabs, epsrel \cdot ||I|| ). \f]
  !!
  !! @param [in] f The function to integrate
  !! @param [in] a (real) lower limit of integration
  !! @param [in] b (real) upper limit of integration
  !! @param [in] args (real, array, **optional**) extra arguments (if needed) to be passed to the function `f`
  !! @param [out] IntVal (same kind as `f`) Approximation to integral
  !! @param [in] epsabs (real, **optional**)    Absolute accuracy requested. Default = 1.e-7
  !! @param [in] epsrel (real, **optional**)    Relative accuracy requested. Default = 1.e-5
  !! @param [out] abserr (real, **optional**)   Estimation of absolute error achieved
  !! @param [out] neval (integer, **optional**)  Number of function evaluations performed
  !! @param [out] ier (integer **optional**)  Error code
  !!
  !! The routine is a simple non-adaptive automatic integrator, based on a sequence of rules with
  !! increasing degree of algebraic precision (Patterson, 1968).  It applies the Gauss-Kronrod
  !! 10-point, 21-point, 43-point and 87-point integration rules in succession until an estimate of
  !! the integral is achieved within the desired absolute and relative error limits.
  !!
  !! Example:
  !! --------
  !! @snippet ex_integr_qng.f90 Simple
  !!
  interface qng
    module procedure :: d_qng, d_qng_arg, d_qng_f
    module procedure :: c_qng, c_qng_arg, c_qng_f
  end interface qng

  !> Globally adaptive integration routine without weights
  !!
  !! The routine calculates an approximation \f$J\f$ to a definite integral
  !! \f[ J \approx I =\int_{a}^{b} f(x, args) dx \f]
  !! hopefully satisfying
  !!   \f[ || I - J || \le \max ( epsabs, epsrel \cdot ||I|| ). \f]
  !!
  !! @param [in] f  The function to integrate
  !! @param [in] a (real) lower limit of integration
  !! @param [in] b (real) upper limit of integration
  !! @param [in] args (real, array, **optional**) extra arguments (if needed) to be passed to the function `f`
  !! @param [out] IntVal (same kind as `f`) Approximation to integral
  !! @param [in] epsabs (real, **optional**)    Absolute accuracy requested. Default = 1.e-7
  !! @param [in] epsrel (real, **optional**)    Relative accuracy requested. Default = 1.e-5
  !! @param [in] rule (character, **optional**)  Choice of integration rule. Default = 'qk21'. Options are:
  !!  - "qk15" ->  7 Gauss points, 15 Gauss-Kronrod points,
  !!  - "qk21" -> 10 Gauss points, 21 Gauss-Kronrod points,
  !!  - "qk31" -> 15 Gauss points, 31 Gauss-Kronrod points,
  !!  - "qk41" -> 20 Gauss points, 41 Gauss-Kronrod points,
  !!  - "qk51" -> 25 Gauss points, 51 Gauss-Kronrod points,
  !!  - "qk61" -> 30 Gauss points, 61 Gauss-Kronrod points.
  !! @param [out] abserr (real, **optional**)   Estimation of absolute error achieved
  !! @param [out] neval (integer, **optional**)  Number of function evaluations performed
  !! @param [out] ier (integer, **optional**)  Error code
  !! @param [in,out] info (**optional**)  Information and workspace.
  !! Must be of type [d_qp_extra](@ref quadpack::d_qp_extra) for integration of real functions and
  !! of type [c_qp_extra](@ref quadpack::c_qp_extra) for integration of complex functions.
  !!
  !! This is a simple globally adaptive integrator using the strategy of Aind (Piessens,
  !! 1973).  The QAG algorithm is a simple adaptive integration procedure. The integration region is
  !! divided into subintervals, and on each iteration the subinterval with the largest estimated error
  !! is bisected. This reduces the overall error rapidly, as the subintervals become concentrated
  !! around local difficulties in the integrand.  It is possible to choose between 6 pairs of
  !! Gauss-Kronrod quadrature formulae for the rule evaluation component.  The pairs of high degree of
  !! precision are suitable for handling integration difficulties due to a strongly oscillating
  !! integrand.
  !!
  !! @note If present, `info` will store and use information about subdivisions, limits, errors, and status messages.
  !!
  !! Example:
  !! --------
  !! @snippet ex_integr_qag.f90 Simple
  !!
  interface qag
    module procedure :: d_qag, d_qag_arg, d_qag_f
    module procedure :: d_qage, d_qage_arg
    module procedure :: c_qag, c_qag_arg, c_qag_f
    module procedure :: c_qage, c_qage_arg
  end interface qag

  !> Subroutine qags is a globally adaptive, automatic interval subdivisions with epsilon extrapolation without weights
  !!
  !! The routine calculates an approximation \f$J\f$ to a definite integral applying adaptively the Gauss-Kronrod 21-point integration rule.
  !!
  !! @param [in] f  The function to integrate
  !! @param [in] a (real) lower limit of integration
  !! @param [in] b (real) upper limit of integration
  !! @param [in] args (real, array, **optional**) extra arguments (if needed) to be passed to the function `f`
  !! @param [out] IntVal (same kind as `f`) Approximation to integral
  !! @param [in] epsabs (real, **optional**)    Absolute accuracy requested. Default = 1.e-7
  !! @param [in] epsrel (real, **optional**)    Relative accuracy requested. Default = 1.e-5
  !! @param [in] gkrule (char(len=4), **optional**)  Choice of integration rule
  !! Possible values are: 'qk15', 'qk21', 'qk31', 'qk41', 'qk51', 'qk61'. Default = 'qk21'.
  !! @param [out] abserr (real, **optional**)   Estimation of absolute error achieved
  !! @param [out] neval (integer, **optional**)  Number of function evaluations performed
  !! @param [out] ier (integer, **optional**)  Error code
  !! @param [in,out] info (**optional**)  Information and workspace.
  !! Must be of type [d_qp_extra](@ref quadpack::d_qp_extra) for integration of real functions and
  !! of type [c_qp_extra](@ref quadpack::c_qp_extra) for integration of complex functions
  !!
  !! @remarks (from
  !! [gsl](https://www.gnu.org/software/gsl/doc/html/integration.html#qags-adaptive-integration-with-singularities))
  !! The presence of an integrable singularity in the integration region causes an adaptive routine to
  !! concentrate new subintervals around the singularity. As the subintervals decrease in size the
  !! successive approximations to the integral converge in a limiting fashion. This approach to the
  !! limit can be accelerated using an extrapolation procedure. The QAGS algorithm combines adaptive
  !! bisection with the Wynn epsilon-algorithm to speed-up the integration of many types of integrable
  !! singularities.
  !!
  !! Either or both limits of integration may be infinite.
  !! By default it uses a Gauss-Kronrod rule of 21 points for finite limits.
  !!
  !! For infinite limits use [nf_minf](@ref basic::nf_minf) for "minus infinite" and [nf_inf](@ref basic::nf_inf) for "plus infinite". A 15-points Gauss-Kronrod rule will be used by default but it will be overriden by the argument `gkrule`.
  !!
  !! The integration over the semi-infinite interval \f$(a,+\infty)\f$ is performed by means of the mapping
  !! \f$x = a + (1-t)/t\f$, onto the semi-open interval \f$(0,1]\f$.
  !!
  !! \f[
  !! \int_{a}^{+\infty} dx f(x) = \int_0^1 dt f(a + (1-t)/t)/t^2
  !! \f]
  !!
  !! The integration over the semi-infinite interval \f$(-\infty, b)\f$ is performed by means of the mapping
  !! \f$x = b - (1-t)/t\f$, onto the semi-open interval \f$(0,1]\f$.
  !!
  !! \f[
  !! \int_{-\infty}^{b} dx f(x) = \int_0^1 dt f(b - (1-t)/t)/t^2
  !! \f]
  !!
  !! The integral over the infinite interval \f$(-\infty,+\infty)\f$ is performed combining the two methods before, i.e: mapping \f$x = (1-t)/t\f$ onto the semi-open interval \f$(0,1]\f$
  !! \f[\int_{-\infty}^{+\infty}dx f(x)=\int_0^1 dt\Big(f((1 - t) / t) + f(-(1 - t) / t) \Big) / t^2.\f]
  !!
  !! Example:
  !! --------
  !! @snippet ex_integr_qags.f90 Simple
  !!
  interface qags
    module procedure :: d_qags, d_qags_arg, d_qags_f
    module procedure :: d_qagse, d_qagse_arg
    module procedure :: c_qags, c_qags_arg, c_qags_f
    module procedure :: c_qagse, c_qagse_arg
  end interface qags

  !> Subroutine qagi is a globally adaptive, automatic interval subdivisions with epsilon
  !> extrapolation without weights for an **infinite** interval
  !!
  !! \copydetails qags
  !!
  !! Either or both limits of integration may be infinite.
  !! For infinite limits use [nf_minf](@ref basic::nf_minf) for "minus infinite" and [nf_inf](@ref basic::nf_inf) for "plus infinite"
  !! @remarks The integral is mapped onto the semi-open interval (0,1] using the transformation \f$x = (1-t)/t\f$ using by default a 15-points Gauss Kronrod rule,
  ! interface qagi
  !   module procedure :: d_qagi, d_qagi_arg, d_qagi_f
  !   module procedure :: d_qagie, d_qagie_arg
  !   module procedure :: c_qagi, c_qagi_arg, c_qagi_f
  !   module procedure :: c_qagie, c_qagie_arg
  ! end interface qagi

  !> Subroutine qagp is a globally adaptive, automatic interval subdivisions with epsilon extrapolation
  !! when a known number of points present singularities. Uses a 21-point Gauss-Kronrod rule
  !!
  !! @param [in] f  The function to integrate
  !! @param [in] a (real) lower limit of integration
  !! @param [in] b (real) upper limit of integration
  !! @param [in] points (real, array) Break points where the integrand may present difficulties, such as singularities or discontinuities.
  !! @param [in] args (real, array, **optional**) extra arguments (if needed) to be passed to the function `f`
  !! @param [out] IntVal (same kind as `f`) Approximation to integral
  !! @param [in] epsabs (real, **optional**)    Absolute accuracy requested. Default = 1.e-7
  !! @param [in] epsrel (real, **optional**)    Relative accuracy requested. Default = 1.e-5
  !! @param [out] abserr (real, **optional**)   Estimation of absolute error achieved
  !! @param [out] neval (integer, **optional**)  Number of function evaluations performed
  !! @param [out] ier (integer, **optional**)  Error code
  !! @param [in,out] info (**optional**)  Information and workspace.
  !! Must be of type [d_qp_extra](@ref quadpack::d_qp_extra) for integration of real functions and
  !! of type [c_qp_extra](@ref quadpack::c_qp_extra) for integration of complex functions
  !!
  interface qagp
    module procedure :: d_qagp, d_qagp_arg, d_qagp_f
    module procedure :: c_qagp, c_qagp_arg, c_qagp_f
    module procedure :: d_qagpe, d_qagpe_arg
    module procedure :: c_qagpe, c_qagpe_arg
  end interface qagp

  interface qfour
    module procedure :: d_qfour, d_qfour_arg, d_qfour_f
    module procedure :: c_qfour, c_qfour_arg, c_qfour_f
  end interface qfour

  !> Subroutine qawo is designed for integrands with an oscillatory factor,
  !! \f[ I= \int_{a}^{b} f(x) W(x) dx \f]
  !! where
  !! \f$ W(x) = \sin(\omega x)\f$ or \f$ W(x) = \cos(\omega x)\f$.
  !!
  !! @param [in] f  The function to integrate
  !! @param [in] a (real) lower limit of integration
  !! @param [in] b (real) upper limit of integration
  !! @param [in] omega (real) factor in the weight function
  !! @param [in] flgw (integer) flag indicating if weight is cosine (flgw=1)
  !! @param [in] args (real, array, **optional**) extra arguments (if needed) to be passed to the function `f`
  !! @param [out] IntVal (same kind as `f`) Approximation to integral
  !! @param [in] epsabs (real, **optional**)    Absolute accuracy requested. Default = 1.e-7
  !! @param [in] epsrel (real, **optional**)    Relative accuracy requested. Default = 1.e-5
  !! @param [out] abserr (real, **optional**)   Estimation of absolute error achieved
  !! @param [out] neval (integer, **optional**)  Number of function evaluations performed
  !! @param [out] ier (integer, **optional**)  Error code
  !! @param [in,out] info (**optional**)  Information and workspace.
  !! Must be of type [d_qp_extra](@ref quadpack::d_qp_extra) for integration of real functions and
  !! of type [c_qp_extra](@ref quadpack::c_qp_extra) for integration of complex functions
  !!
  !! @note An automatic subdivision of the interval is preformed and results are extrapolated using the epsilon-algorithm to accelerate the convergence of the integral.
  interface qawo
    module procedure :: d_qawo, d_qawo_arg, d_qawo_f
    module procedure :: c_qawo, c_qawo_arg, c_qawo_f
  end interface qawo

  !> qawf computes Fourier integrals over the interval [ A, +Infinity ).
  !!
  !! \f[ I= \int_{a}^{\infty} f(x) W(x) dx \f]
  !! where
  !! \f$ W(x) = \sin(\omega x)\f$ or \f$ W(x) = \cos(\omega x)\f$.
  !!
  !! @param [in] f  The function to integrate
  !! @param [in] a (real) lower limit of integration
  !! @param [in] omega (real) factor in the weight function
  !! @param [in] flgw (integer) flag indicating if weight is cosine (flgw=1)
  !! @param [in] args (real, array, **optional**) extra arguments (if needed) to be passed to the function `f`
  !! @param [out] IntVal (same kind as `f`) Approximation to integral
  !! @param [in] epsabs (real, **optional**)    Absolute accuracy requested. Default = 1.e-7
  !! @param [out] abserr (real, **optional**)   Estimation of absolute error achieved
  !! @param [out] neval (integer, **optional**)  Number of function evaluations performed
  !! @param [out] ier (integer, **optional**)  Error code
  !! @param [in,out] info (**optional**)  Information and workspace.
  !! Must be of type [d_qp_extra](@ref quadpack::d_qp_extra) for integration of real functions and
  !! of type [c_qp_extra](@ref quadpack::c_qp_extra) for integration of complex functions.
  !!
  !! @note The integral is computed using the `qawo` subroutine over each of the subintervals:
  !! \f{align*}{
  !!     C_1 &= [a,a+c] \\
  !!     C_2 &= [a+c,a+2c] \\
  !!     \vdots & \\
  !!     C_k &= [a+(k-1)c,a+kc]
  !!     \f}
  !! where \f$ c = (2 \, \mathrm{floor}(|\omega|) + 1) \pi/|\omega| \f$ is chosen to span an odd number of periods.
  !!
  !! Example:
  !! --------
  !! @snippet ex_integr_qawf.f90 Simple
  !!
  interface qawf
    module procedure :: d_qawf, d_qawf_arg, d_qawf_f
    module procedure :: d_qawfe, d_qawfe_arg
    module procedure :: c_qawf, c_qawf_arg, c_qawf_f
    module procedure :: c_qawfe, c_qawfe_arg
  end interface qawf

  !> Subroutine qaws estimates integrals with algebraico-logarithmic endpoint singularities given by the weight function \f$ W(x) \f$
  !! \f[ I= \int_{a}^{b} f(x) W(x) dx \f]
  !!
  !! @param [in] f  The function to integrate
  !! @param [in] a (real) lower limit of integration
  !! @param [in] b (real) upper limit of integration
  !! @param [in] alfa (real) parameter of the weight function
  !! @param [in] beta (real) parameter of the weight function
  !! @param [in] flgw (integer) flag indicating the kind of weight function used:
  !! - flgw = 1  \f$ \Rightarrow W(x) = (x-a)^\alpha (b-x)^\beta \f$
  !! - flgw = 2  \f$ \Rightarrow W(x) = (x-a)^\alpha (b-x)^\beta \log(x-a) \f$
  !! - flgw = 3  \f$ \Rightarrow W(x) = (x-a)^\alpha (b-x)^\beta \log(b-x) \f$
  !! - flgw = 4  \f$ \Rightarrow W(x) = (x-a)^\alpha (b-x)^\beta \log(x-a) \log(b-x) \f$
  !! @param [in] args (real, array, **optional**) extra arguments (if needed) to be passed to the function `f`
  !! @param [out] IntVal (same kind as `f`) Approximation to integral
  !! @param [in] epsabs (real, **optional**)    Absolute accuracy requested. Default = 1.e-7
  !! @param [in] epsrel (real, **optional**)    Relative accuracy requested. Default = 1.e-5
  !! @param [out] abserr (real, **optional**)   Estimation of absolute error achieved
  !! @param [out] neval (integer, **optional**)  Number of function evaluations performed
  !! @param [out] ier (integer, **optional**)  Error code
  !! @param [in,out] info (**optional**)  Information and workspace.
  !! Must be of type [d_qp_extra](@ref quadpack::d_qp_extra) for integration of real functions and
  !! of type [c_qp_extra](@ref quadpack::c_qp_extra) for integration of complex functions.
  !!
  !! Example:
  !! --------
  !! @snippet ex_integr_qaws.f90 Simple
  !!
  interface qaws
    module procedure :: d_qaws, d_qaws_arg, d_qaws_f
    module procedure :: d_qawse, d_qawse_arg
    module procedure :: c_qaws, c_qaws_arg, c_qaws_f
    module procedure :: c_qawse, c_qawse_arg
  end interface qaws

  !> Subroutine qawc computes the Cauchy principal value.
  !!
  !! \f[ I= \int_{a}^{b} f(x) W(x) dx \f]
  !! where
  !! \f$ W(x) = \frac{1}{x-c} \f$.
  !! @param [in] f  The function to integrate
  !! @param [in] a (real) lower limit of integration
  !! @param [in] b (real) upper limit of integration
  !! @param [in] c (real) Point where lies the singularity
  !! @param [in] args (real, array, **optional**) extra arguments (if needed) to be passed to the function `f`
  !! @param [out] IntVal (same kind as `f`) Approximation to integral
  !! @param [in] epsabs (real, **optional**)    Absolute accuracy requested. Default = 1.e-7
  !! @param [in] epsrel (real, **optional**)    Relative accuracy requested. Default = 1.e-5
  !! @param [out] abserr (real, **optional**)   Estimation of absolute error achieved
  !! @param [out] neval (integer, **optional**)  Number of function evaluations performed
  !! @param [out] ier (integer, **optional**)  Error code
  !! @param [in,out] info (**optional**)  Information and workspace.
  !! Must be of type [d_qp_extra](@ref quadpack::d_qp_extra) for integration of real functions and
  !! of type [c_qp_extra](@ref quadpack::c_qp_extra) for integration of complex functions.
  !!
  !!@note As explained in
  !! [gsl](https://www.gnu.org/software/gsl/doc/html/integration.html#qags-adaptive-integration-with-singularities),
  !! the adaptive bisection algorithm of QAG is used, with modifications to ensure that subdivisions
  !! do not occur at the singular point \f$x = c\f$.  When a subinterval contains the point \f$x =
  !! c\f$ or is close to it then a special 25-point modified Clenshaw-Curtis rule is used to control
  !! the singularity.  Further away from the singularity the algorithm uses an ordinary 15-point
  !! Gauss-Kronrod integration rule.
  !!
  !! Example:
  !! --------
  !! @snippet ex_integr_qawc.f90 Simple
  !!
  interface qawc
    module procedure :: d_qawc, d_qawc_arg, d_qawc_f
    module procedure :: d_qawce, d_qawce_arg
    module procedure :: c_qawc, c_qawc_arg, c_qawc_f
    module procedure :: c_qawce, c_qawce_arg
  end interface qawc

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!    Non automatic routines    !!!!!!!!!!!!!!!!!!!!
  !> Routine to perform the integration of a function by Gauss-Kronrod rule
  !!
  !! This routine is non-automatic and approximates the integral of the function and its absolute value
  !! @param [in] f  The function to integrate
  !! @param [in] a (real) lower limit of integration
  !! @param [in] b (real) upper limit of integration
  !! @param [in] args (real, array, optional) extra arguments (if needed) to be passed to the function `f`
  !! @param [out] IntVal (same kind as `f`) Approximation to integral I = `integ(f(x), a, b)`, i.e:
  !! \f$ I =\int_{a}^{b} f(x) dx \f$
  !! @param [in] rule (character(len=4), **optional**)  Choice of integration rule. Default = 'qk21'.
  !! Options are:
  !!  - "qk15" ->  7 Gauss points, 15 Gauss-Kronrod points,
  !!  - "qk21" -> 10 Gauss points, 21 Gauss-Kronrod points,
  !!  - "qk31" -> 15 Gauss points, 31 Gauss-Kronrod points,
  !!  - "qk41" -> 20 Gauss points, 41 Gauss-Kronrod points,
  !!  - "qk51" -> 25 Gauss points, 51 Gauss-Kronrod points,
  !!  - "qk61" -> 30 Gauss points, 61 Gauss-Kronrod points.
  !! @param [out] abserr (real) Estimation of error
  !! @param [out] resabs (same kind as `f`) Approximation to integral of absolute value of `f`
  !! \f$ I_1 =\int_{a}^{b} |f(x)| dx \f$
  !! @param [out] resacs (real) Approximation to integral
  !! \f$ I_2 =\int_{a}^{b} |f(x)-I/(b-a)| dx \f$
  !!
  !! The routine returns the result of applying the m-point Kronrod (result I) rule given by arg `rule`
  !! by optimal addition of abscissae to the n-point Gauss rule (result J), where \f$n=(m-1)/2\f$ .\n
  !! The absolute error is evaluated as `abserr`=|I-J|.
  interface qgk
    module procedure :: d_qgk, d_qgk_arg, d_qgk_f
    module procedure :: c_qgk, c_qgk_arg, c_qgk_f
  end interface qgk

  interface qk15
    module procedure :: d_qk15, d_qk15_arg, d_qk15_f
    module procedure :: c_qk15, c_qk15_arg, c_qk15_f
  end interface qk15

  ! Routine to perform the integration of a function by 21-points Gauss-Kronrod rule
  interface qk21
    module procedure :: d_qk21, d_qk21_arg, d_qk21_f
    module procedure :: c_qk21, c_qk21_arg, c_qk21_f
  end interface qk21

  ! Routine to perform the integration of a function by 31-points Gauss-Kronrod rule
  interface qk31
    module procedure :: d_qk31, d_qk31_arg, d_qk31_f
    module procedure :: c_qk31, c_qk31_arg, c_qk31_f
  end interface qk31

  ! Routine to perform the integration of a function by 41-points Gauss-Kronrod rule
  interface qk41
    module procedure :: d_qk41, d_qk41_arg, d_qk41_f
    module procedure :: c_qk41, c_qk41_arg, c_qk41_f
  end interface qk41
  ! Routine to perform the integration of a function by 51-points Gauss-Kronrod rule
  interface qk51
    module procedure :: d_qk51, d_qk51_arg, d_qk51_f
    module procedure :: c_qk51, c_qk51_arg, c_qk51_f
  end interface qk51

  ! Routine to perform the integration of a function by 61-points Gauss-Kronrod rule
  interface qk61
    module procedure :: d_qk61, d_qk61_arg, d_qk61_f
    module procedure :: c_qk61, c_qk61_arg, c_qk61_f
  end interface qk61

  !> Routine to perform the integration of a weighted function by 15-points Gauss-Kronrod rule
  !!
  !! This routine is non-automatic and approximates the integral of the function and its absolute value
  !! \f[ I= \int_{a}^{b} f(x) W(x) dx \f]
  !! @param f (in) The function to integrate
  !! @param w (in) The (real) weight function `W(x)`
  !! @param p (in, real, array) extra arguments for the weight function
  !! @param kp (in, integer) flag indicating the type of weight function
  !! @param a (in, real) lower limit of integration
  !! @param b (in, real) upper limit of integration
  !! @param args (in, real, array, optional) extra arguments (if needed) to be passed to the function `f`
  !! @param result (out, same kind as `f`) Approximation to integral I = `integ(f(x), a, b)`, i.e:
  !! \f$ I =\int_{a}^{b} W(x) f(x) dx \f$
  !! @param abserr (out, real) Estimation of error
  !! @param resabs (out, same kind as `f`) Approximation to integral of absolute value of `f`
  !! \f$ I_1 =\int_{a}^{b} W(x) |f(x)| dx \f$
  !! @param resacs (out, real) Approximation to integral
  !! \f$ I_2 =\int_{a}^{b} |W(x) f(x)-I/(b-a)| dx \f$
  !!
  !! The routine returns the result of applying the m-point Kronrod (result I) rule obtained
  !! by optimal addition of abscissae to the n-point Gauss rule (result J), where \f$n=(m-1)/2\f$ .\n
  !! The absolute error is evaluated as `abserr`=|I-J|.
  interface qk15w
    module procedure :: d_qk15w, d_qk15w_arg, d_qk15w_f
    module procedure :: c_qk15w, c_qk15w_arg, c_qk15w_f
  end interface qk15w

  interface qc25o
    module procedure :: d_qc25o, d_qc25o_arg, d_qc25o_f
    module procedure :: c_qc25o, c_qc25o_arg, c_qc25o_f
  end interface qc25o

  interface qc25c
    module procedure :: d_qc25c, d_qc25c_arg, d_qc25c_f
    module procedure :: c_qc25c, c_qc25c_arg, c_qc25c_f
  end interface qc25c

  interface qc25s
    module procedure :: d_qc25s, d_qc25s_arg, d_qc25s_f
    module procedure :: c_qc25s, c_qc25s_arg, c_qc25s_f
  end interface qc25s

  ! ! QUADPACK global interface
  ! !> Subroutine `quad` computes a definite integral of a function *without weights* from `a` to `b`
  ! !! (possibly infinite interval) using a routine from QUADPACK.
  ! !!
  ! !! @param [in] f  The function to integrate (real or complex). See @ref integrable functions
  ! !! @param [in] a (real) lower limit of integration. Use `nf_minf` for \f$-\infty \f$.
  ! !! @param [in] b (real) upper limit of integration.  Use `nf_inf` for \f$+\infty \f$.
  ! !! @param [in] args (real, array, **optional**) extra arguments (if needed) to be passed to the function `f`
  ! !! @param [out] IntVal (same kind as `f`) Approximation to integral
  ! !! @param [out] abserr (real)   Estimation of absolute error achieved
  ! !! @param [in] epsabs (real, **optional**)    Absolute accuracy requested. Default = 1.e-7
  ! !! @param [in] epsrel (real, **optional**)    Relative accuracy requested. Default = 1.e-5
  ! !! @param [in] maxsub (integer, **optional**)  Maximum number of subdivisions allowed
  ! !! @param [in] points (real, array, **optional**) Break points where the integrand may present difficulties, such as singularities or discontinuities.
  ! !! @param [out] neval (integer, **optional**)  Number of function evaluations performed
  ! !! @param [out] ier (integer, **optional**)  Error code
  ! !! @param [in,out] info (**optional**)  Information and workspace.
  ! !! Must be of type [d_qp_extra](@ref quadpack::d_qp_extra) for integration of real functions and
  ! !! of type [c_qp_extra](@ref quadpack::c_qp_extra) for integration of complex functions
  ! !!
  ! interface quad
  !   module procedure :: d_quad, d_quad_arg, d_quad_f
  !   module procedure :: c_quad, c_quad_arg, c_quad_f
  ! end interface quad

  ! Auxiliary routines
  !> Subroutine qcheb computes the coefficients of the Chebyshev series expansion
  !! of degrees 12 and 24 of a function using a fast Fourier transform method
  !! \f[  f(x) = \sum_{k} (\mathrm{Cheb}(k) T_{k-1}(x) , \f]
  !!    where \f$ T_{k}(x)\f$ is the Chebyshev polynomial of degree k.
  !!
  !!  **Reference:**
  !!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
  !!    QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983. (Section 2.2.3.2)
  !! @param [in] x (real, dimension(11)) values of \f$ \cos(k*\pi/24)\f$, for \f$k = 1, \dots, 11 \f$
  !! @param [in,out] fval (real or complex, dimension(25))
  !! @param [out] cheb12 (real or complex) the Chebyshev coefficients for degree 12
  !! @param [out] cheb24 (real or complex) the Chebyshev coefficients for degree 24
  interface qcheb
    module procedure :: d_qcheb, c_qcheb
  end interface qcheb
  !

  !> Subroutine qextr carries out the Epsilon extrapolation algorithm.
  !!
  !! The routine determines the limit of a given sequence of approximations, by means of the epsilon algorithm of P. Wynn.
  !! An estimate of the absolute error is given.  The condensed epsilon table is computed.
  !! Only those elements needed for the computation of the next diagonal are preserved.
  !! @param [in] n (integer) The entry of EPSTAB which contains the new element in the first column of the epsilon table.
  !! @param [in] epstab (real, array) the two lower diagonals of the triangular epsilon table.
  !! The elements are numbered starting at the right-hand corner of the triangle.
  !! @param [out] result the estimated value of the integral
  !! @param [out] res3la (real, array dimension(3)) the last three results
  !! @param [in,out] nres (integer) the number of calls to the routine.  This should be zero on
  !! the first call, and is automatically updated before return.
  interface qextr
    module procedure :: d_qextr, c_qextr
  end interface qextr

  !> Type used to pass extra (optional) information with the integration routines.
  type d_qp_extra
    character(len=:), allocatable :: msg         !< Output status message
    real(dp), dimension(:), allocatable :: alist !< List of left end points of subintervals considered up to now
    real(dp), dimension(:), allocatable :: blist !< List of right end points of subintervals considered up to now
    real(dp), dimension(:), allocatable :: rlist !< Approximation to the integral over (alist(i),blist(i))
    real(dp), dimension(:), allocatable :: elist !< Error estimate applying to rlist
    integer, dimension(:), allocatable :: iord   !< Index of error estimation for the subintervals
    !! elist(iord(1)), ..., elist(iord(k)) with k = last if last <= (limit/2+2), and k = limit+1-last otherwise,
    !! form a decreasing sequence.
    integer :: last                              !< the number of subintervals actually produced
    integer :: size = 0                          !< Maximum number of subdivisions allowed (size of arrays)
  contains
    procedure, pass(y) :: init => d_sinit
    generic :: Assignment(=) => init
  end type d_qp_extra

  interface d_qp_extra
    module procedure :: d_init
  end interface d_qp_extra

  !> Type used to pass extra (optional) information with the integration routines.
  type c_qp_extra
    character(len=:), allocatable :: msg            !< Output status message
    real(dp), dimension(:), allocatable :: alist    !< List of left end points of subintervals considered up to now
    real(dp), dimension(:), allocatable :: blist    !< List of right end points of subintervals considered up to now
    complex(dp), dimension(:), allocatable :: rlist !< Approximation to the integral over (alist(i),blist(i))
    real(dp), dimension(:), allocatable :: elist    !< Error estimate applying to rlist
    integer, dimension(:), allocatable :: iord      !< Index of error estimation for the subintervals
    !! elist(iord(1)), ..., elist(iord(k)) with k = last if last <= (limit/2+2), and k = limit+1-last otherwise,
    !! form a decreasing sequence.
    integer :: last                                 !< the number of subintervals actually produced
    integer :: size = 0                             !< Maximum number of subdivisions allowed (size of arrays)
  contains
    procedure, pass(y) :: init => c_sinit
    generic :: Assignment(=) => init
  end type c_qp_extra
  interface c_qp_extra
    module procedure :: c_init
  end interface c_qp_extra

  private         ! By default everything private
  !
  ! Public interfaces
  public :: d_qp_extra, c_qp_extra
  public :: qp_err_msg
  public :: qgk
  public :: qng ! Non automatic
  public :: qag, qags, qagp ! Generic automatic
  public :: qawo, qawf, qaws, qawc ! Specific automatic

contains

  ! ============================================================================
  !> qp_err_msg Returns the error message correspondig to a given code
  function qp_err_msg(ier) result(S)
    implicit none
    character(len=:), allocatable :: S !< Error message
    integer, intent(IN) :: ier !< Error code
    !! Examples:
    !!
    select case (ier)
    case (0)
      S = "Normal and reliable termination of the routine.&
        &  It is assumed that the requested accuracy has been achieved."
    case (1)
      S = "Maximum number of subdivisions allowed has been achieved. &
        & One can allow more sub-divisions by increasing the data value of limit&
        & (and taking the according dimension adjustments into account) .&
        & However, if this yields no improvement it is advised to analyze the integrand&
        &  in order to determine the integration difficulties.&
        & If the position of a local difficulty can be determined &
        & (e.g.singularity, discontinuity within the interval) one will probably gain from &
        & splitting up the interval at this point and calling the integrator on the sub-ranges.&
        & If possible, an appropriate special-purpose integrator should be used, &
        & which is designed for handling the type of difficulty involved."
    case (2)
      S = "The occurrence of roundoff error is detected, which prevents &
        & the requested tolerance from being achieved. The error may be underestimated."
    case (3)
      S = "Extremely bad integrand behavior occurs at some points of the integration interval."
    case (4)
      S = "The algorithm does not converge. Roundoff error is detected in the extrapolation table.  &
        & It is presumed that the requested tolerance cannot be achieved, &
        & and that the returned result is the best which can be obtained."
    case (5)
      S = "The integral is probably divergent, or slowly convergent.&
        &  It must be noted that divergence can occur with any other value of ier."
    case (6)
      S = "The input is invalid &
        & The result, absolute error and number of function evaluations are set to zero."
    case (7)
      S = "Abnormal termination of the computation of one or more subintegrals"
    case (8)
      S = "Maximum number of cycles allowed has been achieved, &
        & i.e. maximum number of subintervals (a+(k-1)c,a+kc) &
        &where c = (2*int(abs(omega))+1)*pi/abs(omega), for k = 1, 2, ..."
    case (9)
      S = "The extrapolation table constructed for convergence acceleration of the series &
        & formed by the integral contributions over the cycles does not converge&
        & to within the requested accuracy."

    end select
  end function qp_err_msg

  !> init initializes the components
  !!
  subroutine d_sinit(y, length)
    implicit none
    class(d_qp_extra), intent(OUT) :: y
    integer, intent(IN) :: length
    if (allocated(y%alist) .and. length /= y%size) then
      deallocate (y%alist, y%blist, y%elist)
      deallocate (y%rlist)
      deallocate (y%iord)
      y%last = 0
      y%size = 0
    end if
    if (length /= y%size .and. length /= 0) then
      y%last = 0
      y%size = length
      allocate (y%alist(length), y%blist(length), y%elist(length))
      allocate (y%rlist(length))
      allocate (y%iord(length))
    end if
  end subroutine d_sinit
  function d_init(length) result(y)
    implicit none
    type(d_qp_extra) :: y !<
    integer, intent(IN) :: length !<
    if (allocated(y%alist) .and. length /= y%size) then
      deallocate (y%alist, y%blist, y%elist)
      deallocate (y%rlist)
      deallocate (y%iord)
      y%last = 0
      y%size = 0
    end if
    if (length /= y%size .and. length /= 0) then
      y%last = 0
      y%size = length
      allocate (y%alist(length), y%blist(length), y%elist(length))
      allocate (y%rlist(length))
      allocate (y%iord(length))
    end if
  end function d_init

  subroutine c_sinit(y, length)
    implicit none
    class(c_qp_extra), intent(OUT) :: y
    integer, intent(IN) :: length
    if (allocated(y%alist) .and. length /= y%size) then
      deallocate (y%alist, y%blist, y%elist)
      deallocate (y%rlist)
      deallocate (y%iord)
      y%last = 0
      y%size = 0
    end if
    if (length /= y%size .and. length /= 0) then
      y%last = 0
      y%size = length
      allocate (y%alist(length), y%blist(length), y%elist(length))
      allocate (y%rlist(length))
      allocate (y%iord(length))
    end if
  end subroutine c_sinit
  function c_init(length) result(y)
    implicit none
    type(c_qp_extra) :: y !<
    integer, intent(IN) :: length !<
    if (allocated(y%alist) .and. length /= y%size) then
      deallocate (y%alist, y%blist, y%elist)
      deallocate (y%rlist)
      deallocate (y%iord)
      y%last = 0
      y%size = 0
    end if
    if (length /= y%size .and. length /= 0) then
      y%last = 0
      y%size = length
      allocate (y%alist(length), y%blist(length), y%elist(length))
      allocate (y%rlist(length))
      allocate (y%iord(length))
    end if
  end function c_init

  ! ============================================================================
  ! ----------------------------- QUADPACK ROUTINES ----------------------------
  ! ============================================================================
  ! Routines that have only one version.
  ! Auxiliary routines that do not depend on the integrand function (neither signature nor return kind)
#include "./quadpack/qaux1.inc"

  ! ------------------------------------------------------------------------
  !            Routines for integration of real(8) functions
  ! ------------------------------------------------------------------------
#define NUMFOR_KINDR real(8)
#define NUMFOR_QP_INFO type(d_qp_extra)

#define PRNM(a) subroutine JOIN(d_,a)

  ! Auxiliary routines that do not depend on the signature of the integrand function but
  ! depend on their return kind (real or complex)
#include "./quadpack/qaux2.inc"

  ! Functions f(x)
#define NUMFOR_EVAL_F(f,x) f(x)
#define NUMFOR_KINDF procedure(funqr)

  ! Quadrature using Gauss-Kronrod rules
#include "./quadpack/qgk.inc"
  ! Quadrature using Clenshaw-Curtis method
#include "./quadpack/qcc.inc"
  ! Quadpack automatic routines
#include "./quadpack/qp_intege.inc"
#include "./quadpack/qp_integ.inc"

#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef PRNM

  ! Functions f(x, args) via new type
#define PRNM(a) subroutine JOIN(JOIN(d_,a),_arg)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF type(nf_rfunction)

  ! Quadrature using Gauss-Kronrod rules
#include "./quadpack/qgk.inc"
  ! Quadrature using Clenshaw-Curtis method
#include "./quadpack/qcc.inc"
  ! Quadpack automatic routines
#include "./quadpack/qp_intege.inc"
#include "./quadpack/qp_integ.inc"

#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef PRNM

  ! Functions f(x, args)
#define PRNM(a) subroutine JOIN(JOIN(d_,a),_f)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF procedure(funqrarg)
#define NF_FARG  type(nf_rfunction)
  ! Quadrature using Gauss-Kronrod rules
#include "./quadpack/qgk.inc"
  ! Quadrature using Clenshaw-Curtis method
#include "./quadpack/qcc.inc"
  ! Quadpack automatic routines
#include "./quadpack/qp_integ.inc"

#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef NF_FARG
#undef PRNM

#undef NUMFOR_QP_INFO
#undef NUMFOR_KINDR

  ! ------------------------------------------------------------------------
  !          Routines for integration of complex(8) functions
  ! ------------------------------------------------------------------------

#define NUMFOR_KINDR complex(8)
#define NUMFOR_QP_INFO type(c_qp_extra)
#define PRNM(a) subroutine JOIN(c_,a)
  ! Auxiliary routines that do not depend on the signature of the integrand function but
  ! depend on their return kind (real or complex)
#include "./quadpack/qaux2.inc"

  ! Functions f(x)
#define NUMFOR_EVAL_F(f,x) f(x)
#define NUMFOR_KINDF procedure(funqc)

  ! Quadrature using Gauss-Kronrod rules
#include "./quadpack/qgk.inc"
  ! Quadrature using Clenshaw-Curtis method
#include "./quadpack/qcc.inc"
  ! Quadpack automatic routines
#include "./quadpack/qp_intege.inc"
#include "./quadpack/qp_integ.inc"

#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef PRNM

  ! Functions f(x, args) via new type
#define PRNM(a) subroutine JOIN(JOIN(c_,a),_arg)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF type(nf_cfunction)
  ! Quadrature using Gauss-Kronrod rules
#include "./quadpack/qgk.inc"
  ! Quadrature using Clenshaw-Curtis method
#include "./quadpack/qcc.inc"
  ! Quadpack automatic routines
#include "./quadpack/qp_intege.inc"
#include "./quadpack/qp_integ.inc"

#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F

#undef PRNM

  ! Functions f(x, args)
#define PRNM(a) subroutine JOIN(JOIN(c_,a),_f)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF procedure(funqcarg)
#define NF_FARG  type(nf_cfunction)
  !   ! Quadrature using Gauss-Kronrod rules
#include "./quadpack/qgk.inc"
  ! Quadrature using Clenshaw-Curtis method
#include "./quadpack/qcc.inc"
  ! Quadpack automatic routines
#include "./quadpack/qp_integ.inc"

#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef NF_FARG
#undef PRNM

  ! Finishing with complex
#undef NUMFOR_QP_INFO
#undef NUMFOR_KINDR

end module quadpack

