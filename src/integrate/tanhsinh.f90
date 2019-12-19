!> @file tanhsinh.f90 provides routines using the tanh-sinh method
!! @date "2019-12-19 08:39:11"

!> Implementation of tanh-sinh integration method.
!! Description: @ref docintegrate
module qtanhsinh
  USE utils, only: dp, Zero, M_PI_2, str, print_msg
  USE func_integ
  implicit none

#ifdef __GFORTRAN__
#define PASTE(a) a
#define JOIN(a,b) PASTE(a)b
#else
#define PASTE(a,b) a ## b
#define JOIN(a,b) PASTE(a,b)
#endif

  !> Subroutine `qnthsh` implements integration by tanh-sinh method.
  !!
  !! This routine is a simple non-adaptive automatic integrator using trapezoidal rule
  !! to integrate a the transformed integrand using the double exponential (also called tanh-sinh) method.
  !! The tahn-sinh scheme is based on the observation that the trapezoid rule converges
  !! very rapidly for functions on the entire real line that go to zero like exp( - exp(t) ).
  !! The change of variables \f$x = \tanh( \pi \sinh(t) /2)\f$ transforms an integral over [-1, 1]
  !! into an integral with integrand suited to the double exponential rule.
  !!
  !! The transformed integral is infinite, but we truncate the domain of integration to [-a, a] (\f$a \le 3\f$).
  !! The value of \f$a\f$ is estimated from the required precision, and the value `3` was chosen for two reasons:
  !! for \f$ t = 3\f$, the transformed \f$x\f$ values are nearly equal to 1 (up to 12 significant figures).
  !! Also, for \f$ t = 3\f$, the smallest weights are 12 orders of magnitude smaller than the largest weights.
  !!
  !! The integration first applies the trapezoid rule to \f$[-a, a]\f$.  Then it subsequently cuts the step size in half each
  !! time, comparing the results.  Integration stops when subsequent iterations are close enough together or the maximum
  !! integration points have been used. By cutting \f$h\f$ in half, the previous integral can be reused; we only need evaluate the
  !! integrand at the newly added points.
  !!
  !! References:
  !! ----------
  !! - https://en.wikipedia.org/wiki/Tanh-sinh_quadrature
  !! - Bailey, David H, ["Tanh-Sinh High-Precision Quadrature", 2006](http://crd.lbl.gov/~dhbailey/dhbpapers/dhb-tanh-sinh.pdf)
  !! - Cook, John D, ["Double Exponential Integration"](https://www.johndcook.com/blog/double_exponential_integration/)
  !!
  !!  @param [in]  f a real or complex function (as defined in [Integrable functions](\ref integwhichfunction))
  !!  @param [in]  a (real(dp)) are the limits of integration
  !!  @param [in]  b (real(dp)) are the limits of integration
  !!  @param [out] IntValue (real or complex, depending on `f`) Estimated value of the integral
  !!  @param [in] epsabs (real, optional)  absolute desired precision. **Default:** epsabs=1.e-7
  !!  @param [in] epsrel (real, optional) relative desired precision. **Default:** epsrel=1.e-5
  !!  @param [out] abserr (real, optional) absolute estimated error
  !!  @param [out] Neval (integer, optional) Number of function evaluations
  !!
  !! Example:
  !! --------
  !! @snippet ex_integr_thsh.f90 Simple
  !!
  interface qnthsh
    module procedure :: d_qnthsh !! real functions
    module procedure :: d_qnthsh_arg !! real functions with extra arguments
    module procedure :: d_qnthsh_f
    module procedure :: c_qnthsh, c_qnthsh_arg, c_qnthsh_f
  end interface qnthsh
  !
  private
  public :: qnthsh

contains
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !              Routines for real(8) functions
#define NUMFOR_KINDR real(dp)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Functions f(x)
#define PRNM(a) subroutine JOIN(d_,a)
#define NUMFOR_EVAL_F(f,x) f(x)
#define NUMFOR_KINDF procedure(funqr)
  ! Quadrature using tanh-sinh scheme
#include "qnthsh.inc"
#undef NUMFOR_EVAL_F
#undef NUMFOR_KINDF
#undef PRNM

  ! Functions f(x, args) via new type
#define PRNM(a) subroutine JOIN(JOIN(d_,a),_arg)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF type(nf_rfunction)
#include "qnthsh.inc"
#undef NUMFOR_EVAL_F
#undef NUMFOR_KINDF
#undef PRNM

  ! Functions f(x, args)
#define PRNM(a) subroutine JOIN(JOIN(d_,a),_f)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF procedure(funqrarg)
#define NF_FARG  type(nf_rfunction)
#include "qnthsh.inc"
#undef NF_FARG
#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef PRNM

#undef NUMFOR_KINDR

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !              Routines for complex functions
#define NUMFOR_KINDR complex(dp)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Functions f(x)
#define PRNM(a) subroutine JOIN(c_,a)
#define NUMFOR_EVAL_F(f,x) f(x)
#define NUMFOR_KINDF procedure(funqc)
#include "qnthsh.inc"
#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef PRNM

  ! Functions f(x, args) via new type
#define PRNM(a) subroutine JOIN(JOIN(c_,a),_arg)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF type(nf_cfunction)
#include "qnthsh.inc"
#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef PRNM

  ! Functions f(x, args)
#define PRNM(a) subroutine JOIN(JOIN(c_,a),_f)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF procedure(funqcarg)
#define NF_FARG  type(nf_cfunction)
#include "qnthsh.inc"
#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef NF_FARG
#undef PRNM

#undef NUMFOR_KINDR

end module qtanhsinh
