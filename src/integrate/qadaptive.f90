!> @file qadaptive.f90 provides routines for integration with adaptive schemes
!! @date "2020-01-09 12:47:37 fiol"

!> @ingroup integrate
!> Globally adaptive Simpson integrator.
!! Description: @ref docintegrate
!!
!! Routines for Globally adaptive integration via Simpson method.
module qadaptive
  USE utils, only: dp, Zero, SMALL, str, print_msg
  USE func_integ
  implicit none

#define DEFAULT_EPSABS 1.e-6_dp
#define DEFAULT_EPSREL 1.e-4_dp

#ifdef __GFORTRAN__
#define PASTE(a) a
#define JOIN(a,b) PASTE(a)b
#else
#define PASTE(a,b) a ## b
#define JOIN(a,b) PASTE(a,b)
#endif

  !> Integration by Adaptive Simpson method of a function
  !! by a globally adaptive strategy, using a Simpson rule.
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
  !! @param [out] abserr (real, **optional**)   Estimation of absolute error achieved
  !! @param [out] neval (integer, **optional**)  Number of function evaluations performed
  !! @param [out] ier (integer, **optional**)  Error code
  interface iads
    module procedure :: d_iads, d_iads_arg, d_iads_f
    module procedure :: c_iads, c_iads_arg, c_iads_f
  end interface iads

  !> Integration by Adaptive Simpson method of a function on a semi-infinite interval,
  !! based on iads()
  !!
  !! The routine calculates an approximation \f$J\f$ to a definite integral
  !! \f[ J \approx I =\int_{a}^{\infty} f(x, args) dx \f]
  !! hopefully satisfying
  !!   \f[ || I - J || \le \max ( epsabs, epsrel \cdot ||I|| ). \f]
  !!
  !! @param [in] f  The function to integrate
  !! @param [in] a (real) lower limit of integration
  !! @param [in] brkpts (real, array) Break points where the integration domain will be split
  !! @param [in] args (real, array, **optional**) extra arguments (if needed) to be passed to the function `f`
  !! @param [out] IntVal (same kind as `f`) Approximation to integral
  !! @param [in] epsabs (real, **optional**)    Absolute accuracy requested. Default = 1.e-7
  !! @param [in] epsrel (real, **optional**)    Relative accuracy requested. Default = 1.e-5
  !! @param [out] abserr (real, **optional**)   Estimation of absolute error achieved
  !! @param [out] neval (integer, **optional**)  Number of function evaluations performed
  !! @param [out] ier (integer, **optional**)  Error code
  interface iadsi
    module procedure :: d_iadsi, d_iadsi_arg, d_iadsi_f
    module procedure :: c_iadsi, c_iadsi_arg, c_iadsi_f
  end interface iadsi
  !
  private
  public :: iads, iadsi

contains
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !              Routines for real(8) functions
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define NUMFOR_KINDR real(dp)

#define PRNM(a) subroutine JOIN(d_,a)

  ! Functions f(x)
#define NUMFOR_EVAL_F(f,x) f(x)
#define NUMFOR_KINDF procedure(funqr)
  ! Quadrature using adaptive-simpson scheme
#include "iads1d.inc"
#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef PRNM

  ! Functions f(x, args) via new type
#define PRNM(a) subroutine JOIN(JOIN(d_,a),_arg)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF type(nf_rfunction)
  ! Quadrature using adaptive-simpson scheme
#include "iads1d.inc"
#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef PRNM

  ! Functions f(x, args)
#define PRNM(a) subroutine JOIN(JOIN(d_,a),_f)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF procedure(funqrarg)
#define NF_FARG  type(nf_rfunction)
  ! Quadrature using adaptive-simpson scheme
#include "iads1d.inc"
#undef NF_FARG
#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef PRNM

#undef NUMFOR_KINDR

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !              Routines for complex functions
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define NUMFOR_KINDR complex(dp)

#define PRNM(a) subroutine JOIN(c_,a)

  ! Functions f(x)
#define NUMFOR_EVAL_F(f,x) f(x)
#define NUMFOR_KINDF procedure(funqc)
  ! Quadrature using adaptive-simpson scheme
#include "iads1d.inc"
#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef PRNM

  ! Functions f(x, args) via new type
#define PRNM(a) subroutine JOIN(JOIN(c_,a),_arg)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF type(nf_cfunction)
  ! Quadrature using adaptive-simpson scheme
#include "iads1d.inc"
#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef PRNM

  ! Functions f(x, args)
#define PRNM(a) subroutine JOIN(JOIN(c_,a),_f)
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)
#define NUMFOR_KINDF procedure(funqcarg)
#define NF_FARG  type(nf_cfunction)
  ! Quadrature using adaptive-simpson scheme
#include "iads1d.inc"
#undef NUMFOR_KINDF
#undef NUMFOR_EVAL_F
#undef NF_FARG
#undef PRNM

#undef NUMFOR_KINDR

end module qadaptive

