!> @file qsimpson.f90 provides routines for Trapezoide and Simpson integration
!! @date "2024-02-21 15:01:01"

!> @ingroup integrate
!> The qsimpson module provides routines for trapezoid and simpson integration of both sampled data and functions.
!! Description: @ref submodule-integrate
module qsimpson
  USE utils, only: dp, Zero, M_PI_2, SMALL, str, print_msg
  USE grids, only: linspace
  USE func_integ
  implicit none

#ifdef __GFORTRAN__
#define PASTE(a) a
#define JOIN(a,b) PASTE(a)b
#else
#define PASTE(a,b) a ## b
#define JOIN(a,b) PASTE(a,b)
#endif

  !> Routine for integration of sampled values or functions by using the trapezoid rule.
  !!
  !! For integration of a real-valued and complex-valued function `f` the following calls are valid
  !!
  !!```{.f90}
  !! call trapz(f)           ! Using 100 equispaced points between 0 and 1
  !! call trapz(f, a, b)     ! Using 100 equispaced points between a and b
  !! call trapz(f, a, b, N)  ! Using N equispaced points between 0 an 1
  !! call trapz(f, x)        ! Using evaluations on values of array x
  !! ```
  !!
  !! Here:
  !! @param [in] f is a real or complex function
  !! @param [in] a (real) Lower limit of integration
  !! @param [in] b (real) Upper limit of integration
  !! @param [in] N (integer) Number of equally-spaced points to use
  !! @param [in] x (real, array) Points where the function will be evaluated
  !!
  !! For sampled values it accepts several input formats,
  !!
  !! ```{.f90}
  !! call trapz(y)           ! Using values in array y, corresponding to x ∈ [0,1]
  !! call trapz(y, a, b)     ! Using values in array y, corresponding to x ∈ [a,b]
  !! call trapz(y, dh)       ! Using values in y, for x values equispaced in dh
  !! call trapz(y, x)        ! Using values in y, corresponding to x
  !! ```
  !!
  !! Here:
  !! @param [in] y is a real or complex array
  !! @param [in] a (real) Lower limit of integration
  !! @param [in] b (real) Upper limit of integration
  !! @param [in] dh (real) Space between points in the x-axis
  !! @param [in] x (real, array) Points where the function is evaluated
  !!
  interface trapz
    !SAMPLE:
    module procedure d_trapz_sample, c_trapz_sample
    module procedure d_trapz_dh_sample, c_trapz_dh_sample
    module procedure d_trapz_nonlin_sample, c_trapz_nonlin_sample
    !FUNCTION
    ! f(x)
    module procedure d_trapz_func, d_trapz_nonlin_func
    module procedure c_trapz_func, c_trapz_nonlin_func
    ! f(x, args)
    module procedure d_trapz_func_arg, d_trapz_nonlin_func_arg
    module procedure c_trapz_func_arg, c_trapz_nonlin_func_arg
    ! ------------------------------------------------------
    ! 2D
    module procedure d_trapz2d_func
    module procedure c_trapz2d_func
    module procedure d_trapz2d_func_eps
    module procedure c_trapz2d_func_eps
    module procedure d_trapz2d_func_arg
    module procedure c_trapz2d_func_arg
    module procedure d_trapz2d_func_eps_arg
    module procedure c_trapz2d_func_eps_arg
  end interface trapz
  !
  interface simpson_dx_sample
    module procedure d_simpson_dx_sample, c_simpson_dx_sample
  end interface simpson_dx_sample

  !SIMPSON'S RULE
  !> Routines for integration of sampled values or functions by using Simpson rule
  !!
  !! The use of this routine is simlar to the use of `trapz`.
  !! For integration of a real-valued and complex-valued function `f` the following calls are valid
  !!
  !!```{.f90}
  !! call simps(f)           ! Using 100 equispaced points between 0 and 1
  !! call simps(f, a, b)     ! Using 100 equispaced points between a and b
  !! call simps(f, a, b, N)  ! Using N equispaced points between 0 an 1
  !! call simps(f, x)        ! Using evaluations on values of array x
  !! ```
  !!
  !! Here:
  !! @param [in] f is a real or complex function
  !! @param [in] a (real) Lower limit of integration
  !! @param [in] b (real) Upper limit of integration
  !! @param [in] N (integer) Number of equally-spaced points to use
  !! @param [in] x (real, array) Points where the function will be evaluated
  !!
  !! For sampled values it accepts several input formats,
  !!
  !! ```{.f90}
  !! call simps(y)           ! Using values in array y, corresponding to x ∈ [0,1]
  !! call simps(y, a, b)     ! Using values in array y, corresponding to x ∈ [a,b]
  !! call simps(y, dh)       ! Using values in y, for x values equispaced in dh
  !! call simps(y, x)        ! Using values in y, corresponding to x
  !! ```
  !!
  !! Here:
  !! @param [in] y is a real or complex array
  !! @param [in] a (real) Lower limit of integration
  !! @param [in] b (real) Upper limit of integration
  !! @param [in] dh (real) Space between points in the x-axis
  !! @param [in] x (real, array) Points where the function is evaluated
  !!
  !! Example:
  !! --------
  !! @snippet ex_integr_trapz.f90 Simple
  !!
  interface simps
    ! SAMPLES:
    module procedure d_simpson_dh_sample, c_simpson_dh_sample
    module procedure d_simpson_sample, c_simpson_sample
    module procedure d_simpson_nonlin_sample, c_simpson_nonlin_sample
    ! FUNCTIONS
    module procedure d_simps_func, c_simps_func
    module procedure d_simps_nonlin_func, c_simps_nonlin_func
    module procedure d_simps_func_arg, c_simps_func_arg
    module procedure d_simps_nonlin_func_arg, c_simps_nonlin_func_arg
    ! interface simps2d
    module procedure d_simps2d_func, c_simps2d_func
    module procedure d_simps2d_func_eps, c_simps2d_func_eps
    module procedure d_simps2d_func_arg, c_simps2d_func_arg
    module procedure d_simps2d_func_eps_arg, c_simps2d_func_eps_arg
  end interface simps

  !
  private
  public :: trapz, simps

contains
  !
  !----------------------------------------------------------------------------
  ! Eval 1D and 2D integrals using trapz (2nd order) and simps (4th order) rule.
  !----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------

  ! Define routines for integrations of functions f(x) -------------------------
#define NUMFOR_EVAL_F(f,x) f(x)
  ! Define real(8) functions and integrals  ----------------
#define NUMFOR_KINDR real(dp)
#define PRNM(a) function JOIN(d_,a)

#define NUMFOR_KINDF procedure(funqr)
#include "qtrs1d.inc"
#undef NUMFOR_KINDF

#define NUMFOR_KINDF procedure(funqrnd)
#include "qtrs2d.inc"
#undef NUMFOR_KINDF

#undef PRNM
#undef NUMFOR_KINDR

  ! Define complex(dp) functions and integrals  -------------
#define NUMFOR_KINDR complex(dp)
#define PRNM(a) function JOIN(c_,a)

#define NUMFOR_KINDF procedure(funqc)
#include "qtrs1d.inc"
#undef NUMFOR_KINDF
#define NUMFOR_KINDF procedure(funqcnd)
#include "qtrs2d.inc"
#undef NUMFOR_KINDF

#undef PRNM
#undef NUMFOR_KINDR

#undef NUMFOR_EVAL_F

  ! ----------------------------------------------------------------------------
  ! Define routines for integrations of functions f(x, args) -------------------
#define NUMFOR_EVAL_F(F,x) F%f(x, F%args)

  ! Define real(8) functions and integrals  ----------------
#define NUMFOR_KINDR real(dp)
#define PRNM(a) function JOIN(JOIN(d_,a),_arg)

#define NUMFOR_KINDF type(nf_rfunction)
#include "qtrs1d.inc"
#undef NUMFOR_KINDF
#define NUMFOR_KINDF type(nf_rndfunction)
#include "qtrs2d.inc"
#undef NUMFOR_KINDF

#undef PRNM
#undef NUMFOR_KINDR

  ! Define complex(dp) functions and integrals  -------------

#define NUMFOR_KINDR complex(dp)
#define PRNM(a) function JOIN(JOIN(c_,a),_arg)

#define NUMFOR_KINDF type(nf_cfunction)
#include "qtrs1d.inc"
#undef NUMFOR_KINDF
#define NUMFOR_KINDF type(nf_cndfunction)
#include "qtrs2d.inc"
#undef NUMFOR_KINDF

#undef PRNM
#undef NUMFOR_KINDR

#undef NUMFOR_EVAL_F

  !> get_simpson_weights sets several quadrature weights for integration of order 2,4
  subroutine get_simpson_weights(wt)
    implicit none
    real(dp), dimension(:), intent(OUT) :: wt !< weights
    integer :: n
    n = size(wt)
    select case (n)           ! Should be N >= 3
    case (1)
      wt = 1._dp
    case (2)
      wt = 0.5_dp
    case (4)                  ! Simpson's 3/8 rule
      wt(1) = 3._dp / 8._dp
      wt(2:3) = 9._dp / 8._dp
      wt(4) = 3._dp / 8._dp
    case default            !Simpson's rule n>=5
      if (mod(n - 1, 2) == 0) then
        wt(1) = 1._dp / 3._dp
        wt(n) = 1._dp / 3._dp
        wt(2:n - 1:2) = 4._dp / 3._dp
        wt(3:n - 2:2) = 2._dp / 3._dp
      else
        wt(1) = 1._dp / 3._dp
        wt(2:n - 4:2) = 4._dp / 3._dp
        wt(3:n - 5:2) = 2._dp / 3._dp
        wt(n - 3) = 17._dp / 24._dp
        wt(n - 2) = 9._dp / 8._dp
        wt(n - 1) = 9._dp / 8._dp
        wt(n) = 3._dp / 8._dp
      end if
    end select
  end subroutine get_simpson_weights

end module qsimpson
