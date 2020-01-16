!> @file uniform.f90 for random number generator (and distributions)
!> @author Juan Fiol <juanfiol@gmail.com> (modifications, see real authors below)
!! @date "2020-01-09 12:39:36"

#ifdef __GFORTRAN__
#define PASTE(a) a
#define JOIN(a,b) PASTE(a)b
#else
#define PASTE(a,b) a ## b
#define JOIN(a,b) PASTE(a,b)
#endif

!> @ingroup randomdist
!! Uniform Random distribution
module uniform
  USE utils, only: dp, i8
  USE random, only: random_real, random_int

  ! !> @ingroup randomdist
  !> Fills a scalar or array with random numbers in the half-open interval [0, 1).
  !!
  !! The signature is `call random_sample(x)` where `x` may be a scalar, vector (1D array), matrix (2D array), 3D-array, 4D-array, or 5D-array.
  !!
  !! To sample uniformly over \f$[a, b)\f$ with \f$b > a\f$,  multiply the output of random_sample by \f$(b-a)\f$ and add \f$a\f$:
  !! ```{.f90}
  !! call random_sample(x)
  !! y =(b - a) * x + a
  !! ```
  !!
  !! @note  The subroutine is implemented to work with arrays of up to five dimensions.
  !! Versions for arrays of larger dimensionality may be implemented straightforwardly.
  !!
  !! Example:
  !! --------
  !! @snippet ex_rand_unif.f90 usesample
  interface random_sample
    module procedure :: random_sample0d, random_sample1d, random_sample2d
    module procedure :: random_sample3d, random_sample4d, random_sample5d
  end interface random_sample

  ! !> @ingroup randomdist
  !> Fills a scalar or array with random numbers in the half-open interval [low, high).
  !!
  !! The signature is `call random_uniform(low, high, x)` where `x` may be a scalar, vector (1D array), matrix (2D array), 3D-array, 4D-array, or 5D-array.
  !!
  !! Example:
  !! --------
  !! @snippet ex_rand_unif.f90 useuniform
  interface random_uniform
    module procedure :: random_unif0d, random_unif1d, random_unif2d
    module procedure :: random_unif3d, random_unif4d, random_unif5d
  end interface random_uniform

  private
  public :: random_sample, random_uniform

contains

#define RNG_NAME random_real

#define ROUTINE_NAME random_sample
#include "rand_tmpl_1.inc"
#undef ROUTINE_NAME

#define ROUTINE_NAME random_unif
#include "rand_tmpl_3.inc"
#undef ROUTINE_NAME

#undef RNG_NAME

end module uniform
! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:
