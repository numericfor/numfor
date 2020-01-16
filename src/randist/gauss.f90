!> @file uniform.f90 for random number generator (and distributions)
!> @author Juan Fiol <juanfiol@gmail.com> (modifications, see real authors below)
!! @date "2020-01-09 13:00:01"

#ifdef __GFORTRAN__
#define PASTE(a) a
#define JOIN(a,b) PASTE(a)b
#else
#define PASTE(a,b) a ## b
#define JOIN(a,b) PASTE(a,b)
#endif

!> @ingroup randomdist
!! Normal random distribution
module gauss
  USE utils, only: dp, i8, M_DPI
  USE random, only: random_real, random_real_pos

  !> Functions returning a variate x with normal probability
  !!
  !! Valid uses:
  !! ```{.f90}
  !! USE numfor, only rng_normal
  !! x = rng_normal()
  !! x = rng_normal(scale)
  !! x = rng_normal(loc, scale)
  !! ```
  interface rng_normal
    module procedure :: ran_ugaussian, ran_gaussian_sym, ran_gaussian
  end interface rng_normal

  !> Functions returning a variate x with normal probability
  !!
  !! Valid uses:
  !! ```{.f90}
  !! USE numfor, only rng_normal
  !! x = rng_normal2()
  !! x = rng_normal2(scale)
  !! x = rng_normal2(loc, scale)
  !! ```
  !! This version is faster than rng_normal() but is not thread-safe
  interface rng_normal2
    module procedure :: ran_ugaussian2, ran_gaussian_sym2, ran_gaussian2
  end interface rng_normal2

  !> @brief Convenience routine. Fills a scalar or array with random numbers following a standard normal (gaussian) distribution. Equivalently to random_normal() with `scale=1` and `loc=0`.
  interface random_standard_normal
    module procedure :: std_norm0d, std_norm1d, std_norm2d
    module procedure :: std_norm3d, std_norm4d, std_norm5d
  end interface random_standard_normal

  ! !> @ingroup randomdist

  !> Fills a scalar or array with random numbers following a normal (gaussian) distribution
  !!
  !! The normal probability distribution located at \f$x_{0}\f$ and standard deviation \f$\sigma\f$, is given by
  !! \f[ p(x) dx = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-(x-x_{0})^2 / 2\sigma^2) dx \f]
  !! The general signature is `call random_normal([[loc,] scale,] x)` where the center `loc`=\f$x_{0}\f$ and the `scale`=\f$\sigma\f$ are optional, and `x` may be a scalar, vector (1D array), matrix (2D array), 3D-array, 4D-array, or 5D-array.
  !!
  !! The following uses are all valid, except the last one
  !! ```{.f90}
  !! USE numfor, only: random_normal
  !! call random_normal(x)               ! Valid. Standard normal (loc=0, scale=1)
  !! call random_normal(scale, x)        ! Valid. loc=0
  !! call random_normal(loc, scale, x)   ! Valid
  !! call random_normal(loc, x)          ! INVALID
  !! !
  !! ```
  !! Example:
  !! --------
  !! @include ex_rand_gauss.f90
  interface random_normal
    module procedure :: sym_norm0d, sym_norm1d, sym_norm2d
    module procedure :: sym_norm3d, sym_norm4d, sym_norm5d
    module procedure :: norm0d, norm1d, norm2d
    module procedure :: norm3d, norm4d, norm5d
    module procedure :: std_norm0d, std_norm1d, std_norm2d
    module procedure :: std_norm3d, std_norm4d, std_norm5d
  end interface random_normal

  !> This function computes the probability density p(x) at x for a Gaussian distribution
  ! with standard deviation `scale`, using the formula given above.
  interface ran_gaussian_pdf
    module procedure :: ran_gauss_pdf_s, ran_gauss_pdf_v
  end interface ran_gaussian_pdf

  private
  public :: random_normal, random_standard_normal, ran_ugaussian2
  public :: ran_gaussian_pdf

contains

  !> @ingroup randomdist
  !> Computes a Gaussian random variate, with mean zero and standard deviation `sigma = 1`
  function ran_ugaussian() result(r)
    implicit none
    real(dp) :: r !< Random number
    real(dp) :: x, y
    real(dp) :: r2

    ! Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122
    ! Copied from gsl (gauss.c).
    ! Differences with gsl: seed is not concerned. Independent setting is supposed
    r2 = 3
    do while (r2 > 1 .OR. r2 == 0._dp) ! See if it is in the unit circle
      x = -1 + 2 * random_real_pos()
      y = -1 + 2 * random_real_pos()
      r2 = x * x + y * y
    end do
    r = y * sqrt(-2 * log(r2) / r2)
  end function ran_ugaussian

  !> computes a Gaussian random variate, with mean zero and standard deviation `sigma`
  function ran_gaussian_sym(scale) result(r)
    implicit none
    real(dp) :: r !< Random number
    real(dp), intent(IN) :: scale !< standard deviation

    r = scale * ran_ugaussian()
  end function ran_gaussian_sym

  !> ran_gaussian computes a Gaussian random variate, with mean zero and standard deviation `sigma`
  function ran_gaussian(loc, scale) result(r)
    implicit none
    real(dp) :: r !< Random number
    real(dp), intent(IN) :: loc !< Center of distribution
    real(dp), intent(IN) :: scale !< standard deviation

    r = scale * ran_ugaussian() + loc
  end function ran_gaussian

  !> Computes a Gaussian random variate, with mean zero and standard deviation `sigma=1`.
  !! @note This version caches the last unused value.
  !! It should be faster but will have problems in multi-threaded programs
  function ran_ugaussian2() result(r)
    implicit none
    real(dp) :: r !< Random number
    real(dp) :: x, y
    real(dp) :: r2
    real(dp), save :: rr
    logical, save :: has_gauss = .False.

    ! Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122
    ! Copied from gsl (gauss.c).
    ! Differences with gsl: seed is not concerned. Independent setting is supposed
    if (has_gauss) then
      r = rr
      has_gauss = .False.
    else
      r2 = 3
      do while (r2 > 1._dp .OR. r2 == 0._dp) ! See if it is in the unit circle
        x = -1 + 2 * random_real_pos()
        y = -1 + 2 * random_real_pos()
        r2 = x * x + y * y
      end do
      r2 = sqrt(-2 * log(r2) / r2)
      r = x * r2
      rr = y * r2
      has_gauss = .True.
    end if
  end function ran_ugaussian2

  function ran_gaussian_sym2(scale) result(r)
    implicit none
    real(dp) :: r !< Random number
    real(dp), intent(IN) :: scale !< standard deviation

    r = scale * ran_ugaussian2()
  end function ran_gaussian_sym2

  !> ran_gaussian computes a Gaussian random variate, with mean zero and standard deviation `sigma`
  function ran_gaussian2(loc, scale) result(r)
    implicit none
    real(dp) :: r !< Random number
    real(dp), intent(IN) :: loc !< Center of distribution
    real(dp), intent(IN) :: scale !< standard deviation

    r = scale * ran_ugaussian2() + loc
  end function ran_gaussian2

  !> ran_gauss_pdf Computes
  !!
  function ran_gauss_pdf_s(x, sigma) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), intent(IN) :: sigma !<
    real(dp) :: u
    u = x / abs(sigma)
    y = (1 / (sqrt(M_DPI) * abs(sigma))) * exp(-u * u * 0.5_dp)
  end function ran_gauss_pdf_s

  !> ran_gauss_pdf Computes
  !!
  function ran_gauss_pdf_v(x, sigma) result(y)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !<
    real(dp), dimension(size(x)) :: y
    real(dp), intent(IN) :: sigma !<
    y = (1 / (sqrt(M_DPI) * abs(sigma))) * exp(-x**2 / (2 * sigma**2))
  end function ran_gauss_pdf_v

  ! JF: This is thread-safe but takes almost twice
! #define RNG_NAME ran_ugaussian2
  ! JF: We are using the non-thread safe version
#define RNG_NAME ran_ugaussian2
#define ROUTINE_NAME std_norm
#include "rand_tmpl_1.inc"
#undef ROUTINE_NAME
#undef RNG_NAME

! Name of arguments
#define A1 loc
#define A2 scale

#define RNG_NAME ran_gaussian_sym2
#define ROUTINE_NAME sym_norm
#include "rand_tmpl_2.inc"
#undef ROUTINE_NAME
#undef RNG_NAME

#define RNG_NAME ran_gaussian2
#define ROUTINE_NAME norm
#include "rand_tmpl_3.inc"
#undef ROUTINE_NAME
#undef RNG_NAME

#undef A1
#undef A2

end module gauss

! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:
