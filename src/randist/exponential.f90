!> @file uniform.f90 for random number generator (and distributions)
!> @author Juan Fiol <juanfiol@gmail.com> (modifications, see real authors below)
!! @date "2020-01-09 14:51:35"

#ifdef __GFORTRAN__
#define PASTE(a) a
#define JOIN(a,b) PASTE(a)b
#else
#define PASTE(a,b) a ## b
#define JOIN(a,b) PASTE(a,b)
#endif

!> @ingroup randomdist
!! Exponential random distribution
module exponential
  USE utils, only: dp, i8, M_DPI
  USE random, only: random_real

  !> Fills a scalar or array with random numbers following an exponential distribution
  !!
  !! The exponential probability distribution located at \f$x_{0}\f$ and mean \f$\mu\f$, is given by
  !! \f[
  !!    p(x) dx = {1 \over \mu} \exp(-x/\mu) dx \quad, \qquad x \ge 0
  !! \f]
  !! The general signature is `call random_exponential(mu, x)` where `x` may be a scalar, vector (1D array), matrix (2D array), 3D-array, 4D-array, or 5D-array.
  !!
  !! The following uses are valid.
  !! ```{.f90}
  !! USE numfor, only: random_exponential
  !! real(dp) :: x
  !! real(dp), dimension(2,3,3) :: xx
  !! call random_exponential(mu, x)           ! Valid.
  !! call random_exponential(mu, xx)          ! Valid.
  !! !
  !! ```
  !! Example:
  !! --------
  !! @include ex_rand_expon.f90
  interface random_exponential
    module procedure :: expon0d, expon1d, expon2d
    module procedure :: expon3d, expon4d, expon5d
  end interface random_exponential

  !> This function computes the probability density p(x) at x for an exponential distribution
  ! with mean `mu`, using the formula given above.
  interface ran_exponential_pdf
    module procedure :: ran_expon_pdf_s, ran_expon_pdf_v
  end interface ran_exponential_pdf

  private
  public :: random_exponential, ran_exponential_pdf

contains

  function rng_expon(mu) result(r)
    implicit none
    real(dp) :: r !< Random number
    real(dp), intent(IN) :: mu
    real(dp) :: log1p
    r = random_real()

    if (r < 1.e-5_dp) then
      log1p = -r * (1._dp + r * (0.5_dp + r / 3))
    else
      log1p = log(1._dp - r)
    endif
    r = -mu * log1p
  end function rng_expon

  function ran_expon_pdf_s(mu, x) result(y)
    implicit none
    real(dp), intent(IN) :: x !<
    real(dp), intent(IN) :: mu !<
    real(dp) :: y !<
    y = 0
    IF (x >= 0) y = exp(-x / mu) / mu
  end function ran_expon_pdf_s

  function ran_expon_pdf_v(mu, x) result(y)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !<
    real(dp), intent(IN) :: mu !<
    real(dp), dimension(size(x)) :: y
    y = 0
    where (x >= 0) y = exp(-x / mu) / mu
  end function ran_expon_pdf_v

! Name of arguments
#define A2 mu
#define RNG_NAME rng_expon
#define ROUTINE_NAME expon
#include "rand_tmpl_2.inc"
#undef ROUTINE_NAME
#undef RNG_NAME
#undef A2

end module exponential

! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:
