program test_qagp
  use utils, only: dp, M_PI, M_PI_2, timer, str, center
  use integrate, only: qagp, nf_rfunction, nf_cfunction
  implicit none

  real(dp) :: Int1, result
  complex(dp), dimension(5) :: results
  real(dp) :: a, b
  real(dp) :: epsabs, epsrel, abserr

  integer :: Neval, ier
  ! integer, parameter :: npts = 4
  real(dp), dimension(2) :: points
  type(timer) T1
  character(len=:), allocatable :: title

  integer :: i
  integer :: width = 70

  ! type(nf_rfunction) :: fun1
  real(dp), dimension(2) :: w
  type(nf_cfunction) :: fun2
  complex(dp) :: Int2

  ! Results for integrate(exp(i w sin(x)), x, 0, pi ) from maxima
  results = [cmplx((0.3024310191244698, -0.3198921069678047), kind=dp), &
    & cmplx((0.1753395985854673, -0.2680962123081692), kind=dp),&
    & cmplx((0.1088370651016893, -0.2415328394565208), kind=dp), &
    & cmplx((0.06278740049149098, -0.2226721656038109), kind=dp), &
    & cmplx((0.02699336268296957, -0.2065688529554112), kind=dp)]

  call T1%start()

  epsabs = 1.e-5_8
  epsrel = 1.e-3_8

  a = 0._dp; b = 3._dp
  points(:2) = [1._8, sqrt(2._8)]

  result = 52.740748_8
  title = "| integrate(x^3 \log{|(x^2-1)(x^2-2)|}, x, 0, 3 ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), width, ' ')
  print "(A)", center(title, width)
  print "(A)", center(repeat('-', len(title) - 2), width, ' ')
  call qagp(fquad454, a, b, points, Int1, epsabs, epsrel, abserr, Neval, ier)
  print "(A)", ' qagp  = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  fun2%f => fquad454wc
  title = "| complex version with additional arguments |"
  print "(A)", center(" (Complex, using a variable argument w:) ", 70, '=')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = [real(i - 1, kind=dp), real(i, kind=dp)]
    fun2%args = [w]
    call qagp(fun2, a, b, points, Int2, epsabs, epsrel, abserr, Neval)
    print "(A)", '  qagp(w='//str(w)//') = '//str(Int2)//" (Error="//str(abserr)//", N="//str(Neval)//")"
  end do

  call T1%stop()
  call T1%show()

contains
  !> Example for QAGP \f$\log{|(x^2-1)(x^2-2)|} \f$. Result = 52.740748
  function fquad454(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = x**3 * log(abs((x**2 - 1._dp) * (x**2 - 2._dp)))
  end function fquad454

  function fquad454wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = cmplx(w(1), w(2), kind=dp) * x**3 * log(abs((x**2 - 1._dp) * (x**2 - 2._dp)))
  end function fquad454wc
end program test_qagp

