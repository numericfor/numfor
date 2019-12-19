program test_qag
  use utils, only: dp, M_PI, M_PI_2, timer, str, center
  use integrate, only: qag, nf_rfunction, nf_cfunction
  ! use integrate, only: trapz, simps, qags!, d_trapz_ab_func
  implicit none

  real(dp) :: Int1, result
  complex(dp), dimension(5) :: results
  real(dp) :: a, b
  real(dp) :: epsabs, epsrel, abserr

  integer :: Neval, ier
  type(timer) T1
  character(len=:), allocatable :: title
  integer :: i

  type(nf_rfunction) :: fun1
  real(dp) :: w
  real(dp) :: err
  type(nf_cfunction) :: fun2
  complex(dp) :: Int2

  ! Results for integrate(exp(i w sin(x)), x, 0, pi ) from maxima
  results = [cmplx((0.3024310191244698, -0.3198921069678047), kind=dp), &
    & cmplx((0.1753395985854673, -0.2680962123081692), kind=dp),&
    & cmplx((0.1088370651016893, -0.2415328394565208), kind=dp), &
    & cmplx((0.06278740049149098, -0.2226721656038109), kind=dp), &
    & cmplx((0.02699336268296957, -0.2065688529554112), kind=dp)]

  call T1%start()

  epsabs = 1.e-8_8
  epsrel = 1.e-3_8

  a = 0._dp; b = 1._dp
  result = -4._dp / 9
  title = "| integrate( sqrt(x) log(x), x, 0, 1 ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  print "(A)", center(title, 60)
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  call qag(fquad451, a, b, Int1, epsabs, epsrel, "qk15", abserr, Neval, ier)
  print "(A)", ' qag          = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  result = 0.06278740_8
  title = "| integrate( cos(100 * sin(x)), x, 0, pi ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  print "(A)", center(title, 60)
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  a = 0._dp; b = M_PI
  call qag(fquad452, a, b, Int1, epsabs, epsrel, "qk15", abserr, Neval, ier)
  print "(A)", ' qag          = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  title = "| integrate( cos(w * sin(x)), x, 0, pi ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  print "(A)", center(title, 60)
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  a = 0._dp; b = M_PI
  fun1%f => fquad452w
  do i = 1, 5
    w = 25 * i
    fun1%args = [w]
    call qag(fun1, a, b, Int1, epsabs, epsrel, 'qk21', abserr, Neval, ier)
    err = abs(Int1 - real(results(i), kind=dp))
    print "(A)", ' maxima qags = '//str(Int1)//" (Diff="//str(err)//") with N = "//str(Neval)
    print "(A)", ' qag (k=21)  = '//str(Int1)//" (Error="//str(abserr)//") with w = "//str(w)
    call qag(fun1, a, b, Int1, epsabs, epsrel, 'qk31', abserr, Neval, ier)
    err = abs(Int1 - real(results(i), kind=dp))
    print "(A)", ' maxima qags = '//str(Int1)//" (Diff="//str(err)//") with N = "//str(Neval)
    print "(A)", ' qag (k=31)   = '//str(Int1)//" (Error="//str(abserr)//") with w = "//str(w)
  end do

  title = "| integrate( sin(w * sin(x)), x, 0, pi ) = "//str(result)//" |"
  print "(A)", center(" (Using a variable argument:) ", 70, '-')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  fun1%f => fquad452ws
  do i = 1, 5
    w = 25._dp * i
    result = aimag(results(i))
    fun1%args = [w]
    call qag(fun1, a, b, Int1, epsabs, epsrel, 'qk21', abserr, Neval, ier)
    err = abs((result - Int1) / result)
    print "(A)", 'maxima qags  = '//str(result)//" - Error="//str(err)
    print "(A)", ' qag (k=21, w='//str(w)//') = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
    call qag(fun1, a, b, Int1, epsabs, epsrel, 'qk31', abserr, Neval, ier)
    err = abs((result - Int1) / result)
    print "(A)", 'maxima qags  = '//str(result)//" - Error="//str(err)
    print "(A)", ' qag (k=31, w='//str(w)//') = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  end do

  fun2%f => fquad452wc
  title = "| integrate( exp( i(w * sin(x)) ), x, 0, pi ) |"
  print "(A)", center(" (Complex, using a variable argument w:) ", 70, '=')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = 25._dp * i
    fun2%args = [w]
    call qag(fun2, a, b, Int2, epsabs, epsrel, "qk41", abserr, Neval, ier)
    err = abs((results(i) - Int2) / results(i))
    print "(A)", '  qag(w='//str(w)//') = '//str(Int2)//" (Error="//str(err)//", Estimated="//str(abserr)//")"
  end do

  fun1%f => fquad452w
  title = "| integrate( cos(w * sin(x)), x, 0, pi ) |"
  print "(A)", center(" (Integrating the Real part) ", 70, '=')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = 25._dp * i
    result = real(results(i))
    fun1%args = [w]
    call qag(fun1, a, b, Int1, epsabs, epsrel, "qk41", abserr, Neval, ier)
    err = abs((result - Int1) / result)
    print "(A)", '  qag(w='//str(w)//') = '//str(Int1)//" (Error="//str(err)//", Estimated="//str(abserr)//")"
  end do

  fun1%f => fquad452ws
  title = "| integrate( sin(w * sin(x)), x, 0, pi ) |"
  print "(A)", center(" (Integrating the Imag part) ", 70, '=')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = 25._dp * i
    result = aimag(results(i))
    fun1%args = [w]
    call qag(fun1, a, b, Int1, epsabs, epsrel, "qk41", abserr, Neval, ier)
    err = abs((result - Int1) / result)
    print "(A)", '  qag(w='//str(w)//') = '//str(Int1)//" (Error="//str(err)//", Estimated="//str(abserr)//")"
  end do

  call T1%stop()
  call T1%show()

contains
  !> \f$\sqrt(x) \log(x) \f$
  function fquad451(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = sqrt(x) * log(x)
  end function fquad451

  !> \f$\cos(100 \sin(x)) \f$
  function fquad452(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = cos(100 * sin(x))
  end function fquad452

  !> \f$\cos(w \sin(x)) \f$
  function fquad452wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = cmplx(cos(w(1) * sin(x)), sin(w(1) * sin(x)), kind=dp)
  end function fquad452wc

  !> \f$\cos(w \sin(x)) \f$
  function fquad452w(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = cos(w(1) * sin(x))
  end function fquad452w

  !> \f$\cos(w \sin(x)) \f$
  function fquad452ws(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = sin(w(1) * sin(x))
  end function fquad452ws

end program test_qag

