program test_qag
  use utils, only: dp, M_PI, M_PI_2, timer, str, center
  use integrate, only: quad, nf_rfunction, nf_cfunction, d_qp_extra, c_qp_extra
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

  type(d_qp_extra) :: info
  type(c_qp_extra) :: infoc
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

  a = 0._dp; b = 1._dp
  result = -4._dp / 9
  title = "| integrate( sqrt(x) log(x), x, 0, 1 ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  print "(A)", center(title, 60)
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  call quad(fquad451, a, b, Int1, abserr, epsabs, epsrel, maxsub=100, neval=Neval, ier=ier)
  print "(A)", 'quad semifull = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  call quad(fquad451, a, b, Int1, abserr, neval=Neval)
  print "(A)", 'quad default  = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  call quad(fquad451, a, b, Int1, abserr, neval=Neval, info=info)
  print "(A)", 'quad info(0)  = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  print "(A)", "INFO: last,size = "//str(info%last)//", "//str(size(info%alist))
  print "(A)", "INFO: rlist(:last) = "//str(info%rlist(:info%last))
  info = 0
  info = 50
  ! info = d_qp_extra(50)
  ! call info%init(50)
  call quad(fquad451, a, b, Int1, abserr, neval=Neval, info=info)
  print "(A)", 'quad info(50) = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  print "(A)", "INFO: last,size = "//str(info%last)//", "//str(size(info%alist))
  print "(A)", "INFO: rlist(:last) = "//str(info%rlist(:info%last))

  title = "| integrate( cos(w * sin(x)), x, 0, pi ) |"
  print "(A)", center(" (Integrating the Real part) ", 70, '=')
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  print "(A)", center(title, 60)
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  a = 0._dp; b = M_PI
  fun1%f => fquad452w
  do i = 1, 5
    result = real(results(i), kind=dp)
    w = 25 * i
    fun1%args = [w]
    call quad(fun1, a, b, Int1, abserr, epsrel=1.49e-8_8, neval=Neval)
    err = abs(Int1 - result)
    print "(A)", ' maxima qags = '//str(result)//" (Diff="//str(err)//") with w = "//str(w)
    print "(A)", ' numfor quad = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
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
    call quad(fun2, a, b, Int2, abserr, epsabs, epsrel, neval=Neval, info=infoc)
    err = abs((results(i) - Int2) / results(i))
    print "(A)", '  quad(w='//str(w)//') = '//str(Int2)//" (Error="//str(err)//", N="//str(Neval)//")"
  end do
  print "(A)", 'quad infoc(0) = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  print "(A)", "INFOC: last,size = "//str(infoc%last)//", "//str(size(infoc%alist))
  print "(A)", "INFOC: rlist(:last) = "//str(infoc%rlist(:infoc%last))

  title = "| integrate( exp( i(w * sin(x)) ), x, 0, pi ) |"
  print "(A)", center(" (Complex, using a variable argument w:) ", 70, '=')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = 25._dp * i
    call quad(fquad452wc, a, b, [w], Int2, abserr, epsabs, epsrel, neval=Neval, info=infoc)
    err = abs((results(i) - Int2) / results(i))
    print "(A)", '  quad(w='//str(w)//') = '//str(Int2)//" (Error="//str(err)//", N="//str(Neval)//")"
  end do
  print "(A)", 'quad infoc(0) = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  print "(A)", "INFOC: last,size = "//str(infoc%last)//", "//str(size(infoc%alist))
  print "(A)", "INFOC: rlist(:last) = "//str(infoc%rlist(:infoc%last))

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

  !> \f$\log(x)/\sqrt{(x)} \f$
  function fquad453(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = log(x) / sqrt(x)
  end function fquad453

  !> \f$\log(x)/(1 + ln(x)^{2})^{2} \f$
  function fquad458(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = log10(x) / (1 + log(x)**2)**2
  end function fquad458

  !> \f$\log(x)/(1 + ln(x)^{2})^{2} \f$
  function fquad458w(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = log10(x) / (w(1) + log(x)**2)**2
  end function fquad458w

  !> \f$\cos(w \sin(x)) \f$
  function fquad452ws(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = sin(w(1) * sin(x))
  end function fquad452ws

end program test_qag

