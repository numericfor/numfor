program test_qthsh
  use utils, only: dp, M_PI, M_PI_2, timer, str, center
  use integrate, only: qag, qags, qng, qnthsh, nf_rfunction, nf_cfunction
  implicit none

  complex(dp), dimension(5) :: results
  real(dp) :: Int1, result
  real(dp) :: a, b
  real(dp) :: epsabs, epsrel, abserr

  integer :: Neval

  type(timer) T1
  character(len=:), allocatable :: title
  ! integer :: i
  integer, parameter :: width = 70

  integer :: i
  ! type(nf_rfunction) :: fun1
  real(dp) :: w
  real(dp) :: err
  intrinsic dsin

  ! type(nf_cfunction) :: fun2
  ! complex(dp) :: Int2

  ! Results for integrate(exp(i w sin(x)), x, 0, pi ) from maxima
  results = [cmplx((0.3024310191244698, -0.3198921069678047), kind=dp), &
    & cmplx((0.1753395985854673, -0.2680962123081692), kind=dp),&
    & cmplx((0.1088370651016893, -0.2415328394565208), kind=dp), &
    & cmplx((0.06278740049149098, -0.2226721656038109), kind=dp), &
    & cmplx((0.02699336268296957, -0.2065688529554112), kind=dp)]

  call T1%start()

  epsabs = 1.e-7_8
  epsrel = 1.e-3_8

  a = 0._dp; b = 1._dp
  result = -4._dp / 9
  title = "| integrate( sqrt(x) log(x), x, 0, 1 ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), width, ' ')
  print "(A)", center(title, width)
  print "(A)", center(repeat('-', len(title) - 2), width, ' ')
  call qag(fquad451, a, b, Int1, epsabs, epsrel, 1, abserr, Neval)
  print "(A)", ' qag          = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  call qnthsh(fquad451, a, b, Int1, epsabs, epsrel, abserr, Neval)
  print "(A)", ' th-sh        = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  a = 0; b = M_PI
  result = 2.403939430634413_8
  title = "| integrate( cos(sin(x)), x, 0, pi ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), width, ' ')
  print "(A)", center(title, width)
  print "(A)", center(repeat('-', len(title) - 2), width, ' ')
  call qnthsh(fquad452, a, b, Int1, epsabs, epsrel, abserr, Neval)
  print "(A)", ' th-sh = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  ! call qag(fquad452, a, b, Int1, epsabs, epsrel, 4, abserr, Neval)
  ! print "(A)", ' qag   = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  ! call qng(fquad452, a, b, Int1, epsabs, epsrel, abserr, Neval)
  ! print "(A)", ' qng   = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  result = real(results(1), kind=dp)
  title = "| integrate( cos(w * sin(x)), x, 0, pi ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), width, ' ')
  print "(A)", center(title, width)
  print "(A)", center(repeat('-', len(title) - 2), width, ' ')
  a = 0._dp; b = M_PI
  ! ! fun1%f => fquad452w
  ! ! do i = 1, 5
  ! !   w = 25 * i
  ! !   result = real(results(i), kind=dp)
  ! !   fun1%args = [w]
  ! !   call qag(fun1, a, b, Int1, epsabs, epsrel, 5, abserr, Neval)
  ! !   err = abs(Int1 - result)
  ! !   print "(A)", ' maxima qags = '//str(result)//" (Diff="//str(err)//") with w = "//str(w)
  ! !   print "(A)", ' qag (k=51)  = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  ! !   call qnthsh(fun1, a, b, Int1, epsabs, epsrel, abserr, Neval)
  ! !   err = abs(Int1 - result)
  ! !   ! print "(A)", ' maxima qags = '//str(result)//" (Diff="//str(err)//") with N = "//str(Neval)
  ! !   print "(A)", '      qnthsh = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  ! ! end do

  print "(A)", center(" Calling the function f(x, args) directly ", width, '*')
  do i = 1, 5
    w = 25 * i
    result = real(results(i), kind=dp)
    call qag(fquad452w, a, b, [w], Int1, epsabs, epsrel, 5, abserr, Neval)
    err = abs(Int1 - result)
    print "(A)", ' maxima qags = '//str(result)//" (Diff="//str(err)//") with w = "//str(w)
    print "(A)", ' qag (k=51)  = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
    call qnthsh(fquad452w, a, b, [w], Int1, epsabs, epsrel, abserr, Neval)
    err = abs(Int1 - result)
    print "(A)", '      qnthsh = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)
  end do
  call qnthsh(dsin, a, b, Int1, epsabs, epsrel, abserr, Neval)
  print "(A)", '\int sin(x) dx = '//str(Int1)//" (Error="//str(abs(cos(a) - cos(b) - Int1))//") with N = "//str(Neval)

  call T1%stop()
  call T1%show()

contains
  !> \f$\cos(1 \sin(x)) \f$
  function fquad452(x) result(y)
    implicit none
    real(8) :: y !<
    real(8), intent(IN) :: x !<
    y = cos(sin(x))
  end function fquad452
  !> \f$\sqrt(x) \log(x) \f$
  function fquad451(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = sqrt(x) * log(x)
  end function fquad451

  !> \f$\cos(w \sin(x)) \f$
  function fquad452w(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = cos(w(1) * sin(x))
  end function fquad452w

end program test_qthsh

