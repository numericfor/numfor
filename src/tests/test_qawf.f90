program test_qawf
  use utils, only: dp, M_PI, M_PI_2, timer, str, center
  use integrate, only: qawf, nf_rfunction, nf_cfunction
  implicit none

  real(dp) :: Int1, result
  real(dp) :: a, b
  real(dp) :: omega
  real(dp) :: epsabs, epsrel, abserr

  integer :: Neval, ier
  type(timer) T1
  character(len=:), allocatable :: title

  integer :: i
  integer :: flgw
  ! type(nf_rfunction) :: fun1
  real(dp), dimension(2) :: w
  type(nf_cfunction) :: fun2
  complex(dp) :: Int2

  call T1%start()

  epsabs = 1.e-5_8
  epsrel = 1.e-3_8

  a = 0._dp; b = 1._dp
  omega = M_PI / 2
  flgw = 1

  result = 1.0_8
  title = "| integrate(cos(pi x /2 ) / sqrt(x), x, 0, inf ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  print "(A)", center(title, 60)
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')

  call qawf(fquad457, a, omega, flgw, Int1, epsabs, abserr, Neval, ier)
  print "(A)", ' qawf  = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  fun2%f => fquad457wc
  title = "| integrate((k,k-1) cos(pi x /2 ) / sqrt(x), x, 0, inf )  "
  print "(A)", center(" (Complex, using a variable argument w:) ", 70, '=')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = [real(i - 1, kind=dp), real(i, kind=dp)]
    fun2%args = [w]
    Int2 = cmplx(0._8, 0._8, kind=dp)
    call qawf(fun2, a, omega, flgw, Int2, epsabs, abserr, Neval, ier)
    print "(A)", ' qawf(w='//str(w)//') = '//str(Int2)//" (Error="//str(abserr)//", N="//str(Neval)//")"
  end do

  call T1%stop()
  call T1%show()

contains

  !> Example for QAWF \f$\sqrt(x) \cos(\pi x/2) \f$. Result(w=10 pi) = 1
  function fquad457(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = 0._8
    IF (x > 0) y = 1.0_8 / sqrt(x)
  end function fquad457
  function fquad457wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = 0._8
    IF (x > 0) y = cmplx(w(1), w(2), kind=dp) / sqrt(x)
  end function fquad457wc

end program test_qawf

