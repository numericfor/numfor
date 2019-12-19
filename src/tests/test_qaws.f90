program test_qaws
  use utils, only: dp, M_PI, M_PI_2, timer, str, center
  use integrate, only: qaws, nf_rfunction, nf_cfunction
  implicit none

  real(dp) :: Int1, result
  real(dp) :: a, b
  real(dp) :: alfa, beta
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

  epsabs = 0.e-5_8
  epsrel = 1.e-5_8

  a = 0._dp; b = 1._dp
  alfa = 0._8; beta = 0._8
  flgw = 2

  result = -0.1892752_8

  title = "| integrate(ln(x)/(1 + ln(x)^{2})^{2}, x, 0, 1 ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')

  call qaws(fquad458, a, b, alfa, beta, flgw, Int1, epsabs, epsrel, abserr, Neval, ier)
  print "(A)", ' qaws  = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  fun2%f => fquad458wc
  title = "| integrate((k,k-1) ln(x)/(1 + ln(x)^{2})^{2}, x, 0, 1 )  "
  print "(A)", center(" (Complex, using a variable argument w:) ", 70, '=')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = [real(i - 1, kind=dp), real(i, kind=dp)]
    fun2%args = [w]
    Int2 = cmplx(0._8, 0._8, kind=dp)
    call qaws(fun2, a, b, alfa, beta, flgw, Int2, epsabs, epsrel, abserr, Neval)
    print "(A)", ' qaws(w='//str(w)//') = '//str(Int2)//" (Error="//str(abserr)//", N="//str(Neval)//")"
  end do

  call T1%stop()
  call T1%show()

contains

  !> Example for QAWS \f$\ln(x)/(1 + ln(x)^{2})^{2} \f$. Result = 0.1892752
  function fquad458(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = 0._8
    IF (x > 0._8) y = 1._8 / (1._8 + log(x)**2)**2
  end function fquad458
  function fquad458wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), dimension(:), intent(IN) :: w !<
    real(dp), intent(IN) :: x !<
    y = 0._8
    IF (x > 0) y = cmplx(w(1), w(2), kind=dp) / (1 + log(x)**2)**2
  end function fquad458wc

end program test_qaws

