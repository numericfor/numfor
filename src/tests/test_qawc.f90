program test_qawc
  use utils, only: dp, M_PI, M_PI_2, timer, str, center
  use integrate, only: qawc, nf_rfunction, nf_cfunction
  implicit none

  real(dp) :: Int1, result
  real(dp) :: a, b
  real(dp) :: c
  real(dp) :: epsabs, epsrel, abserr

  integer :: Neval, ier
  type(timer) T1
  character(len=:), allocatable :: title

  integer :: i
  ! type(nf_rfunction) :: fun1
  real(dp), dimension(2) :: w
  type(nf_cfunction) :: fun2
  complex(dp) :: Int2

  call T1%start()

  epsabs = 0._8
  epsrel = 1.e-3_8

  a = -1._dp; b = 5._dp
  c = 0._dp

  result = -0.08994401_8

  title = "| integrate(1/(x (5 x^3 + 6)), x, 0, 1 ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')

  call qawc(fquad459, a, b, c, Int1, epsabs, epsrel, abserr, Neval, ier)
  print "(A)", ' qawc  = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  fun2%f => fquad459wc
  title = "| integrate(1/(x ((5 + k-1) x^3 + 6 + j k)), x, 0, 1 )  "
  print "(A)", center(" (Complex, using a variable argument w:) ", 70, '=')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = [real(i - 1, kind=dp), real(i, kind=dp)]
    fun2%args = [w]
    Int2 = cmplx(0._8, 0._8, kind=dp)
    call qawc(fun2, a, b, c, Int2, epsabs, epsrel, abserr, Neval, ier)
    print "(A)", ' qawc(w='//str(w)//') = '//str(Int2)//" (Error="//str(abserr)//", N="//str(Neval)//")"
  end do

  call T1%stop()
  call T1%show()

contains

  !> Example for QAWC \f$ \frac{1}{x (w(1) x^3+ (w(2)))} \f$. Result(w=[5,6]) = -0.08994401
  function fquad459(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = 1 / (5 * x**3 + 6)
  end function fquad459
  function fquad459wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = 1 / ((5 + w(1)) * x**3 + 6 + cmplx(0._8, 1._8, kind=dp) * w(2))
  end function fquad459wc

end program test_qawc

