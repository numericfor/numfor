program test_qawo
  use utils, only: dp, M_PI, M_PI_2, timer, str, center
  use integrate, only: qawo, nf_rfunction, nf_cfunction
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
  omega = 10 * M_PI
  flgw = 2

  result = -0.12181316_8
  title = "| integrate(log(x) sin(10 pi x), x, 0, 1 ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  print "(A)", center(title, 60)
  call qawo(fquad456, a, b, omega, flgw, Int1, epsabs, epsrel, abserr, Neval, ier)
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')

  print "(A)", ' qawo  = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  fun2%f => fquad456wc
  title = "| integrate((k,k-1) log(x) sin(10 pi x), x, 0,1 )  "
  print "(A)", center(" (Complex, using a variable argument w:) ", 70, '=')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = [real(i - 1, kind=dp), real(i, kind=dp)]
    fun2%args = [w]
    call qawo(fun2, a, b, omega, flgw, Int2, epsabs, epsrel, abserr, Neval, ier)
    print "(A)", ' qawo(w='//str(w)//') = '//str(Int2)//" (Error="//str(abserr)//", N="//str(Neval)//")"
  end do

  call T1%stop()
  call T1%show()

contains
  !> Example for QAWO \f$\log(x) \sin(10 \pi x) \f$. Result(w=10 pi) = -0.12181316
  function fquad456(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = 0._8
    IF (x > 0) y = log(x)
  end function fquad456
  function fquad456wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = cmplx(0._8, 0._8, kind=8)
    IF (x > 0) y = cmplx(1 + w(1), w(1), kind=dp) * log(w(2) * x)
  end function fquad456wc

end program test_qawo

