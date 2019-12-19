program test_qagi
  use utils, only: dp, Zero, M_PI, M_PI_2, timer, str, center, nf_minf, nf_inf
  use integrate, only: qagi, nf_rfunction, nf_cfunction
  implicit none

  real(dp) :: Int1, result
  real(dp) :: a
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

  epsabs = 1.e-5_8
  epsrel = 1.e-3_8

  a = Zero; 
  result = -0.3616892_8
  title = "| integrate(log(x)/(1 + w x^2), x, 0, \inf ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  print "(A)", center(title, 60)
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  call qagi(fquad455, Zero, nf_inf, Int1, epsabs, epsrel, 'qk21', abserr, Neval, ier)
  print "(A)", ' qagi  = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  call qagi(fquad455, Zero, 20._dp, Int1, epsabs, epsrel, 'qk21', abserr, Neval, ier)
  print "(A)", ' qagi  = '//str(Int1)//" (Error="//str(abserr)//") with N = "//str(Neval)

  fun2%f => fquad455wc
  title = "| integrate(log(abs(x))/(1 + w x^2), x, -\infty, \inf )  "
  print "(A)", center(" (Complex, using a variable argument w:) ", 70, '=')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = [real(i - 1, kind=dp), real(i, kind=dp)]
    fun2%args = [w]
    ! call qagi(fquad455wc, nf_minf, nf_inf, [w], Int2, epsabs, epsrel, 'qk21', abserr, Neval, ier)
    call qagi(fun2, nf_minf, nf_inf, Int2, epsabs, epsrel, 'qk21', abserr, Neval, ier)
    print "(A)", ' qags(w='//str(w)//') = '//str(Int2)//" (Error="//str(abserr)//", N="//str(Neval)//")"
  end do

  call T1%stop()
  call T1%show()

contains
  !> Example for QAGI \f$\log(x)/(1 + w x^2) \f$. Result(w=100) = -0.3616892
  function fquad455(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = log(x) / (1 + 100 * x**2)
  end function fquad455
  function fquad455w(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = log(x) / (1 + w(1) * x**2)
  end function fquad455w

  function fquad455wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = log(abs(x)) / (1 + cmplx(w(1), w(2), kind=dp) * x**2)
  end function fquad455wc

end program test_qagi

