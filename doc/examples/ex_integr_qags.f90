program ex_qags
  USE numfor, only: dp, Zero, M_PI, nf_inf, nf_minf, str, center, qags
  implicit none

  real(dp) :: err
  integer :: N
  real(dp) :: w
  !< [Simple]
  real(dp) :: Integ1

  intrinsic dsin
  print "(A)", center(" Integrate sin(x) between 0 and pi ", 70, '-')
  call qags(dsin, Zero, M_PI, Integ1)
  print "(A)", 'integrate(sin(x), 0, pi) = '//str(Integ1)//" (Difference="//str(abs(2 - Integ1))//")"
  ! that prints:
  !------------------ Integrate sin(x) between 0 and pi -----------------
  ! integrate(sin(x), 0, pi) = 1.9999999827467 (Difference=1.725330678326e-08)
  !< [Simple]
  w = 100._dp

  !  Get estimation of error and number of evaluations
  call qags(fquad453, Zero, 1._dp, Integ1, gkrule='qk15', abserr=err, Neval=N)
  print "(A)", 'integrate( \log(x)/\sqrt(x), 0, 1) = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)
  ! integrate( \log(x)/\sqrt(x), 0, 1) = -3.9999999999999 (Error=5.1247894816697e-13) with N = 225

  call qags(fquad455w, Zero, nf_inf, [w], Integ1, gkrule='qk21', abserr=err, Neval=N)
  print "(A)", 'integrate( \log(x)/(1 + w x^2), 0, +∞) = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)
  ! integrate( \log(x)/(1 + w x^2), 0, +∞) = -0.3616892188416 (Error=2.170181735946e-07) with N = 399
contains

  !> \f$\log(x)/\sqrt{(x)} \f$
  function fquad453(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = log(x) / sqrt(x)
  end function fquad453

  !> \f$\log(x)/(1 + w x^2) \f$. Result(w=100) = -0.3616892
  function fquad455w(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = log(x) / (1 + w(1) * x**2)
  end function fquad455w
end program ex_qags
