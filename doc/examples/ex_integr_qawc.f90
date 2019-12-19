program ex_qaws
  USE numfor, only: dp, str, qawc, nf_minf
  implicit none
  real(dp) :: err
  integer :: N, ier
  !< [Simple]
  real(dp) :: Integ1
  call qawc(fquad459, -1._dp, 5._dp, 0._dp, Integ1)
  print "(A)", 'integrate(1/(x(5 x^3 + 6), -1, 5))  = '//str(Integ1)
  ! integrate(1/(x(5 x^3 + 6), -1, 5))  = -0.0899440069576
  !< [Simple]

  call qawc(fquad459, -1._dp, 5._dp, 0._dp, Integ1, epsrel=1.e-3_dp, abserr=err, Neval=N, ier=ier)
  print "(A)", 'integrate(1/(x(5 x^3 + 6), -1, 5)) = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)
  ! integrate(1/(x(5 x^3 + 6), -1, 5)) = -0.0899440069584 (Error=1.1852922354142e-06) with N = 215
contains

  !> Example for QAWC \f$ \frac{1}{x (w(1) x^3+ (w(2)))} \f$. Result = -0.08994401
  function fquad459(x) result(y)
    implicit none
    real(dp) :: y !< function value
    real(dp), intent(IN) :: x !< variable
    y = 1 / (5 * x**3 + 6)
  end function fquad459

end program ex_qaws
