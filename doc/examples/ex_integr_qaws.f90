program ex_qaws
  USE numfor, only: dp, Zero, str, qaws
  implicit none

  real(dp) :: err
  integer :: N
  !< [Simple]
  real(dp) :: Integ1
  real(dp) :: alfa, beta
  integer :: flgw
  alfa = 0._dp
  beta = 0._dp
  flgw = 2                 ! weight function => log(x-a)^alfa = log(x)
  ! fquad458 = 1/(1 + \ln(x)^{2})^{2}
  call qaws(fquad458, Zero, 1._dp, alfa, beta, flgw, Integ1)
  print "(A)", 'integrate(log(x)/(1 + ln(x)^{2})^{2}, 0, 1) = '//str(Integ1)
  ! integrate(log(x)/(1 + ln(x)^{2})^{2}, 0, 1) = -0.1892750041577
  !< [Simple]

  !  Get estimation of error and number of evaluations
  call qaws(fquad458, Zero, 1._dp, alfa, beta, flgw, Integ1, epsrel=1.e-3_dp, abserr=err, Neval=N)
  print "(A)", 'integrate(log(x)/(1 + ln(x)^{2})^{2}, 0, 1) = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)
  ! integrate(log(x)/(1 + ln(x)^{2})^{2}, 0, 1) = -0.1892747444915 (Error=2.2226716553438e-06) with N = 40
contains

  !> \f$1/(1 + \ln(x)^{2})^{2} \f$
  function fquad458(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = 1 / (1 + log(x)**2)**2
  end function fquad458

end program ex_qaws
