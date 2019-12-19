program ex_iads
  USE numfor, only: dp, Zero, M_PI, str, iads
  implicit none

  real(dp) :: err
  integer :: N

  !< [Simple]
  real(dp) :: Integ1
  intrinsic dsin

  call iads(dsin, Zero, M_PI, Integ1)
  print "(A)", 'integrate(sin(x), 0, pi) = '//str(Integ1)//" (Difference="//str(abs(2 - Integ1))//")"
  ! integrate(sin(x), 0, pi) = 1.9999998386284 (Difference=1.6137164227104e-07)
  ! < [Simple]
  call iads(fquad451, 1.e-12_dp, M_PI, Integ1, epsrel=1.e-7_dp, abserr=err, Neval=N)
  print "(A)", 'integrate( \log(x) \sqrt(x), 1.e-12, pi) = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(abs(N))
  ! integrate( \log(x) \sqrt(x), 1.e-12, pi) = 1.7746775593141 (Error=1.3697319985965e-05) with N = 89

contains
  !> \f$\sqrt(x) \log(x) \f$
  function fquad451(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = sqrt(x) * log(x)
  end function fquad451

end program ex_iads

