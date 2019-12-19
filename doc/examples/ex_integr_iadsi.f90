program ex_iadsi
  USE numfor, only: dp, Zero, M_PI, str, iadsi, nf_inf, qags
  implicit none

  real(dp) :: err
  integer :: N

  real(dp) :: Integ1
  intrinsic dsin

  call iadsi(f, Zero, [1._dp], Integ1, epsrel=1.e-7_dp, abserr=err, Neval=N)
  print "(A)", 'integrate(log(1+x)/(1 + 100*x^2), 0, inf) = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(abs(N))
  ! integrate(log(1+x)/(1 + 100*x^2), 0, inf) = 0.0335216255386 (Error=1.259162853146e-12) with N = 165
contains

  !> \f$\log(x)/(1 + x^2) \f$.
  function f(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = log(1 + x) / (1 + 100 * x**2)
  end function f

end program ex_iadsi

