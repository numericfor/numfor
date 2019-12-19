program ex_qawf
  USE numfor, only: dp, Zero, M_PI, str, qawf
  implicit none
  real(dp) :: err
  integer :: N, ier
  !< [Simple]
  real(dp) :: Integ1
  ! Fourth argument, flgw=1 => weight function is cosine
  ! fquad457 => 1/\sqrt(x)
  call qawf(fquad457, Zero, M_PI / 2, 1, Integ1)
  print "(A)", 'integrate(cos(pi x /2 ) / sqrt(x), x, 0, inf ) = '//str(Integ1)
  ! integrate(cos(pi x /2 ) / sqrt(x), x, 0, inf ) = 0.9999999999321
  !< [Simple]

  call qawf(fquad457, Zero, M_PI / 2, 1, Integ1, epsabs=1.e-3_dp, abserr=err, Neval=N, ier=ier)
  print "(A)", 'integrate(cos(pi x /2 ) / sqrt(x), x, 0, inf ) = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)
  ! integrate(cos(pi x /2 ) / sqrt(x), x, 0, inf ) = 0.9999969531148 (Error=0.0005923419178) with N = 380
contains
  !> Example for QAWF \f$ 1/\sqrt(x) \f$. Result = 1
  function fquad457(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = 0._8
    IF (x > 0) y = 1.0_8 / sqrt(x)
  end function fquad457

end program ex_qawf
