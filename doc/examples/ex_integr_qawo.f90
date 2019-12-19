program ex_qawo
  USE numfor, only: dp, Zero, M_PI, str, qawo
  implicit none
  real(dp) :: err
  integer :: N
  !< [Simple]
  real(dp) :: Integ1
  ! Fifth argument, flgw=1 => weight function is cosine
  ! fquad456 => log(x)
  call qawo(fquad456, Zero, 1._dp, 10 * M_PI, 2, Integ1)
  print "(A)", 'integrate(sin(10 pi x) log(x), x, 0, 1) = '//str(Integ1)
  ! integrate(sin(10 pi x) log(x), x, 0, 1) = -0.1281368491001
  !< [Simple]

  call qawo(fquad456, Zero, 1._dp, 10 * M_PI, 1, Integ1, epsabs=Zero, epsrel=1.e-3_dp, abserr=err, Neval=N)

  print "(A)", 'integrate(cos(10 pi x) log(x), x, 0, 1) = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)
  ! integrate(cos(10 pi x) log(x), x, 0, 1) = -0.0489888171135 (Error=4.8212107917056e-10) with N = 305
contains
  !> Example for QAWO \f$\log(x) \f$. Result = -0.12181316
  function fquad456(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = Zero
    IF (x > 0) y = log(x)
  end function fquad456

end program ex_qawo
