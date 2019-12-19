program ex_qagp
  USE numfor, only: dp, Zero, str, qagp
  implicit none

  real(dp) :: err
  integer :: N
  real(dp), parameter :: epsabs = 0.0_8
  real(dp), parameter :: epsrel = 1.e-3_8

  !< [Simple]
  real(dp), dimension(2) :: points
  real(dp) :: Integ1

  points(:2) = [1._8, sqrt(2._8)]

  call qagp(fquad454, Zero, 3._dp, points, Integ1, epsabs, epsrel, err, Neval=N)
  print "(A)", 'integrate(x^3 \log{|(x^2-1)(x^2-2)|}, x, 0, 3) = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)
  ! integrate(x^3 \log{|(x^2-1)(x^2-2)|}, x, 0, 3) = 52.7408058906746 (Error=0.0001754605945) with N = 777

contains

  !> Example for QAGP \f$\log{|(x^2-1)(x^2-2)|} \f$. Result = 52.740748
  function fquad454(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = x**3 * log(abs((x**2 - 1._dp) * (x**2 - 2._dp)))
  end function fquad454

end program ex_qagp
