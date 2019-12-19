program ex_trapz
  USE numfor, only: dp, Zero, M_PI, str, trapz, linspace
  implicit none

  real(dp), dimension(20) :: x
  real(dp), dimension(20) :: y
  !< [Simple]
  real(dp) :: Integ1
  intrinsic dsin

  Integ1 = trapz(dsin, Zero, M_PI)
  print "(A)", '\int \sin(x) dx = '//str(Integ1)//" (Difference="//str(abs(2 - Integ1))//")"
  ! \int \sin(x) dx = 1.999832163894 (Difference=0.000167836106)
  !< [Simple]
  !  Integration of sampled data
  x = linspace(Zero, M_PI, 20)
  y = dsin(x)
  Integ1 = trapz(y, Zero, M_PI)
  print "(A)", '\int \sin(x) dx = '//str(Integ1)//" (Difference="//str(abs(2 - Integ1))//")"
  ! \int \sin(x) dx = 1.9980436909706 (Difference=0.0019563090294)

  ! Integration of a function
  Integ1 = trapz(dsin, Zero, M_PI, N=size(y))
  print "(A)", '\int \sin(x) dx = '//str(Integ1)//" (Difference="//str(abs(2 - Integ1))//")"
  ! \int \sin(x) dx = 1.9980436909706 (Difference=0.0019563090294)

end program ex_trapz

