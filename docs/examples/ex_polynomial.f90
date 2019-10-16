program ex_polynomial
  USE numfor, only: dp, Zero, M_PI, center, arange, str
  USE numfor, only: polyval, polyder, polyint
  implicit none

  real(dp), dimension(5) :: p1

  p1 = arange(5, 1, -1)

  print "(A)", center(' Polynomials ', 40, '*')
  print "(A)", 'Coefficients:'//str(p1)
  print "(A)", center('Evaluations', 40, '-')
  !< [evaluate]
  print "(A)", "P(-1)= "//str(polyval(p1, -1._dp)) ! P(-1)= 3
  print "(A)", "P(0) = "//str(polyval(p1, 0._dp))  ! P(0) = 1
  print "(A)", "P(1) = "//str(polyval(p1, 1._dp))  ! P(1) = 15
  print "(A)", "P([-1,0,1]) = "//str(polyval(p1, [-1._dp, 0._dp, 1._dp]))
  ! P([-1,0,1]) = [ 3, 1, 15]
  !< [evaluate]
  print "(A)", center('Derivatives', 40, '-')
  !< [derivative]
  print "(A)", "Coefficients of P^(1)(x) = "//str(polyder(p1, 1))
  ! Coefficients of P^(1)(x) = [ 20, 12, 6, 2]
  print "(A)", "Coefficients of P^(2)(x) = "//str(polyder(p1, 2))
  ! Coefficients of P^(2)(x) = [ 60, 24, 6]
  print "(A)", "Coefficients of P^(3)(x) = "//str(polyder(p1, 3))
  ! Coefficients of P^(3)(x) = [120, 24]
  print "(A)", "P'([-1,0,1]) = "//str(polyval(polyder(p1, 1), [-1._dp, 0._dp, 1._dp]))
  ! P'([-1,0,1]) = [-12, 2, 40]
  !< [derivative]
  print "(A)", center('Integration', 40, '-')
  !< [integrate]
  print "(A)", "Coefficients of first antiderivative  (cte=0) ="//str(polyint(p1, 1))
  ! Coefficients of first antiderivative  (cte=0) =[ 1, 1, 1, 1, 1, 0]
  print "(A)", "Coefficients of first antiderivative  (cte=3) ="//str(polyint(p1, 1, 3._dp))
  ! Coefficients of first antiderivative  (cte=3) =[ 1, 1, 1, 1, 1, 3]
  print "(A)", "Coefficients of double antiderivative (ctes=0)="//str(polyint(p1, 2))
  ! Coefficients of double antiderivative (ctes=0)=[ 0.1666666666667, 0.2, 0.25, 0.3333333333333, 0.5, 0, 0]
  print "(A)", "(\int P(x) dx)([-1,0,1]) = "//str(polyval(polyint(p1, 1), [-1._dp, 0._dp, 1._dp]))
  ! (\int P(x) dx)([-1,0,1]) = [ -1, 0, 5]
  !< [integrate]
end program ex_polynomial

