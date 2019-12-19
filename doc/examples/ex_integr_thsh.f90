program ex_qthsh1
  USE numfor, only: dp, Zero, M_PI, str, center, qnthsh

  real(dp) :: err
  integer :: N
  !< [Simple]
  real(dp) :: Integ1
  intrinsic dsin
  print "(A)", center(" Integrate sin(x) between 0 and pi ", 70, '-')
  call qnthsh(dsin, Zero, M_PI, Integ1)
  print "(A)", '\int \sin(x) dx = '//str(Integ1)//" (Difference="//str(abs(2 - Integ1))//")"
  ! that prints:
  !------------------ Integrate sin(x) between 0 and pi -----------------
  ! \int \sin(x) dx = 2 (Difference=0)
  !< [Simple]
  print "(A)", center(" Get estimation of error and number of evaluations ", 70, '-')
  call qnthsh(dsin, Zero, M_PI, Integ1, abserr=err, Neval=N)
  print "(A)", '\int \sin(x) dx = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)
  ! that prints:
  !---------- Get estimation of error and number of evaluations ---------
  ! \int \sin(x) dx = 2 (Error=4.8361314952672e-12) with N = 65
end program ex_qthsh1
