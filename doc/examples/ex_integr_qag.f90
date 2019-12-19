program ex_qag
  USE numfor, only: dp, Zero, M_PI, str, center, qag
  implicit none

  real(dp) :: err
  integer :: N
  real(dp) :: w
  !< [Simple]
  real(dp) :: Integ1
  complex(dp) :: Integ2
  real(dp), parameter :: epsrel = 1.e-3_8

  intrinsic dsin
  print "(A)", center(" Integrate sin(x) between 0 and pi ", 70, '-')
  call qag(dsin, Zero, M_PI, Integ1)
  print "(A)", '\int \sin(x) dx = '//str(Integ1)//" (Difference="//str(abs(2 - Integ1))//")"
  ! that prints:
  !------------------ Integrate sin(x) between 0 and pi -----------------
  ! \int \sin(x) dx = 2 (Difference=0)
  !< [Simple]
  w = 100._dp
  print "(A)", center(" Get estimation of error and number of evaluations ", 70, '-')
  call qag(fquad451, Zero, M_PI, Integ1, abserr=err, Neval=N)
  print "(A)", '\int \sqrt(x) \log(x) dx = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)

  call qag(fquad451, Zero, M_PI, Integ1, rule='qk21', abserr=err, Neval=N)
  print "(A)", '\int \sqrt(x) \log(x) dx = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)
  ! ---------- Get estimation of error and number of evaluations ---------
  ! \int \sqrt(x) \log(x) dx = 1.7746752010439 (Error=7.9254688067762e-06) with N = 375
  ! \int \sqrt(x) \log(x) dx = 1.7746751858103 (Error=1.4425452046185e-05) with N = 441

  print "(A)", center(" Integrate a complex function ", 70, '-')
  call qag(fquad452wc, Zero, M_PI, [w], Integ2, rule='qk61', abserr=err, Neval=N)
  print "(A)", '\int \cos(w \sin(x)) dx = '//str(Integ2)//" (Error="//str(err)//") with N = "//str(N)
  ! that prints:
  !-------------------- Integrate a complex function --------------------
  ! \int \cos(w \sin(x)) dx = (0.062787390876-0.2226721893626j) (Error=1.7375840100087e-08) with N = 427

contains
  !> \f$\sqrt(x) \log(x) \f$
  function fquad451(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = sqrt(x) * log(x)
  end function fquad451

  !> \f$\exp(i w \sin(x)) \f$
  function fquad452wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = cmplx(cos(w(1) * sin(x)), sin(w(1) * sin(x)), kind=dp)
  end function fquad452wc

end program ex_qag
