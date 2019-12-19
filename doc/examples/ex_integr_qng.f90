program ex_qng
  USE numfor, only: dp, Zero, M_PI, str, center, qng
  implicit none

  real(dp) :: err
  integer :: N
  real(dp) :: w
  !< [Simple]
  real(dp) :: Integ1
  complex(dp) :: Integ2
  intrinsic dsin
  print "(A)", center(" Integrate sin(x) between 0 and pi ", 70, '-')
  call qng(dsin, Zero, M_PI, Integ1)
  print "(A)", '\int \sin(x) dx = '//str(Integ1)//" (Difference="//str(abs(2 - Integ1))//")"
  ! that prints:
  !------------------ Integrate sin(x) between 0 and pi -----------------
  ! \int \sin(x) dx = 2 (Difference=0)
  !< [Simple]
  w = 100._dp
  print "(A)", center(" Get estimation of error and number of evaluations ", 70, '-')
  call qng(fquad451, Zero, M_PI, Integ1, abserr=err, Neval=N)
  print "(A)", '\int \sqrt(x) \log(x) dx = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)
  ! that prints:
  !---------- Get estimation of error and number of evaluations ---------
  ! \int \sqrt(x) \log(x) dx = 1.774674495015 (Error=6.9889313343983e-05) with N = 87
  print "(A)", center(" Integrate a complex function ", 70, '-')
  call qng(fquad452wc, Zero, M_PI, [w], Integ2, abserr=err, Neval=N)
  print "(A)", '\int \cos(w \sin(x)) dx = '//str(Integ2)//" (Error="//str(err)//") with N = "//str(N)
  ! that prints:
  !-------------------- Integrate a complex function --------------------
  ! \int \cos(w \sin(x)) dx = (0.062787400402-0.2226721656122j) (Error=2.8652541230259) with N = 87

contains
  !> \f$\sqrt(x) \log(x) \f$
  function fquad451(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = sqrt(x) * log(x)
  end function fquad451

  !> \f$\cos(w \sin(x)) \f$
  function fquad452wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = cmplx(cos(w(1) * sin(x)), sin(w(1) * sin(x)), kind=dp)
  end function fquad452wc

end program ex_qng
