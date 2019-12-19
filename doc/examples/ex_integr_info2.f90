program ex_qp_extra
  USE numfor, only: dp, Zero, M_PI_2, str, d_qp_extra, qagp, qags
  implicit none

  real(dp) :: err
  integer :: ier, N
  integer :: j, k
  type(d_qp_extra) :: info
  real(dp) :: Integ1
  real(dp), dimension(1) :: pts

  ! Initialize allowing up to 50 subintervals
  ! This is optional, by default it will allow for 500 subintervals
  info = 50

 call qagp(fquad4510, Zero, 1._dp, [sqrt(3._dp) - 1], Integ1, epsabs=Zero, epsrel=1.e-4_dp, abserr=err, Neval=N, ier=ier, info=info)

  print "(A)", 'integrate(1/sqrt(|x^2 + 2x -2|), x, 0, 1) = '//str(Integ1)//" (Error="//str(err)//") with N = "//str(N)
  print "(A)", "Difference = "//str(abs(Integ1 - (M_PI_2 - asin(1 / sqrt(3._dp)) + log(3._dp) / 2)))//&
    & " relative to analytical value"
  print "(A)", "Error code = "//str(ier)
  print "(A)", "Meaning:"
  print "(A)", "-------"
  print "(A)", info%msg
  print "(A)", ""
  print "(A)", " j        Limits of subinterval S(j)     Integral over S(j)    Error on S(j)    Error relative "
  do k = 1, info%last
    j = info%iord(k)
    print "(I2, 5(es18.11,1x))", k, info%alist(j), info%blist(j), info%rlist(j), info%elist(j), info%elist(j) / info%rlist(j)
  end do

contains

  !> \f$1/|x^2 - 2x -2| \f$. Result = 1.504622
  function fquad4510(x) result(y)
    implicit none
    real(dp) :: y
    real(dp), intent(IN) :: x
    y = 1._dp / sqrt(abs(x**2 + 2 * x - 2))
  end function fquad4510
end program ex_qp_extra

