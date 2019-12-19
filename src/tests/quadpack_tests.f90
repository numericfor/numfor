!> The routine called *fquadxyz* computes the function from section x.y.z of Quadpack, Piessens et al (1983)
module testquadf
contains
  !> \f$\sqrt(x) \log(x) \f$
  function fquad451(x) result(y)
    implicit none
    real(dp) :: y !< function value
    real(dp), intent(IN) :: x !< input
    y = sqrt(x) * log(x)
  end function fquad451

  !> \f$\cos(100 \sin(x)) \f$
  function fquad452(x) result(y)
    implicit none
    real(dp) :: y !< function value
    real(dp), intent(IN) :: x !< input
    y = cos(100 * sin(x))
  end function fquad452

  !> funct1 Computes
  function funct1(x) result(y)
    implicit none
    real(dp) :: y !< function value
    real(dp), intent(IN) :: x !< input
    y = sin(x) * exp(-x) * log(1 + x)
  end function funct1
end module testquadf
