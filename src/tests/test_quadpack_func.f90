module quadpack_sample_func
  USE numfor, only: dp, M_PI

  interface fquad452
    module procedure :: fquad452_0, fquad452w, fquad452wc, fquad452ws
  end interface fquad452

  public

contains

  !> Example for QNG \f$\log(x) \sqrt{(x)} \f$. Result = -4/9
  function fquad451(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = sqrt(x) * log(x)
  end function fquad451

  !> Example for QAGS \f$\cos(100 \sin(x)) \f$. Result = J_0(100) = 0.0627874
  function fquad452_0(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = cos(100 * sin(x))
  end function fquad452_0

  !> Example for QNG complex \f$\cos(w \sin(x)) \f$
  function fquad452wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = cmplx(cos(w(1) * sin(x)), sin(w(1) * sin(x)), kind=dp)
  end function fquad452wc
  !> Example for QAG \f$\cos(w \sin(x)) \f$
  function fquad452w(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = cos(w(1) * sin(x))
  end function fquad452w
  !> \f$\cos(w \sin(x)) \f$
  function fquad452ws(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = sin(w(1) * sin(x))
  end function fquad452ws

  !> Example for QAGS \f$\log(x)/\sqrt{(x)} \f$. Result = -4
  function fquad453(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = log(x) / sqrt(x)
  end function fquad453

  !> Example for QAGP \f$x^3 \log{|(x^2-1)(x^2-2)|} \f$. Result = 52.740748
  function fquad454(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = x**3 * log(abs((x**2 - 1._dp) * (x**2 - 2._dp)))
  end function fquad454

  !> Example for QAGP \f$\log{|(x^2-1)(x^2-2)|} \f$. Result = 52.740748
  function fquad454cw(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = cmplx(w(1), w(2)) * x**3 * log(abs((x**2 - 1._dp) * (x**2 - 2._dp)))
  end function fquad454cw

  !> Example for QAGI \f$\log(x)/(1 + w x^2) \f$. Result(w=100) = -0.3616892
  function fquad455(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = log(x) / (1 + 100 * x**2)
  end function fquad455
  function fquad455w(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = log(x) / (1 + w(1) * x**2)
  end function fquad455w

  function fquad455wc(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = log(x) / (1 + cmplx(w(1), w(2), kind=dp) * x**2)
  end function fquad455wc

  !> Example for QAWO \f$\log(x) \sin(10 \pi x) \f$. Result(w=10 pi) = -0.12181316
  function fquad456(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = log(x)
  end function fquad456
  function fquad456w(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = log(w(1) * x)
  end function fquad456w

  !> Example for QAWF \f$\sqrt(x) \cos(\pi x/2) \f$. Result(w=10 pi) = 1
  function fquad457(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = 0._8
    IF (x > 0) y = 1.0_8 / sqrt(x)
  end function fquad457
  function fquad457wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = 0._8
    IF (x > 0) y = cmplx(w(1), w(2), kind=dp) / sqrt(x)
  end function fquad457wc

  !> Example for QAWS \f$\log(x)/(1 + ln(x)^{2})^{2} \f$. Result = 0.1892752
  function fquad458(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = 0._8
    IF (x > 0) y = log(x) / (1 + ln(x)**2)**2
  end function fquad458
  !> Example for QAWS \f$\log(x)/(1 + ln(x)^{2})^{2} \f$. Result = 0.1892752
  function fquad458wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), dimension(:), intent(IN) :: w !<
    real(dp), intent(IN) :: x !<
    y = 0._8
    IF (x > 0) y = cmplx(w(1), w(2), kind=dp) / sqrt(x)
  end function fquad458wc

  !> Example for QAWC \f$ \frac{1}{x (w(1) x^3+ (w(2)))} \f$. Result(w=[5,6]) = -0.08994401
  function fquad459(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = 1 / (5 * x**3 + 6)
  end function fquad459
  function fquad459w(x, w) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = 1 / (w(1) * x**3 + w(2))
  end function fquad459w
  function fquad459wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = 1 / ((5 + w(1)) * x**3 + 6 + cmplx(0._8, 1._8, kind=dp) * w(2))
  end function fquad459wc

end module quadpack_sample_func
