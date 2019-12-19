program test_trapz
  use utils, only: dp, M_PI, M_PI_2, timer, str, center
  use arrays, only: linspace
  use integrate, only: trapz, nf_rfunction, nf_cfunction
  implicit none

  integer, parameter :: N = 21
  real(dp) :: Int1, result, err
  complex(dp) :: Int2
  complex(dp), dimension(5) :: results
  real(dp) :: a, b
  real(dp) :: w
  integer :: i
  type(timer) T1
  character(len=:), allocatable :: title
  type(nf_rfunction) :: fun1
  type(nf_cfunction) :: fun2

  ! Results for integrate(exp(i w sin(x)), x, 0, pi ) from maxima
  results = [cmplx((0.3024310191244698, -0.3198921069678047), kind=dp), &
    & cmplx((0.1753395985854673, -0.2680962123081692), kind=dp),&
    & cmplx((0.1088370651016893, -0.2415328394565208), kind=dp), &
    & cmplx((0.06278740049149098, -0.2226721656038109), kind=dp), &
    & cmplx((0.02699336268296957, -0.2065688529554112), kind=dp)]

  call T1%start()

  result = -4._dp / 9
  title = "| integrate( sqrt(x) log(x), x, 0, 1 ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  a = 1.e-5_dp; b = 1._dp
  do i = 1, 5
    Int1 = trapz(fquad451, a, b, i * N)
    err = abs((result - Int1) / result)
    print "(A)", 'trapz = '//str(Int1)//" (Error="//str(err)//") with N = "//str(i * N)
  end do

  result = 0.06278740_8
  title = "| integrate( cos(100 * sin(x)), x, 0, pi ) = "//str(result)//" |"
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  print "(A)", center(title, 60)
  print "(A)", center(repeat('-', len(title) - 2), 60, ' ')
  a = 0._dp; b = M_PI
  do i = 1, 5
    Int1 = trapz(fquad452, a, b, i * N)
    err = abs((result - Int1) / result)
    print "(A)", 'trapz = '//str(Int1)//" (Error="//str(err)//") with N = "//str(i * N)
  end do

  fun1%f => fquad452w
  title = "| integrate( cos(w * sin(x)), x, 0, pi ) = "//str(result)//" |"
  print "(A)", center(" (Using a variable argument:) ", 70, '-')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = 25._dp * i
    result = real(results(i))
    fun1%args = [w]
    Int1 = trapz(fun1, a, b, 4 * N)
    err = abs((result - Int1) / result)
    print "(A)", 'trapz = '//str(Int1)//" (w="//str(w)//", Error="//str(err)//") with N = "//str(4 * N)
  end do

  fun1%f => fquad452ws
  title = "| integrate( sin(w * sin(x)), x, 0, pi ) = "//str(result)//" |"
  print "(A)", center(" (Using a variable argument:) ", 70, '-')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  do i = 1, 5
    w = 25._dp * i
    result = aimag(results(i))
    fun1%args = [w]
    Int1 = trapz(fun1, a, b, 4 * N)
    err = abs((result - Int1) / result)
    print "(A)", 'trapz = '//str(Int1)//" (w="//str(w)//", Error="//str(err)//") with N = "//str(4 * N)
  end do

  fun2%f => fquad452wc
  title = "| integrate( exp( i(w * sin(x)) ), x, 0, pi ) = "//str(result)//" |"
  print "(A)", center(" (Complex, using a variable argument w:) ", 70, '-')
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')
  print "(A)", center(title, 70)
  print "(A)", center(repeat('-', len(title) - 2), 70, ' ')

  do i = 1, 5
    w = 25_dp * i
    fun2%args = [w]
    Int2 = trapz(fun2, a, b, 10 * N)
    err = abs((results(i) - Int2) / results(i))
    print "(A)", "maxima qags = "//str(results(i))//"  w = "//str(w)
    print "(A)", 'Aprox trapz = '//str(Int2)//" Error="//str(err)//"   N = "//str(10 * N)
    ! print *, Int2
  end do

  call T1%stop()
  call T1%show()

contains
  !> \f$\sqrt(x) \log(x) \f$
  function fquad451(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = sqrt(x) * log(x)
  end function fquad451

  !> \f$\cos(100 \sin(x)) \f$
  function fquad452(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    y = cos(100 * sin(x))
  end function fquad452

  !> \f$\cos(w \sin(x)) \f$
  function fquad452wc(x, w) result(y)
    implicit none
    complex(dp) :: y !<
    real(dp), intent(IN) :: x !<
    real(dp), dimension(:), intent(IN) :: w !<
    y = cmplx(cos(w(1) * sin(x)), sin(w(1) * sin(x)), kind=dp)
  end function fquad452wc

  !> \f$\cos(w \sin(x)) \f$
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

end program test_trapz
