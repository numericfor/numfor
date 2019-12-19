program test_quadrature
  use utils, only: dp, M_PI_2, timer, str
  use arrays, only: linspace
  use integrate, only: trapz, simps, qags, qag, qk31
  ! use integrate, only: trapz, simps, qags!, d_trapz_ab_func
  implicit none

  integer, parameter :: N = 20
  real(dp), dimension(N) :: x
  real(dp), dimension(N) :: f
  real(dp) :: Int1, Int2
  real(dp) :: a, b
  real(dp) :: epsabs, epsrel, abserr
  integer :: i, Neval
  integer :: ier
  type(timer) T1
  intrinsic ::  dcos

  abstract interface
    function func(z)
      real(8) :: func
      real(8), intent(in) :: z
    end function func
  end interface
  procedure(func), pointer :: f_ptr => null()

  call T1%start()
  a = 0._dp; b = 1._dp
  x = linspace(a, b, N)
  do i = 1, N
    f(i) = funct1(x(i))
  end do

  do i = 1, 7
    Int1 = trapz(funct1, a, b, i * i * N)
    print "(A)", 'trapz = '//str(Int1)//" with N = "//str(i * i * N)
    Int2 = simps(funct1, a, b, i * i * N)
    print "(A)", 'simps = '//str(Int2)//" with N = "//str(i * i * N)
  end do

  epsabs = 1.e-8_8
  epsrel = 1.e-16_8
  call qag(funct1, a, b, epsabs, epsrel, 1, Int1, abserr, Neval, ier)
  print "(A)", 'qag  = '//str(Int1)//" with N="//str(Neval)//",  abserr="//str(abserr)
  call qag(funct1, a, b, epsabs, epsrel, 2, Int1, abserr, Neval, ier)
  print "(A)", 'qag  = '//str(Int1)//" with N="//str(Neval)//",  abserr="//str(abserr)
  call qags(funct1, a, b, epsabs, epsrel, 'qk15', Int1, abserr, Neval, ier)
  print "(A)", 'qags  = '//str(Int1)//" with N="//str(Neval)//",  abserr="//str(abserr)
  call qags(funct1, a, b, epsabs, epsrel, Int1, abserr, Neval, ier)
  print "(A)", 'qags  = '//str(Int1)//" with N="//str(Neval)//",  abserr="//str(abserr)
  call qags(funct1, a, b, epsabs, epsrel, 'qk31', Int1, abserr, Neval, ier)
  print "(A)", 'qags  = '//str(Int1)//" with N="//str(Neval)//",  abserr="//str(abserr)

  Int1 = trapz(f, x)
  Int2 = trapz(f, a, b)
  print "(A)", 'sample trapz = '//str(Int1)//" with N="//str(size(x))
  print "(A)", 'sample trapz = '//str(Int2)//" with N="//str(size(x))

  f_ptr => dcos
  Int2 = trapz(f_ptr, 0._dp, M_PI_2, N)
  print "(A)", 'trapz (cos(x)) = '//str(Int2)//" with N="//str(N)
  Int2 = trapz(f_ptr, 0._dp, M_PI_2, 5 * N)
  print "(A)", 'trapz (cos(x)) = '//str(Int2)//" with N="//str(5 * N)
  Int2 = trapz(f_ptr, 0._dp, M_PI_2, 10 * N)
  print "(A)", 'trapz (cos(x)) = '//str(Int2)//" with N="//str(10 * N)

  call T1%stop()
  call T1%show()

contains
  !> funct1 Computes
  !!
  function funct1(x) result(y)
    implicit none
    real(dp) :: y !<
    real(dp), intent(IN) :: x !<
    !! Examples:
    !!
    y = sin(x) * exp(-x) * log(1 + x)
  end function funct1

end program test_quadrature
