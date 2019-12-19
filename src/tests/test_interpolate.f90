!> test_interpolate
program test_interpolate
  USE utils, only: dp, Zero, M_PI
  USE arrays, only: linspace, geomspace, allclose
  USE interpolate
  implicit none
  call test_splines_1()
  ! call test_splines_2()

contains
  subroutine test_splines_1()
    integer, parameter :: Ndim = 20, Ndimnew = 4 * Ndim
    type(CubicSpline) :: csp
    type(CubicSpline) :: csd
    real(dp), dimension(Ndim) :: x, y, Err
    real(dp), dimension(Ndimnew) :: xnew, ynew, yint, yder
    integer :: i1

    x = linspace(1.e-6_dp, 20._dp, Ndim)
    y = sin(x)
    xnew = linspace(1.e-6_dp, 19.99_dp, Ndimnew)

    ! Notice that we know exactly the second derivative
    ! call csplrep(x, y, -y(1), -y(Ndim), csp)
    csp = CubicSpline(x, y, -y(1), -y(Ndim))
    ynew = csp%evaluate(xnew)
    print *, 'single value', csp%evaluate(1.24_dp)

    ! First-order derivative
    csd = csplder(csp, 1)
    yder = csd%evaluate(xnew)
    print *, 'All derivatives equal', allclose(yder, csp%derivative(xnew))

    write (*, "(A)") ' '//repeat("*-*", 21)
    write (*, '(A)') '# Test 1 of splines module'
    write (*, '(A)') ' function: y= sin(x) '
    write (*, "(A)") repeat("*", 65)

    open (UNIT=9, file='test_splines1.dat')
    write (9, '(A)') '#    x          y=sin(x)        spl(y)      I=-cos(x) Integral(spl(y))'
    do i1 = 1, Ndimnew
      ! yint(i1) = -1 + csplint(xnew(1), xnew(i1), csp) ! Integrate
      yint(i1) = -1 + csp%integrate(xnew(1), xnew(i1)) ! Integrate
      write (9, '(6(ES13.5,1x))') xnew(i1), sin(xnew(i1)), ynew(i1), -cos(xnew(i1)), yint(i1), -yder(i1)
    end do
    close (9)

    write (*, "((f0.5,1x))") csplint_square(Zero, M_PI, csp)
    call spleps(x, y, Err)
    open (UNIT=9, file='error_splines1.dat')
    write (9, '(3(f10.5,1x))') (x(i1), y(i1), abs(Err(i1)), i1=2, Ndim - 1)
    close (9)
  end subroutine test_splines_1

  subroutine test_splines_2()
    integer, parameter :: Ndim = 150, Ndimnew = 4 * Ndim
    type(CubicSpline) :: tck
    real(dp), dimension(Ndim) :: x, y
    real(dp), dimension(Ndim) :: ylog
    real(dp), dimension(Ndim) :: Err
    real(dp), dimension(Ndimnew) :: xnew, ynew, ylognew
    integer :: i1
    integer ::  NN
    real(dp) :: frac, xmed, last

    frac = 0.9_8
    NN = Ndim / 2
    last = 15.65_dp
    xmed = frac * last
    x(:NN) = linspace(1.e-4_dp, xmed, NN, endpoint=.False.)
    ! x(:NN) = geomspace(1.e-2_dp, xmed, NN, endpoint=.False.)
    x(NN + 1:) = linspace(xmed, last, Ndim - NN)

    y = f(x)
    ylog = log(y)

    ! print "(6(g0.4,1x))", x
    ! print "(6(g0.4,1x))", y

    xnew = linspace(0.015_dp, 0.98_dp * last, Ndimnew)

    call csplrep(x, y, Zero, Zero, tck)
    ynew = csplev(xnew, tck)

    call csplrep(x, ylog, Zero, Zero, tck)
    ylognew = csplev(xnew, tck)

    write (*, "(A)") repeat("*", 65)
    write (*, '(A)') '# Test 2 of splines module'
    write (*, '(A)') ' function: y= x**2/sqrt(2.*(-0.58 + (9/x)*(1 - x + x**2 + x**9)))'
    write (*, "(A)") repeat("*", 65)
    open (UNIT=9, file='test_splines2.dat')
    write (9, *) '#    x         y= V(x)      spl(y)     exp(spl(log(y)))'
    do i1 = 1, Ndimnew
      write (9, '(4(1PE13.5,1x))') xnew(i1), f(xnew(i1)), ynew(i1), exp(ylognew(i1))
    end do
    close (9)

    call spleps(x, y, Err)
    open (UNIT=9, file='error_splines2.dat')
    write (9, '(2(1PE13.5,1x))') (x(i1), abs(Err(i1)), i1=1, Ndim)
    close (9)
  end subroutine test_splines_2

  elemental function f(x1) result(y1)
    real(dp), intent(IN) :: x1
    real(dp) :: y1
    y1 = abs(sin(x1))
    ! y1 = x1**2 / sqrt(2.*(-0.58 + (9._8 / x1) * (1._8 - x1 + x1**2 + x1**9)))
  end function f

end program test_interpolate
