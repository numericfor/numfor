!> test_interpolate
program test_interpolate
  USE utils, only: dp, Zero, M_PI
  USE arrays, only: linspace
  USE interpolate
  implicit none
  call test_splines_1()
  call test_splines_2()

contains
  subroutine test_splines_1()
    integer, parameter :: Ndim = 20, Ndimnew = 4 * Ndim
    type(cspl_rep) :: tck
    real(dp), dimension(Ndim) :: x, y, Err
    real(dp), dimension(Ndimnew) :: xnew, ynew, yint
    integer :: i1

    x = linspace(Zero, 20._dp, Ndim)
    y = sin(x)
    xnew = linspace(1.e-6_dp, 19.99_dp, Ndimnew)

    call csplrep(x, y, Zero, Zero, tck)
    ynew = csplev(xnew, tck)
    write (*, "(A)") ' '//repeat("*-*", 21)
    write (*, '(A)') '# Test 1 of splines module'
    write (*, '(A)') ' function: y= sin(x) '
    write (*, "(A)") repeat("*", 65)
    open (UNIT=9, file='test_splines1.dat')
    write (9, '(A)') '#    x          y=sin(x)        spl(y)      I=-cos(x) Integral(spl(y))'
    do i1 = 1, Ndimnew
      yint(i1) = -1 + splint(xnew(1), xnew(i1), tck) ! Integrate
      write (9, '(5(ES13.5,1x))') xnew(i1), sin(xnew(i1)), ynew(i1), -cos(xnew(i1)), yint(i1)
    end do
    close (9)

    write (*, "((f0.5,1x))") splint_square(Zero, M_PI, tck)
    call spleps(x, y, Err)
    open (UNIT=9, file='error_splines1.dat')
    write (9, '(3(f10.5,1x))') (x(i1), y(i1), abs(Err(i1)), i1=2, Ndim - 1)
    close (9)
  end subroutine test_splines_1

  subroutine test_splines_2()
    integer, parameter :: Ndim = 250, Ndimnew = 2 * Ndim
    type(cspl_rep) :: tck
    real(dp), dimension(Ndim) :: x, y
    real(dp), dimension(Ndim) :: ylog
    real(dp), dimension(Ndim) :: Err
    real(dp), dimension(Ndimnew) :: xnew, ynew, ylognew
    integer :: i1
    integer ::  NN
    real(dp) :: Dx, Dxn
    real(dp) :: frac, xmed

    frac = 0.9_8
    NN = Ndim / 2
    xmed = 15.65 * frac
    Dx = (xmed - Zero) / (NN - 1)
    do i1 = 1, NN
      x(i1) = Dx * (i1 - 1) + 0.0001_8
      y(i1) = f(x(i1))
    end do
    Dx = (15.65 - xmed) / (NN)
    do i1 = NN, Ndim
      x(i1) = xmed + Dx * (i1 - NN) + 0.0001_8
      y(i1) = f(x(i1))
    end do

    ylog = log(y)
    Dxn = 15.65*.98 / (Ndimnew - 1)
    do i1 = 1, Ndimnew
      xnew(i1) = Dxn * (i1 - 1) + 0.015_8
    end do

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

  function f(x1) result(y1)
    real(dp), intent(IN) :: x1
    real(dp) :: y1
    y1 = x1**2 / sqrt(2.*(-.58 + (9._8 / x1) * (1._8 + exp(-x1 / 3.5))))
  end function f

end program test_interpolate
