subroutine fpcuro(a, b, c, d, x, n)
  !  subroutine fpcuro finds the real zeros of a cubic polynomial
  !  p(x) = a*x**3+b*x**2+c*x+d.
  !
  !  calling sequence:
  !     call fpcuro(a,b,c,d,x,n)
  !
  !  input parameters:
  !    a,b,c,d: real values, containing the coefficients of p(x).
  !
  !  output parameters:
  !    x      : real array,length 3, which contains the real zeros of p(x)
  !    n      : integer, giving the number of real zeros of p(x).
  !  ..
  !  ..scalar arguments..
  real(8) :: a, b, c, d
  integer :: n
  !  ..array argument..
  real(8) :: x(3)
  !  ..local scalars..
  integer :: i
  real(8) :: a1, b1, c1, df, disc, d1, e3, f, four, half, ovfl, pi3, p3, q, r
  real(8) :: step, tent, three, two, u, u1, u2, y
  !  ..function references..
  real(8) :: abs, max, datan, atan2, cos, sign, sqrt
  !  set constants
  two = 2._8
  three = 3._8
  four = 4._8
  ovfl = 1.e+04_8
  half = 0.5_8
  tent = 0.1_8
  e3 = tent / 0.3_8
  pi3 = datan(1._8) / 0.75_8
  a1 = abs(a)
  b1 = abs(b)
  c1 = abs(c)
  d1 = abs(d)

  if (max(b1, c1, d1) < a1 * ovfl) then
    !  p(x) is a third degree polynomial.
    b1 = b / a * e3
    c1 = c / a
    d1 = d / a
    q = c1 * e3 - b1 * b1
    r = b1 * b1 * b1 + (d1 - b1 * c1) * half
    disc = q * q * q + r * r
    if (disc > 0._8) then
      u = sqrt(disc)
      u1 = -r + u
      u2 = -r - u
      n = 1
      x(1) = sign(abs(u1)**e3, u1) + sign(abs(u2)**e3, u2) - b1
    else
      u = sqrt(abs(q))
      if (r < 0.) u = -u
      p3 = atan2(sqrt(-disc), abs(r)) * e3
      u2 = u + u
      n = 3
      x(1) = -u2 * cos(p3) - b1
      x(2) = u2 * cos(pi3 - p3) - b1
      x(3) = u2 * cos(pi3 + p3) - b1
    end if

  else if (max(c1, d1) < b1 * ovfl) then
    !   p(x) is a second degree polynomial.
    disc = c * c - four * b * d
    n = 0
    if (disc < 0._8) return
    n = 2
    u = sqrt(disc)
    b1 = b + b
    x(1) = (-c + u) / b1
    x(2) = (-c - u) / b1

  else if (d1 < c1 * ovfl) then
    !   p(x) is a first degree polynomial.
    n = 1
    x(1) = -d / c
  else
    !  p(x) is a constant function.
    n = 0
    return
  end if

  !
  !  apply a newton iteration to improve the accuracy of the roots.
  do i = 1, n
    y = x(i)
    f = ((a * y + b) * y + c) * y + d
    df = (three * a * y + two * b) * y + c
    step = 0.
    if (abs(f) < abs(df) * tent) step = f / df
    x(i) = y - step
  end do

end subroutine fpcuro
