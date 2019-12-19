function fprati(p1, f1, p2, f2, p3, f3) result(y)
  !  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
  !  gives the value of p such that the rational interpolating function
  !  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
  !  ..
  !  ..scalar arguments..
  real(8) :: p1, f1, p2, f2, p3, f3
  real(8) :: y
  !  ..local scalars..
  real(8) :: h1, h2, h3, p
  !  ..
  if (p3 <= 0._8) then
    !  value of p in case p3 = infinity.
    p = (p1 * (f1 - f3) * f2 - p2 * (f2 - f3) * f1) / ((f1 - f2) * f3)
  else
    !  value of p in case p3 ^= infinity.
    h1 = f1 * (f2 - f3)
    h2 = f2 * (f3 - f1)
    h3 = f3 * (f1 - f2)
    p = -(p1 * p2 * h3 + p2 * p3 * h1 + p3 * p1 * h2) / (p1 * h1 + p2 * h2 + p3 * h3)
  end if

  !  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  if (f2 >= 0._8) then
    p1 = p2
    f1 = f2
  else
    p3 = p2
    f3 = f2
  end if
  y = p
end function fprati
