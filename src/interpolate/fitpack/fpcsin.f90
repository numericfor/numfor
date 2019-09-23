subroutine fpcsin(a, b, par, sia, coa, sib, cob, ress, resc)
  !  fpcsin calculates the integrals ress=integral((b-x)**3*sin(par*x))
  !  and resc=integral((b-x)**3*cos(par*x)) over the interval (a,b),
  !  given sia=sin(par*a),coa=cos(par*a),sib=sin(par*b) and cob=cos(par*b)
  !  ..
  !  ..scalar arguments..
  real(8) :: a, b, par, sia, coa, sib, cob, ress, resc
  !  ..local scalars..
  integer :: i, j
  real(8) :: ab, ab4, ai, alfa, beta, b2, b4, eps, fac, f1, f2
  real(8) :: one, quart, six, three, two
  !  ..function references..
  real(8) :: abs
  !  ..
  one = 1._8
  two = 2._8
  three = 3._8
  six = 6._8
  quart = 0.25_8
  eps = 0.1e-09_8
  ab = b - a
  ab4 = ab**4
  alfa = ab * par
  ! the way of calculating the integrals ress and resc depends on
  ! the value of alfa = (b-a)*par.
  if (abs(alfa) <= one) then
    ! ress and resc are found by evaluating a series expansion.
    fac = quart
    f1 = fac
    f2 = 0.
    i = 4
    do j = 1, 5
      i = i + 1
      ai = i
      fac = fac * alfa / ai
      f2 = f2 + fac
      if (abs(fac) <= eps) exit
      i = i + 1
      ai = i
      fac = -fac * alfa / ai
      f1 = f1 + fac
      if (abs(fac) <= eps) exit
    end do

    ress = ab4 * (coa * f2 + sia * f1)
    resc = ab4 * (coa * f1 - sia * f2)
  else
    ! integration by parts.
    beta = one / alfa
    b2 = beta**2
    b4 = six * b2**2
    f1 = three * b2 * (one - two * b2)
    f2 = beta * (one - six * b2)
    ress = ab4 * (coa * f2 + sia * f1 + sib * b4)
    resc = ab4 * (coa * f1 - sia * f2 + cob * b4)
  end if

end subroutine fpcsin
