subroutine fpbfou(t, n, par, ress, resc)
  !  subroutine fpbfou calculates the integrals
  !                    /t(n-3)
  !    ress(j) =      !        nj,4(x)*sin(par*x) dx    and
  !              t(4)/
  !                    /t(n-3)
  !    resc(j) =      !        nj,4(x)*cos(par*x) dx ,  j=1,2,...n-4
  !              t(4)/
  !  where nj,4(x) denotes the cubic b-spline defined on the knots
  !  t(j),t(j+1),...,t(j+4).
  !
  !  calling sequence:
  !     call fpbfou(t,n,par,ress,resc)
  !
  !  input parameters:
  !    t    : real array,length n, containing the knots.
  !    n    : integer, containing the number of knots.
  !    par  : real, containing the value of the parameter par.
  !
  !  output parameters:
  !    ress  : real array,length n, containing the integrals ress(j).
  !    resc  : real array,length n, containing the integrals resc(j).
  !
  !  restrictions:
  !    n >= 10, t(4) < t(5) < ... < t(n-4) < t(n-3).
  !  ..
  !  ..scalar arguments..
  integer :: n
  real(8) :: par
  !  ..array arguments..
  real(8) :: t(n), ress(n), resc(n)
  !  ..local scalars..
  integer :: i, ic, ipj, is, j, jj, jp1, jp4, k, li, lj, ll, nmj, nm3, nm7
  real(8) :: ak, beta, con1, con2, c1, c2, delta, eps, fac, f1, f2, f3, one, quart, sign, six, s1, s2, term
  !  ..local arrays..
  real(8) :: co(5), si(5), hs(5), hc(5), rs(3), rc(3)
  !  ..function references..
  real(8) :: cos, sin, abs
  !  ..
  !  initialization.
  one = 1._8
  six = 6._8
  eps = 0.1e-07_8
  quart = 0.25_8
  con1 = 0.5e-01_8
  con2 = 0.12e+03_8
  nm3 = n - 3
  nm7 = n - 7
  if (par /= 0.) term = six / par
  beta = par * t(4)
  co(1) = cos(beta)
  si(1) = sin(beta)
  !  calculate the integrals ress(j) and resc(j), j=1,2,3 by setting up
  !  a divided difference table.
  do j = 1, 3
    jp1 = j + 1
    jp4 = j + 4
    beta = par * t(jp4)
    co(jp1) = cos(beta)
    si(jp1) = sin(beta)
    call fpcsin(t(4), t(jp4), par, si(1), co(1), si(jp1), co(jp1), rs(j), rc(j))
    i = 5 - j
    hs(i) = 0.
    hc(i) = 0.
    do jj = 1, j
      ipj = i + jj
      hs(ipj) = rs(jj)
      hc(ipj) = rc(jj)
    end do

    do jj = 1, 3
      if (i < jj) i = jj
      k = 5
      li = jp4
      do ll = i, 4
        lj = li - jj
        fac = t(li) - t(lj)
        hs(k) = (hs(k) - hs(k - 1)) / fac
        hc(k) = (hc(k) - hc(k - 1)) / fac
        k = k - 1
        li = li - 1
      end do
    end do

    ress(j) = hs(5) - hs(4)
    resc(j) = hc(5) - hc(4)
  end do

  if (nm7 < 4) go to 160
  !  calculate the integrals ress(j) and resc(j),j=4,5,...,n-7.
  do j = 4, nm7
    jp4 = j + 4
    beta = par * t(jp4)
    co(5) = cos(beta)
    si(5) = sin(beta)
    delta = t(jp4) - t(j)
    !  the way of computing ress(j) and resc(j) depends on the value of
    !  beta = par*(t(j+4)-t(j)).
    beta = delta * par
    if (abs(beta) <= one) go to 60
    !  if !beta! > 1 the integrals are calculated by setting up a divided
    !  difference table.
    do k = 1, 5
      hs(k) = si(k)
      hc(k) = co(k)
    end do

    do jj = 1, 3
      k = 5
      li = jp4
      do ll = jj, 4
        lj = li - jj
        fac = par * (t(li) - t(lj))
        hs(k) = (hs(k) - hs(k - 1)) / fac
        hc(k) = (hc(k) - hc(k - 1)) / fac
        k = k - 1
        li = li - 1
      end do
    end do

    s2 = (hs(5) - hs(4)) * term
    c2 = (hc(5) - hc(4)) * term
    go to 130
    !  if !beta! <= 1 the integrals are calculated by evaluating a series
    !  expansion.
60  f3 = 0.
    do i = 1, 4
      ipj = i + j
      hs(i) = par * (t(ipj) - t(j))
      hc(i) = hs(i)
      f3 = f3 + hs(i)
    end do

    f3 = f3 * con1
    c1 = quart
    s1 = f3
    if (abs(f3) <= eps) go to 120
    sign = one
    fac = con2
    k = 5
    is = 0
    do ic = 1, 20
      k = k + 1
      ak = k
      fac = fac * ak
      f1 = 0.
      f3 = 0.
      do i = 1, 4
        f1 = f1 + hc(i)
        f2 = f1 * hs(i)
        hc(i) = f2
        f3 = f3 + f2
      end do

      f3 = f3 * six / fac
      if (is == 0) go to 90
      is = 0
      s1 = s1 + f3 * sign
      go to 100
90    sign = -sign
      is = 1
      c1 = c1 + f3 * sign
100   if (abs(f3) <= eps) go to 120
    end do

120 s2 = delta * (co(1) * s1 + si(1) * c1)
    c2 = delta * (co(1) * c1 - si(1) * s1)
130 ress(j) = s2
    resc(j) = c2
    do i = 1, 4
      co(i) = co(i + 1)
      si(i) = si(i + 1)
    end do

  end do

  !  calculate the integrals ress(j) and resc(j),j=n-6,n-5,n-4 by setting
  !  up a divided difference table.
160 do j = 1, 3
    nmj = nm3 - j
    i = 5 - j
    call fpcsin(t(nm3), t(nmj), par, si(4), co(4), si(i - 1), co(i - 1), rs(j), rc(j))
    hs(i) = 0.
    hc(i) = 0.
    do jj = 1, j
      ipj = i + jj
      hc(ipj) = rc(jj)
      hs(ipj) = rs(jj)
    end do

    do jj = 1, 3
      if (i < jj) i = jj
      k = 5
      li = nmj
      do ll = i, 4
        lj = li + jj
        fac = t(lj) - t(li)
        hs(k) = (hs(k - 1) - hs(k)) / fac
        hc(k) = (hc(k - 1) - hc(k)) / fac
        k = k - 1
        li = li + 1
      end do

    end do

    ress(nmj) = hs(4) - hs(5)
    resc(nmj) = hc(4) - hc(5)
  end do

end subroutine fpbfou
