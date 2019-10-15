subroutine fpintb(t, n, bint, nk1, x, y)
  implicit none
  !  subroutine fpintb calculates integrals of the normalized b-splines
  !  nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n.
  !  it makes use of the formulae of gaffney for the calculation of
  !  indefinite integrals of b-splines.
  !
  !  calling sequence:
  !     call fpintb(t,n,bint,nk1,x,y)
  !
  !  input parameters:
  !    t    : real array,length n, containing the position of the knots.
  !    n    : integer value, giving the number of knots.
  !    nk1  : integer value, giving the number of b-splines of degree k,
  !           defined on the set of knots ,i.e. nk1 = n-k-1.
  !    x,y  : real values, containing the end points of the integration
  !           interval.
  !  output parameter:
  !    bint : array,length nk1, containing the integrals of the b-splines.
  !  ..
  !  ..scalars arguments..
  integer :: n, nk1
  real(8) :: x, y
  !  ..array arguments..
  real(8) :: t(n), bint(nk1)
  !  ..local scalars..
  integer :: i, ia, ib, it, j, j1, k, k1, l, li, lj, lk
  real(8) :: a, ak, arg, b, f, one
  !  ..local arrays..
  real(8) :: aint(6), h(6), h1(6)
  logical :: minus
  !  initialization.
  one = 1._8
  k1 = n - nk1
  ak = k1
  k = k1 - 1

  bint(:nk1) = 0._8

  !  the integration limits are arranged in increasing order.
  if (x < y) then
    a = x; b = y
    minus = .False.
  else if (x > y) then
    a = y; b = x
    minus = .True.
  else
    return
  end if

  if (a < t(k1)) a = t(k1)
  if (b > t(nk1 + 1)) b = t(nk1 + 1)
  if (a > b) return
  !  using the expression of gaffney for the indefinite integral of a
  !  b-spline we find that
  !  bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
  !    where for t(l) <= x < t(l+1)
  !    res(j,x) = 0, j=1,2,...,l-k-1
  !             = 1, j=l+1,l+2,...,nk1
  !             = aint(j+k-l+1), j=l-k,l-k+1,...,l
  !               = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
  !                 i=0,1,...,k
  l = k1
  ! l0 = l + 1
  !  set arg = a.
  arg = a
  do it = 1, 2
    !  search for the knot interval t(l) <= arg < t(l+1).
    do while (arg >= t(l + 1) .and. l /= nk1)
      ! l = l0
      ! l0 = l + 1
      l = l + 1
    end do

    !  calculation of aint(j), j=1,2,...,k+1.
    !  initialization.
    aint(:k1) = 0._8
    aint(1) = (arg - t(l)) / (t(l + 1) - t(l))
    h1(1) = one
    do j = 1, k
      !  evaluation of the non-zero b-splines of degree j at arg,i.e.
      !    h(i+1) = nl-j+i,j(arg), i=0,1,...,j.
      h(1) = 0._8
      do i = 1, j
        li = l + i
        lj = li - j
        f = h1(i) / (t(li) - t(lj))
        h(i) = h(i) + f * (t(li) - arg)
        h(i + 1) = f * (arg - t(lj))
      end do

      !  updating of the integrals aint.
      j1 = j + 1
      do i = 1, j1
        li = l + i
        lj = li - j1
        aint(i) = aint(i) + h(i) * (arg - t(lj)) / (t(li) - t(lj))
        h1(i) = h(i)
      end do

    end do

    if (it == 2) exit
    !  updating of the integrals bint
    lk = l - k
    ia = lk
    do i = 1, k1
      bint(lk) = -aint(i)
      lk = lk + 1
    end do
    arg = b
  end do

  !  updating of the integrals bint.
  lk = l - k
  ib = lk - 1
  do i = 1, k1
    bint(lk) = bint(lk) + aint(i)
    lk = lk + 1
  end do

  if (ib >= ia) bint(ia:ib) = bint(ia:ib) + one

  !  the scaling factors are taken into account.
  f = one / ak
  do i = 1, nk1
    j = i + k1
    bint(i) = bint(i) * (t(j) - t(i)) * f
  end do

  !  the order of the integration limits is taken into account.
  IF (minus) bint(:nk1) = -bint(:nk1)
end subroutine fpintb

