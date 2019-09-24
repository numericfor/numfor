subroutine fpchec(x, m, t, n, k, ier)
  !  subroutine fpchec verifies the number and the position of the knots
  !  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
  !  and the position of the data points x(i),i=1,2,...,m. if all of the
  !  following conditions are fulfilled, the error parameter ier is set
  !  to zero. if one of the conditions is violated ier is set to ten.
  !      1) k+1 <= n-k-1 <= m
  !      2) t(1) <= t(2) <= ... <= t(k+1)
  !         t(n-k) <= t(n-k+1) <= ... <= t(n)
  !      3) t(k+1) < t(k+2) < ... < t(n-k)
  !      4) t(k+1) <= x(i) <= t(n-k)
  !      5) the conditions specified by schoenberg and whitney must hold
  !         for at least one subset of data points, i.e. there must be a
  !         subset of data points y(j) such that
  !             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
  !  ..
  !  ..scalar arguments..
  integer :: m, n, k, ier
  !  ..array arguments..
  real(8) :: x(m), t(n)
  !  ..local scalars..
  integer :: i, j, k1, k2, l, nk1, nk2, nk3
  real(8) :: tj, tl
  !  ..
  k1 = k + 1
  k2 = k1 + 1
  nk1 = n - k1
  nk2 = nk1 + 1
  ier = 10
  !  check condition no 1
  if (nk1 < k1 .or. nk1 > m) return
  !  check condition no 2
  j = n
  do i = 1, k
    if (t(i) > t(i + 1)) return
    if (t(j) < t(j - 1)) return
    j = j - 1
  end do

  !  check condition no 3
  do i = k2, nk2
    if (t(i) <= t(i - 1)) return
  end do

  !  check condition no 4
  if (x(1) < t(k1) .or. x(m) > t(nk2)) return
  !  check condition no 5
  if (x(1) >= t(k2) .or. x(m) <= t(nk1)) return
  i = 1
  l = k2
  nk3 = nk1 - 1
  if (nk3 < 2) then
    ier = 0
    return
  end if

  do j = 2, nk3
    tj = t(j)
    l = l + 1
    tl = t(l)
40  i = i + 1
    if (i >= m) return
    if (x(i) <= tj) go to 40
    if (x(i) >= tl) return
  end do

end subroutine fpchec
