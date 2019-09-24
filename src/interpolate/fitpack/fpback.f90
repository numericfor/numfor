subroutine fpback(a, z, n, k, c, nest)
  !  subroutine fpback calculates the solution of the system of
  !  equations a*c = z with a a n x n upper triangular matrix
  !  of bandwidth k.
  !  ..
  !  ..scalar arguments..
  integer :: n, k, nest
  !  ..array arguments..
  real(8) :: a(nest, k), z(n), c(n)
  !  ..local scalars..
  real(8) :: store
  integer :: i, i1, j, k1, l, m
  !  ..
  k1 = k - 1
  c(n) = z(n) / a(n, 1)
  i = n - 1
  if (i == 0) return
  do j = 2, n
    store = z(i)
    i1 = k1
    if (j <= k1) i1 = j - 1
    m = i
    do l = 1, i1
      m = m + 1
      store = store - c(m) * a(i, l + 1)
    end do

    c(i) = store / a(i, 1)
    i = i - 1
  end do
end subroutine fpback
