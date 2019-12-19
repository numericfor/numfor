subroutine fpbacp(a, b, z, n, k, c, k1, nest)
  !  subroutine fpbacp calculates the solution of the system of equations
  !  g * c = z  with g  a n x n upper triangular matrix of the form
  !            ! a '   !
  !        g = !   ' b !
  !            ! 0 '   !
  !  with b a n x k matrix and a a (n-k) x (n-k) upper triangular
  !  matrix of bandwidth k1.
  !  ..
  !  ..scalar arguments..
  integer :: n, k, k1, nest
  !  ..array arguments..
  real(8) :: a(nest, k1), b(nest, k), z(n), c(n)
  !  ..local scalars..
  integer :: i, i1, j, l, l0, l1, n2
  real(8) :: store
  !  ..
  n2 = n - k
  l = n
  do i = 1, k
    store = z(l)
    j = k + 2 - i
    if (i > 1) then
      l0 = l
      do l1 = j, k
        l0 = l0 + 1
        store = store - c(l0) * b(l, l1)
      end do
    end if

    c(l) = store / b(l, j - 1)
    l = l - 1
    if (l == 0) return
  end do

  do i = 1, n2
    store = z(i)
    l = n2
    do j = 1, k
      l = l + 1
      store = store - c(l) * b(i, j)
    end do

    c(i) = store
  end do

  i = n2
  c(i) = c(i) / a(i, 1)

  if (i == 1) return

  do j = 2, n2
    i = i - 1
    store = c(i)
    i1 = k
    if (j <= k) i1 = j - 1
    l = i
    do l0 = 1, i1
      l = l + 1
      store = store - c(l) * a(i, l0 + 1)
    end do
    c(i) = store / a(i, 1)
  end do

end subroutine fpbacp
