subroutine fpdisc(t, n, k2, b, nest)
  !  subroutine fpdisc calculates the discontinuity jumps of the kth
  !  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
  !  ..scalar arguments..
  integer :: n, k2, nest
  !  ..array arguments..
  real(8) :: t(n), b(nest, k2)
  !  ..local scalars..
  real(8) :: an, fac, prod
  integer :: i, ik, j, jk, k, k1, l, lj, lk, lmk, lp, nk1, nrint
  !  ..local array..
  real(8) :: h(12)
  !  ..
  k1 = k2 - 1
  k = k1 - 1
  nk1 = n - k1
  nrint = nk1 - k
  an = nrint
  fac = an / (t(nk1 + 1) - t(k1))
  do l = k2, nk1
    lmk = l - k1
    do j = 1, k1
      ik = j + k1
      lj = l + j
      lk = lj - k2
      h(j) = t(l) - t(lk)
      h(ik) = t(l) - t(lj)
    end do

    lp = lmk
    do j = 1, k2
      jk = j
      prod = h(j)
      do i = 1, k
        jk = jk + 1
        prod = prod * h(jk) * fac
      end do

      lk = lp + k1
      b(lmk, j) = (t(lk) - t(lp)) / prod
      lp = lp + 1
    end do

  end do

  return
end subroutine fpdisc
