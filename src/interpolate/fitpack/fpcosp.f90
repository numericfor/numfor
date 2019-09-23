subroutine fpcosp(m, x, y, w, n, t, e, maxtr, maxbin, c, sq, sx, bind, nm, mb, a,&
  &b, const, z, zz, u, q, info, up, left, right, jbind, ibind, ier)
  !  ..
  !  ..scalar arguments..
  real(8) :: sq
  integer :: m, n, maxtr, maxbin, nm, mb, ier
  !  ..array arguments..
  real(8) :: x(m), y(m), w(m), t(n), e(n), c(n), sx(m), a(n, 4), b(nm, maxbin),&
    &const(n), z(n), zz(n), u(maxbin), q(m, 4)
  integer :: info(maxtr), up(maxtr), left(maxtr), right(maxtr), jbind(mb), ibind(mb)
  logical bind(n)
  !  ..local scalars..
  integer :: count, i, i1, j, j1, j2, j3, k, kdim, k1, k2, k3, k4, k5, k6,&
    &l, lp1, l1, l2, l3, merk, nbind, number, n1, n4, n6
  real(8) :: f, wi, xi
  !  ..local array..
  real(8) :: h(4)
  !  ..subroutine references..
  !    fpbspl,fpadno,fpdeno,fpfrno,fpseno
  !  ..
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  if we use the b-spline representation of s(x) our approximation     c
  !  problem results in a quadratic programming problem:                 c
  !    find the b-spline coefficients c(j),j=1,2,...n-4 such that        c
  !        (1) sumi((wi*(yi-sumj(cj*nj(xi))))**2),i=1,2,...m is minimal  c
  !        (2) sumj(cj*n''j(t(l+3)))*e(l) <= 0, l=1,2,...n-6.            c
  !  to solve this problem we use the theil-van de panne procedure.      c
  !  if the inequality constraints (2) are numbered from 1 to n-6,       c
  !  this algorithm finds a subset of constraints ibind(1)..ibind(nbind) c
  !  such that the solution of the minimization problem (1) with these   c
  !  constraints in equality form, satisfies all constraints. such a     c
  !  feasible solution is optimal if the lagrange parameters associated  c
  !  with that problem with equality constraints, are all positive.      c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  determine n6, the number of inequality constraints.
  n6 = n - 6
  !  fix the parameters which determine these constraints.
  do i = 1, n6
    const(i) = e(i) * (t(i + 4) - t(i + 1)) / (t(i + 5) - t(i + 2))
  end do

  !  initialize the triply linked tree which is used to find the subset
  !  of constraints ibind(1),...ibind(nbind).
  count = 1
  info(1) = 0
  left(1) = 0
  right(1) = 0
  up(1) = 1
  merk = 1
  !  set up the normal equations  n'nc=n'y  where n denotes the m x (n-4)
  !  observation matrix with elements ni,j = wi*nj(xi)  and y is the
  !  column vector with elements yi*wi.
  !  from the properties of the b-splines nj(x),j=1,2,...n-4, it follows
  !  that  n'n  is a (n-4) x (n-4)  positive definite bandmatrix of
  !  bandwidth 7. the matrices n'n and n'y are built up in a and z.
  n4 = n - 4
  !  initialization

  z(:n4) = 0._8
  a(:n4, :4) = 0._8

  l = 4
  lp1 = l + 1
  do i = 1, m
    !  fetch the current row of the observation matrix.
    xi = x(i)
    wi = w(i)**2
    !  search for knot interval  t(l) <= xi < t(l+1)
30  if (xi < t(lp1) .or. l == n4) go to 40
    l = lp1
    lp1 = l + 1
    go to 30
    !  evaluate the four non-zero cubic b-splines nj(xi),j=l-3,...l.
40  call fpbspl(t, n, 3, xi, l, h)
    !  store in q these values h(1),h(2),...h(4).
    q(i, 1:4) = h(1:4)
    !  add the contribution of the current row of the observation matrix
    !  n to the normal equations.
    l3 = l - 3
    k1 = 0
    do j1 = l3, l
      k1 = k1 + 1
      f = h(k1)
      z(j1) = z(j1) + f * wi * y(i)
      k2 = k1
      j2 = 4
      do j3 = j1, l
        a(j3, j2) = a(j3, j2) + f * wi * h(k2)
        k2 = k2 + 1
        j2 = j2 - 1
      end do
    end do
  end do

  !  since n'n is a symmetric matrix it can be factorized as
  !        (3)  n'n = (r1)'(d1)(r1)
  !  with d1 a diagonal matrix and r1 an (n-4) x (n-4)  unit upper
  !  triangular matrix of bandwidth 4. the matrices r1 and d1 are built
  !  up in a. at the same time we solve the systems of equations
  !        (4)  (r1)'(z2) = n'y
  !        (5)  (d1) (z1) = (z2)
  !  the vectors z2 and z1 are kept in zz and z.
  do i = 1, n4
    k1 = 1
    if (i < 4) k1 = 5 - i
    k2 = i - 4 + k1
    k3 = k2
    do j = k1, 4
      k4 = j - 1
      k5 = 4 - j + k1
      f = a(i, j)
      if (k1 <= k4) then
        k6 = k2
        do k = k1, k4
          f = f - a(i, k) * a(k3, k5) * a(k6, 4)
          k5 = k5 + 1
          k6 = k6 + 1
        end do
      end if

      if (j == 4) exit
      a(i, j) = f / a(k3, 4)
      k3 = k3 + 1
    end do

    a(i, 4) = f
    f = z(i)
    if (i /= 1) then
      k4 = i
      do j = k1, 3
        k = k1 + 3 - j
        k4 = k4 - 1
        f = f - a(i, k) * z(k4) * a(k4, 4)
      end do
    end if

    z(i) = f / a(i, 4)
    zz(i) = f
  end do

  !  start computing the least-squares cubic spline without taking account
  !  of any constraint.
  nbind = 0
  n1 = 1
  ibind(1) = 0
  !  main loop for the least-squares problems with different subsets of
  !  the constraints (2) in equality form. the resulting b-spline coeff.
  !  c and lagrange parameters u are the solution of the system
  !            ! n'n  b' ! ! c !   ! n'y !
  !        (6) !         ! !   ! = !     !
  !            !  b   0  ! ! u !   !  0  !
  !  z1 is stored into array c.
150 c(:n4) = z(:n4)

  !  if there are no equality constraints, compute the coeff. c directly.
  if (nbind == 0) go to 370
  !  initialization
  kdim = n4 + nbind
  b(:kdim, :nbind) = 0._8

  !  matrix b is built up,expressing that the constraints nrs ibind(1),...
  !  ibind(nbind) must be satisfied in equality form.
  do i = 1, nbind
    l = ibind(i)
    b(l, i) = e(l)
    b(l + 1, i) = -(e(l) + const(l))
    b(l + 2, i) = const(l)
  end do

  !  find the matrix (b1) as the solution of the system of equations
  !        (7)  (r1)'(d1)(b1) = b'
  !  (b1) is built up in the upper part of the array b(rows 1,...n-4).
  do k1 = 1, nbind
    l = ibind(k1)
    do i = l, n4
      f = b(i, k1)
      if (i == 1) then
        b(i, k1) = f / a(i, 4)
      else
        k2 = 3
        if (i < 4) k2 = i - 1
        do k3 = 1, k2
          l1 = i - k3
          l2 = 4 - k3
          f = f - b(l1, k1) * a(i, l2) * a(l1, 4)
        end do
      end if
    end do
  end do

  !  factorization of the symmetric matrix  -(b1)'(d1)(b1)
  !        (8) ::  -(b1)'(d1)(b1) = (r2)'(d2)(r2)
  !  with (d2) a diagonal matrix and (r2) an nbind x nbind unit upper
  !  triangular matrix. the matrices r2 and d2 are built up in the lower
  !  part of the array b (rows n-3,n-2,...n-4+nbind).
  do i = 1, nbind
    i1 = i - 1
    do j = i, nbind
      f = 0.
      do k = 1, n4
        f = f + b(k, i) * b(k, j) * a(k, 4)
      end do

      k1 = n4 + 1
      if (i1 /= 0) then
        do k = 1, i1
          f = f + b(k1, i) * b(k1, j) * b(k1, k)
          k1 = k1 + 1
        end do
      end if

      b(k1, j) = -f
      if (j /= i) b(k1, j) = b(k1, j) / b(k1, i)
    end do
  end do

  !  according to (3),(7) and (8) :: the system of equations (6) becomes
  !         ! (r1)'    0  ! ! (d1)    0  ! ! (r1)  (b1) ! ! c !   ! n'y !
  !    (9)  !             ! !            ! !            ! !   ! = !     !
  !         ! (b1)'  (r2)'! !   0   (d2) ! !   0   (r2) ! ! u !   !  0  !
  !  backward substitution to obtain the b-spline coefficients c(j),j=1,..
  !  n-4 and the lagrange parameters u(j),j=1,2,...nbind.
  !  first step of the backward substitution: solve the system
  !             ! (r1)'(d1)      0     ! ! (c1) !   ! n'y !
  !        (10) !                      ! !      ! = !     !
  !             ! (b1)'(d1)  (r2)'(d2) ! ! (u1) !   !  0  !
  !  from (4) and (5) we know that this is equivalent to
  !        (11)  (c1) = (z1)
  !        (12)  (r2)'(d2)(u1) = -(b1)'(z2)
  do i = 1, nbind
    f = 0.
    do j = 1, n4
      f = f + b(j, i) * zz(j)
    end do

    i1 = i - 1
    k1 = n4 + 1
    if (i1 /= 0) then
      do j = 1, i1
        f = f + u(j) * b(k1, i) * b(k1, j)
        k1 = k1 + 1
      end do
    end if

    u(i) = -f / b(k1, i)
  end do

  !  second step of the backward substitution: solve the system
  !             ! (r1)  (b1) ! ! c !   ! c1 !
  !        (13) !            ! !   ! = !    !
  !             !   0   (r2) ! ! u !   ! u1 !
  k1 = nbind
  k2 = kdim
  !  find the lagrange parameters u.
  do i = 1, nbind
    f = u(k1)
    if (i /= 1) then
      k3 = k1 + 1
      do j = k3, nbind
        f = f - u(j) * b(k2, j)
      end do
    end if

    u(k1) = f
    k1 = k1 - 1
    k2 = k2 - 1
  end do

  !  find the b-spline coefficients c.
  do i = 1, n4
    f = c(i)
    do j = 1, nbind
      f = f - u(j) * b(i, j)
    end do

    c(i) = f
  end do

370 k1 = n4
  do i = 2, n4
    k1 = k1 - 1
    f = c(k1)
    k2 = 1
    if (i < 5) k2 = 5 - i
    k3 = k1
    l = 3
    do j = k2, 3
      k3 = k3 + 1
      f = f - a(k3, l) * c(k3)
      l = l - 1
    end do

    c(k1) = f
  end do

  !  test whether the solution of the least-squares problem with the
  !  constraints ibind(1),...ibind(nbind) in equality form, satisfies
  !  all of the constraints (2).
  k = 1
  !  number counts the number of violated inequality constraints.
  number = 0
  do j = 1, n6
    l = ibind(k)
    k = k + 1
    if (j == l) cycle
    k = k - 1
    !  test whether constraint j is satisfied
    f = e(j) * (c(j) - c(j + 1)) + const(j) * (c(j + 2) - c(j + 1))
    if (f <= 0.) cycle
    !  if constraint j is not satisfied, add a branch of length nbind+1
    !  to the tree. the nodes of this branch contain in their information
    !  field the number of the constraints ibind(1),...ibind(nbind) and j,
    !  arranged in increasing order.
    number = number + 1
    k1 = k - 1
    if (k1 /= 0) jbind(:k1) = ibind(:k1)
    jbind(k) = j
    if (l /= 0) jbind(k + 1:nbind + 1) = ibind(k:nbind)

    call fpadno(maxtr, up, left, right, info, count, merk, jbind, n1, ier)
    !  test whether the storage space which is required for the tree,exceeds
    !  the available storage space.
    if (ier /= 0) go to 560
  end do

  !  test whether the solution of the least-squares problem with equality
  !  constraints is a feasible solution.
  if (number == 0) go to 470
  !  test whether there are still cases with nbind constraints in
  !  equality form to be considered.
450 if (merk > 1) go to 460
  nbind = n1
  !  test whether the number of knots where s''(x)=0 exceeds maxbin.
  if (nbind > maxbin) go to 550
  n1 = n1 + 1
  ibind(n1) = 0
  !  search which cases with nbind constraints in equality form
  !  are going to be considered.
  call fpdeno(maxtr, up, left, right, nbind, merk)
  !  test whether the quadratic programming problem has a solution.
  if (merk == 1) go to 570
  !  find a new case with nbind constraints in equality form.
460 call fpseno(maxtr, up, left, right, info, merk, ibind, nbind)
  go to 150
  !  test whether the feasible solution is optimal.
470 ier = 0
  bind(:n6) = .false.

  if (nbind == 0) go to 500
  do i = 1, nbind
    if (u(i) <= 0.) go to 450
    j = ibind(i)
    bind(j) = .true.
  end do

  !  evaluate s(x) at the data points x(i) and calculate the weighted
  !  sum of squared residual right hand sides sq.
500 sq = 0.
  l = 4
  lp1 = 5
  do i = 1, m
510 if (x(i) < t(lp1) .or. l == n4) go to 520
    l = lp1
    lp1 = l + 1
    go to 510
520 sx(i) = c(l - 3) * q(i, 1) + c(l - 2) * q(i, 2) + c(l - 1) * q(i, 3) + c(l) * q(i, 4)
    sq = sq + (w(i) * (y(i) - sx(i)))**2
  end do

  go to 600
  !  error codes and messages.
550 ier = 1
  go to 600
560 ier = 2
  go to 600
570 ier = 3
600 return
end subroutine fpcosp
