subroutine fprank(a, f, n, m, na, tol, c, sq, rank, aa, ff, h)
  !  subroutine fprank finds the minimum norm solution of a least-
  !  squares problem in case of rank deficiency.
  !
  !  input parameters:
  !    a : array, which contains the non-zero elements of the observation
  !        matrix after triangularization by givens transformations.
  !    f : array, which contains the transformed right hand side.
  !    n : integer,which contains the dimension of a.
  !    m : integer, which denotes the bandwidth of a.
  !  tol : real value, giving a threshold to determine the rank of a.
  !
  !  output parameters:
  !    c : array, which contains the minimum norm solution.
  !   sq : real value, giving the contribution of reducing the rank
  !        to the sum of squared residuals.
  ! rank : integer, which contains the rank of matrix a.
  !
  !  ..scalar arguments..
  integer :: n, m, na, rank
  real(8) :: tol, sq
  !  ..array arguments..
  real(8) :: a(na, m), f(n), c(n), aa(n, m), ff(n), h(m)
  !  ..local scalars..
  integer :: i, ii, ij, i1, i2, j, jj, j1, j2, j3, k, kk, m1, nl
  real(8) :: cos, fac, piv, sin, yi
  double precision store, stor1, stor2, stor3
  !  ..function references..
  integer :: min0
  !  ..subroutine references..
  !    fpgivs,fprota
  !  ..
  m1 = m - 1
  !  the rank deficiency nl is considered to be the number of sufficient
  !  small diagonal elements of a.
  nl = 0
  sq = 0.
  do i = 1, n
    if (a(i, 1) > tol) cycle
    !  if a sufficient small diagonal element is found, we put it to
    !  zero. the remainder of the row corresponding to that zero diagonal
    !  element is then rotated into triangle by givens rotations .
    !  the rank deficiency is increased by one.
    nl = nl + 1
    if (i == n) cycle
    yi = f(i)

    h(:m1) = a(i, 2:m1 + 1)

    h(m) = 0.
    i1 = i + 1
    do ii = i1, n
      i2 = min0(n - ii, m1)
      piv = h(1)
      if (piv == 0.) go to 30
      call fpgivs(piv, a(ii, 1), cos, sin)
      call fprota(cos, sin, yi, f(ii))
      if (i2 == 0) go to 70
      do j = 1, i2
        j1 = j + 1
        call fprota(cos, sin, h(j1), a(ii, j1))
        h(j) = h(j1)
      end do

      go to 50
30    if (i2 == 0) go to 70
      do j = 1, i2
        h(j) = h(j + 1)
      end do

50    h(i2 + 1) = 0.
    end do

    !  add to the sum of squared residuals the contribution of deleting
    !  the row with small diagonal element.
70  sq = sq + yi**2
  end do

  !  rank denotes the rank of a.
  rank = n - nl
  !  let b denote the (rank*n) upper trapezoidal matrix which can be
  !  obtained from the (n*n) upper triangular matrix a by deleting
  !  the rows and interchanging the columns corresponding to a zero
  !  diagonal element. if this matrix is factorized using givens
  !  transformations as  b = (r) (u)  where
  !    r is a (rank*rank) upper triangular matrix,
  !    u is a (rank*n) orthonormal matrix
  !  then the minimal least-squares solution c is given by c = b' v,
  !  where v is the solution of the system  (r) (r)' v = g  and
  !  g denotes the vector obtained from the old right hand side f, by
  !  removing the elements corresponding to a zero diagonal element of a.
  !  initialization.

  aa(:rank, :m) = 0._8

  !  form in aa the upper triangular matrix obtained from a by
  !  removing rows and columns with zero diagonal elements. form in ff
  !  the new right hand side by removing the elements of the old right
  !  hand side corresponding to a deleted row.
  ii = 0
  do i = 1, n
    if (a(i, 1) <= tol) go to 120
    ii = ii + 1
    ff(ii) = f(i)
    aa(ii, 1) = a(i, 1)
    jj = ii
    kk = 1
    j = i
    j1 = min0(j - 1, m1)
    if (j1 == 0) go to 120
    do k = 1, j1
      j = j - 1
      if (a(j, 1) <= tol) go to 110
      kk = kk + 1
      jj = jj - 1
      aa(jj, kk) = a(j, k + 1)
110   continue
    end do

120 continue
  end do

  !  form successively in h the columns of a with a zero diagonal element.
  ii = 0
  do i = 1, n
    ii = ii + 1
    if (a(i, 1) > tol) cycle
    ii = ii - 1
    if (ii == 0) cycle
    jj = 1
    j = i
    j1 = min0(j - 1, m1)
    do k = 1, j1
      j = j - 1
      if (a(j, 1) > tol) then
        h(jj) = a(j, k + 1)
        jj = jj + 1
      end if
    end do

    h(jj:m) = 0._8

    !  rotate this column into aa by givens transformations.
    jj = ii
    do i1 = 1, ii
      j1 = min0(jj - 1, m1)
      piv = h(1)
      if (piv /= 0.) go to 160
      if (j1 == 0) cycle
      do j2 = 1, j1
        j3 = j2 + 1
        h(j2) = h(j3)
      end do

      go to 180
160   call fpgivs(piv, aa(jj, 1), cos, sin)
      if (j1 == 0) cycle
      kk = jj
      do j2 = 1, j1
        j3 = j2 + 1
        kk = kk - 1
        call fprota(cos, sin, h(j3), aa(kk, j3))
        h(j2) = h(j3)
      end do

180   jj = jj - 1
      h(j3) = 0._8
    end do
  end do

  !  solve the system (aa) (f1) = ff
  ff(rank) = ff(rank) / aa(rank, 1)
  i = rank - 1
  if (i == 0) go to 230
  do j = 2, rank
    store = ff(i)
    i1 = min0(j - 1, m1)
    k = i
    do ii = 1, i1
      k = k + 1
      stor1 = ff(k)
      stor2 = aa(i, ii + 1)
      store = store - stor1 * stor2
    end do

    stor1 = aa(i, 1)
    ff(i) = store / stor1
    i = i - 1
  end do

  !  solve the system  (aa)' (f2) = f1
230 ff(1) = ff(1) / aa(1, 1)
  if (rank == 1) go to 260
  do j = 2, rank
    store = ff(j)
    i1 = min0(j - 1, m1)
    k = j
    do ii = 1, i1
      k = k - 1
      stor1 = ff(k)
      stor2 = aa(k, ii + 1)
      store = store - stor1 * stor2
    end do

    stor1 = aa(j, 1)
    ff(j) = store / stor1
  end do

  !  premultiply f2 by the transpoze of a.
260 k = 0
  do i = 1, n
    store = 0.
    if (a(i, 1) > tol) k = k + 1
    j1 = min0(i, m)
    kk = k
    ij = i + 1
    do j = 1, j1
      ij = ij - 1
      if (a(ij, 1) <= tol) cycle
      stor1 = a(ij, j)
      stor2 = ff(kk)
      store = store + stor1 * stor2
      kk = kk - 1
    end do
    c(i) = store
  end do

  !  add to the sum of squared residuals the contribution of putting
  !  to zero the small diagonal elements of matrix (a).
  stor3 = 0.
  do i = 1, n
    if (a(i, 1) <= tol) then
      store = f(i)
      i1 = min0(n - i, m1)
      if (i1 /= 0) then
        do j = 1, i1
          ij = i + j
          stor1 = c(ij)
          stor2 = a(i, j + 1)
          store = store - stor1 * stor2
        end do
      end if

      fac = a(i, 1) * c(i)
      stor1 = a(i, 1)
      stor2 = c(i)
      stor1 = stor1 * stor2
      stor3 = stor3 + stor1 * (stor1 - store - store)
    end if
  end do

  fac = stor3
  sq = sq + fac
end subroutine fprank

