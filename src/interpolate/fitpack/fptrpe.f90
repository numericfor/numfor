subroutine fptrpe(m, mm, idim, n, nr, sp, p, b, z, a, aa, q, right)
!  subroutine fptrpe reduces the (m+n-7) x (n-7) cyclic bandmatrix a
!  to upper triangular form and applies the same givens transformations
!  to the (m) x (mm) x (idim) matrix z to obtain the (n-7) x (mm) x
!  (idim) matrix q.
!  ..
!  ..scalar arguments..
  real(8) :: p
  integer :: m, mm, idim, n
!  ..array arguments..
  real(8) :: sp(m, 4), b(n, 5), z(m * mm * idim), a(n, 5), aa(n, 4), &
    &q((n - 7) * mm * idim), right(mm * idim)
  integer :: nr(m)
!  ..local scalars..
  real(8) :: co, pinv, piv, si, one
  integer :: i, irot, it, ii, i2, i3, j, jj, l, mid, nmd, m2, m3,&
    &nrold, n4, number, n1, n7, n11, m1
!  ..local arrays..
  real(8) :: h(5), h1(5), h2(4)
!  ..subroutine references..
!    fpgivs,fprota
!  ..
  one = 1
  if (p > 0.) pinv = one / p
  n4 = n - 4
  n7 = n - 7
  n11 = n - 11
  mid = mm * idim
  m2 = m * mm
  m3 = n7 * mm
  m1 = m - 1
!  we determine the matrix (a) and then we reduce her to
!  upper triangular form (r) using givens rotations.
!  we apply the same transformations to the rows of matrix
!  z to obtain the (mm) x (n-7) matrix g.
!  we store matrix (r) into a and aa, g into q.
!  the n7 x n7 upper triangular matrix (r) has the form
!             | a1 '     |
!       (r) = |    ' a2  |
!             |  0 '     |
!  with (a2) a n7 x 4 matrix and (a1) a n11 x n11 upper
!  triangular matrix of bandwidth 5.
!  initialization.
  nmd = n7 * mid

  q(:nmd) = 0.

  a(:n4, 5) = 0.

  a(:n4, :4) = 0.
  aa(:n4, :4) = 0.

  jper = 0
  nrold = 0
  do 760 it = 1, m1
    number = nr(it)
120 if (nrold == number) go to 180
    if (p <= 0.) go to 740
!  fetch a new row of matrix (b).
    n1 = nrold + 1

    h(:5) = b(n1, :5) * pinv

!  find the appropriate row of q.

    right(:mid) = 0.

    go to 240
!  fetch a new row of matrix (sp)
180 h(5) = 0.
    h(:4) = sp(it, :4)

!  find the appropriate row of q.
    j = 0
    do ii = 1, idim
      l = (ii - 1) * m2 + (it - 1) * mm
      do jj = 1, mm
        j = j + 1
        l = l + 1
        right(j) = z(l)

      end do
    end do

!  test whether there are non-zero values in the new row of (a)
!  corresponding to the b-splines n(j,*),j=n7+1,...,n4.
240 if (nrold < n11) go to 640
    if (jper /= 0) go to 320
!  initialize the matrix (aa).
    jk = n11 + 1
    do 300 i = 1, 4
      ik = jk
      do 260 j = 1, 5
        if (ik <= 0) go to 280
        aa(ik, i) = a(ik, j)
        ik = ik - 1
260     continue
280     jk = jk + 1
300     continue
        jper = 1
!  if one of the non-zero elements of the new row corresponds to one of
!  the b-splines n(j;*),j=n7+1,...,n4,we take account of the periodicity
!  conditions for setting up this row of (a).
320     do 340 i = 1, 4
          h1(i) = 0.
          h2(i) = 0.
340       continue
          h1(5) = 0.
          j = nrold - n11
          do 420 i = 1, 5
            j = j + 1
            l0 = j
360         l1 = l0 - 4
            if (l1 <= 0) go to 400
            if (l1 <= n11) go to 380
            l0 = l1 - n11
            go to 360
380         h1(l1) = h(i)
            go to 420
400         h2(l0) = h2(l0) + h(i)
420         continue
!  rotate the new row of (a) into triangle.
            if (n11 <= 0) go to 560
!  rotations with the rows 1,2,...,n11 of (a).
            do 540 irot = 1, n11
              piv = h1(1)
              i2 = min0(n11 - irot, 4)
              if (piv == 0.) go to 500
!  calculate the parameters of the givens transformation.
              call fpgivs(piv, a(irot, 1), co, si)
!  apply that transformation to the columns of matrix q.
              j = 0
              do ii = 1, idim
                l = (ii - 1) * m3 + irot
                do jj = 1, mm
                  j = j + 1
                  call fprota(co, si, right(j), q(l))
                  l = l + n7
                end do
              end do

!  apply that transformation to the rows of (a) with respect to aa.
              do 460 i = 1, 4
                call fprota(co, si, h2(i), aa(irot, i))
460             continue
!  apply that transformation to the rows of (a) with respect to a.
                if (i2 == 0) go to 560
                do 480 i = 1, i2
                  i1 = i + 1
                  call fprota(co, si, h1(i1), a(irot, i1))
480               continue
500               do 520 i = 1, i2
                    h1(i) = h1(i + 1)
520                 continue
                    h1(i2 + 1) = 0.
540                 continue
!  rotations with the rows n11+1,...,n7 of a.
560                 do irot = 1, 4
                      ij = n11 + irot
                      if (ij <= 0) go to 620
                      piv = h2(irot)
                      if (piv == 0.) go to 620
!  calculate the parameters of the givens transformation.
                      call fpgivs(piv, aa(ij, irot), co, si)
!  apply that transformation to the columns of matrix q.
                      j = 0
                      do ii = 1, idim
                        l = (ii - 1) * m3 + ij
                        do jj = 1, mm
                          j = j + 1
                          call fprota(co, si, right(j), q(l))
                          l = l + n7
                        end do
                      end do

                      if (irot == 4) go to 620
!  apply that transformation to the rows of (a) with respect to aa.
                      j1 = irot + 1
                      do i = j1, 4
                        call fprota(co, si, h2(i), aa(ij, i))
                      end do

620                   continue
                    end do

                    go to 720
!  rotation into triangle of the new row of (a), in case the elements
!  corresponding to the b-splines n(j;*),j=n7+1,...,n4 are all zero.
640                 irot = nrold
                    do 700 i = 1, 5
                      irot = irot + 1
                      piv = h(i)
                      if (piv == 0.) go to 700
!  calculate the parameters of the givens transformation.
                      call fpgivs(piv, a(irot, 1), co, si)
!  apply that transformation to the columns of matrix g.
                      j = 0
                      do ii = 1, idim
                        l = (ii - 1) * m3 + irot
                        do jj = 1, mm
                          j = j + 1
                          call fprota(co, si, right(j), q(l))
                          l = l + n7
                        end do
                      end do

!  apply that transformation to the rows of (a).
                      if (i == 5) go to 700
                      i2 = 1
                      i3 = i + 1
                      do 680 j = i3, 5
                        i2 = i2 + 1
                        call fprota(co, si, h(j), a(irot, i2))
680                     continue
700                     continue
720                     if (nrold == number) go to 760
740                     nrold = nrold + 1
                        go to 120
760                     continue
                        return
                      end
