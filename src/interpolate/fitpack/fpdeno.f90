      subroutine fpdeno(maxtr, up, left, right, nbind, merk)
!  subroutine fpdeno frees the nodes of all branches of a triply linked
!  tree with length < nbind by putting to zero their up field.
!  on exit the parameter merk points to the terminal node of the
!  most left branch of length nbind or takes the value 1 if there
!  is no such branch.
!  ..
!  ..scalar arguments..
        integer :: maxtr, nbind, merk
!  ..array arguments..
        integer :: up(maxtr), left(maxtr), right(maxtr)
!  ..local scalars ..
        integer :: i, j, k, l, niveau, point
!  ..
        i = 1
        niveau = 0
10      point = i
        i = left(point)
        if (i == 0) go to 20
        niveau = niveau + 1
        go to 10
20      if (niveau == nbind) go to 70
30      i = right(point)
        j = up(point)
        up(point) = 0
        k = left(j)
        if (point /= k) go to 50
        if (i /= 0) go to 40
        niveau = niveau - 1
        if (niveau == 0) go to 80
        point = j
        go to 30
40      left(j) = i
        go to 10
50      l = right(k)
        if (point == l) go to 60
        k = l
        go to 50
60      right(k) = i
        point = k
70      i = right(point)
        if (i /= 0) go to 10
        i = up(point)
        niveau = niveau - 1
        if (niveau == 0) go to 80
        point = i
        go to 70
80      k = 1
        l = left(k)
        if (up(l) == 0) return
90      merk = k
        k = left(k)
        if (k /= 0) go to 90
        return
      end
