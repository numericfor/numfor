      subroutine fpchep(x, m, t, n, k, ier)
!  subroutine fpchep verifies the number and the position of the knots
!  t(j),j=1,2,...,n of a periodic spline of degree k, in relation to
!  the number and the position of the data points x(i),i=1,2,...,m.
!  if all of the following conditions are fulfilled, ier is set
!  to zero. if one of the conditions is violated ier is set to ten.
!      1) k+1 <= n-k-1 <= m+k-1
!      2) t(1) <= t(2) <= ... <= t(k+1)
!         t(n-k) <= t(n-k+1) <= ... <= t(n)
!      3) t(k+1) < t(k+2) < ... < t(n-k)
!      4) t(k+1) <= x(i) <= t(n-k)
!      5) the conditions specified by schoenberg and whitney must hold
!         for at least one subset of data points, i.e. there must be a
!         subset of data points y(j) such that
!             t(j) < y(j) < t(j+k+1), j=k+1,...,n-k-1
!  ..
!  ..scalar arguments..
        integer :: m, n, k, ier
!  ..array arguments..
        real(8) :: x(m), t(n)
!  ..local scalars..
        integer :: i, i1, i2, j, j1, k1, k2, l, l1, l2, mm, m1, nk1, nk2
        real(8) :: per, tj, tl, xi
!  ..
        k1 = k + 1
        k2 = k1 + 1
        nk1 = n - k1
        nk2 = nk1 + 1
        m1 = m - 1
        ier = 10
!  check condition no 1
        if (nk1 < k1 .or. n > m + 2 * k) go to 130
!  check condition no 2
        j = n
        do 20 i = 1, k
          if (t(i) > t(i + 1)) go to 130
          if (t(j) < t(j - 1)) go to 130
          j = j - 1
20        continue
!  check condition no 3
          do 30 i = k2, nk2
            if (t(i) <= t(i - 1)) go to 130
30          continue
!  check condition no 4
            if (x(1) < t(k1) .or. x(m) > t(nk2)) go to 130
!  check condition no 5
            l1 = k1
            l2 = 1
            do 50 l = 1, m
              xi = x(l)
40            if (xi < t(l1 + 1) .or. l == nk1) go to 50
              l1 = l1 + 1
              l2 = l2 + 1
              if (l2 > k1) go to 60
              go to 40
50            continue
              l = m
60            per = t(nk2) - t(k1)
              do 120 i1 = 2, l
                i = i1 - 1
                mm = i + m1
                do 110 j = k1, nk1
                  tj = t(j)
                  j1 = j + k1
                  tl = t(j1)
70                i = i + 1
                  if (i > mm) go to 120
                  i2 = i - m1
                  if (i2 <= 0) go to 80
                  go to 90
80                xi = x(i)
                  go to 100
90                xi = x(i2) + per
100               if (xi <= tj) go to 70
                  if (xi >= tl) go to 120
110               continue
                  ier = 0
                  go to 130
120               continue
130               return
                end
