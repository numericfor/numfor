      subroutine fpinst(iopt, t, n, c, k, x, l, tt, nn, cc, nest)
!  given the b-spline representation (knots t(j),j=1,2,...,n, b-spline
!  coefficients c(j),j=1,2,...,n-k-1) of a spline of degree k, fpinst
!  calculates the b-spline representation (knots tt(j),j=1,2,...,nn,
!  b-spline coefficients cc(j),j=1,2,...,nn-k-1) of the same spline if
!  an additional knot is inserted at the point x situated in the inter-
!  val t(l)<=x<t(l+1). iopt denotes whether (iopt /= 0) or not (iopt=0)
!  the given spline is periodic. in case of a periodic spline at least
!  one of the following conditions must be fulfilled: l>2*k or l<n-2*k.
!
!  ..scalar arguments..
        integer :: k, n, l, nn, iopt, nest
        real(8) :: x
!  ..array arguments..
        real(8) :: t(nest), c(nest), tt(nest), cc(nest)
!  ..local scalars..
        real(8) :: fac, per, one
        integer :: i, i1, j, k1, m, mk, nk, nk1, nl, ll
!  ..
        one = 0.1e+01
        k1 = k + 1
        nk1 = n - k1
!  the new knots
        ll = l + 1
        i = n
        do 10 j = ll, n
          tt(i + 1) = t(i)
          i = i - 1
10        continue
          tt(ll) = x
          do 20 j = 1, l
            tt(j) = t(j)
20          continue
!  the new b-spline coefficients
            i = nk1
            do 30 j = l, nk1
              cc(i + 1) = c(i)
              i = i - 1
30            continue
              i = l
              do 40 j = 1, k
                m = i + k1
                fac = (x - tt(i)) / (tt(m) - tt(i))
                i1 = i - 1
                cc(i) = fac * c(i) + (one - fac) * c(i1)
                i = i1
40              continue
                do 50 j = 1, i
                  cc(j) = c(j)
50                continue
                  nn = n + 1
                  if (iopt == 0) return
!   incorporate the boundary conditions for a periodic spline.
                  nk = nn - k
                  nl = nk - k1
                  per = tt(nk) - tt(k1)
                  i = k1
                  j = nk
                  if (ll <= nl) go to 70
                  do 60 m = 1, k
                    mk = m + nl
                    cc(m) = cc(mk)
                    i = i - 1
                    j = j - 1
                    tt(i) = tt(j) - per
60                  continue
                    return
70                  if (ll > (k1 + k)) return
                    do 80 m = 1, k
                      mk = m + nl
                      cc(mk) = cc(m)
                      i = i + 1
                      j = j + 1
                      tt(j) = tt(i) + per
80                    continue
                      return
                    end
