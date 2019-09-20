      subroutine fpader(t, n, c, k1, x, l, d)
!  subroutine fpader calculates the derivatives
!             (j-1)
!     d(j) = s     (x) , j=1,2,...,k1
!  of a spline of order k1 at the point t(l)<=x<t(l+1), using the
!  stable recurrence scheme of de boor
!  ..
!  ..scalar arguments..
        real(8) :: x
        integer :: n, k1, l
!  ..array arguments..
        real(8) :: t(n), c(n), d(k1)
!  ..local scalars..
        integer :: i, ik, j, jj, j1, j2, ki, kj, li, lj, lk
        real(8) :: ak, fac, one
!  ..local array..
        real(8) :: h(20)
!  ..
        one = 1._8
        lk = l - k1
        do i = 1, k1
          ik = i + lk
          h(i) = c(ik)
        end do

        kj = k1
        fac = one
        do j = 1, k1
          ki = kj
          j1 = j + 1
          if (j == 1) go to 300
          i = k1
          do jj = j, k1
            li = i + lk
            lj = li + kj
            h(i) = (h(i) - h(i - 1)) / (t(lj) - t(li))
            i = i - 1
          end do

300       do i = j, k1
            d(i) = h(i)
          end do

          if (j == k1) go to 600
          do jj = j1, k1
            ki = ki - 1
            i = k1
            do j2 = jj, k1
              li = i + lk
              lj = li + ki
              d(i) = ((x - t(li)) * d(i) + (t(lj) - x) * d(i - 1)) / (t(lj) - t(li))
              i = i - 1
            end do
          end do

600       d(j) = d(k1) * fac
          ak = k1 - j
          fac = fac * ak
          kj = kj - 1
        end do

        return
      end
