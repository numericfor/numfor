subroutine fpfrno(maxtr, up, left, right, info, point, merk, n1, count, ier)
  !  subroutine fpfrno collects the free nodes (up field zero) of the
  !  triply linked tree the information of which is kept in the arrays
  !  up,left,right and info. the maximal length of the branches of the
  !  tree is given by n1. if no free nodes are found, the error flag
  !  ier is set to 1.
  !  ..
  !  ..scalar arguments..
  integer :: maxtr, point, merk, n1, count, ier
  !  ..array arguments..
  integer :: up(maxtr), left(maxtr), right(maxtr), info(maxtr)
  !  ..local scalars
  integer :: i, j, k, l, n, niveau
  !  ..
  ier = 1
  if (n1 == 2) return
  niveau = 1
  count = 2
10 j = 0
  i = 1
20 do while (j < niveau)
    k = 0
    l = left(i)
    if (l == 0) go to 110
    i = l
    j = j + 1
  end do

30 if (i < count) go to 110
  if (i == count) go to 100
  if (up(count) == 0) go to 50
  count = count + 1
  go to 30
50 up(count) = up(i)
  left(count) = left(i)
  right(count) = right(i)
  info(count) = info(i)
  if (merk == i) merk = count
  if (point == i) point = count
  if (k == 0) then
    n = up(i)
    left(n) = count
  else
    right(k) = count
  end if
  l = left(i)
  do while (l /= 0)
    up(l) = count
    l = right(l)
  end do
  up(i) = 0
  i = count
100 count = count + 1
110 l = right(i)
  k = i
  if (l == 0) go to 120
  i = l
  go to 20
120 l = up(i)
  j = j - 1
  if (j == 0) go to 130
  i = l
  go to 110
130 niveau = niveau + 1
  if (niveau <= n1) go to 10
  if (count > maxtr) return
  ier = 0
end subroutine fpfrno
