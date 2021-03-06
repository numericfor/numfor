subroutine fpseno(maxtr, up, left, right, info, merk, ibind, nbind)
  !  subroutine fpseno fetches a branch of a triply linked tree the
  !  information of which is kept in the arrays up,left,right and info.
  !  the branch has a specified length nbind and is determined by the
  !  parameter merk which points to its terminal node. the information
  !  field of the nodes of this branch is stored in the array ibind. on
  !  exit merk points to a new branch of length nbind or takes the value
  !  1 if no such branch was found.
  !  ..
  !  ..scalar arguments..
  integer :: maxtr, merk, nbind
  !  ..array arguments..
  integer :: up(maxtr), left(maxtr), right(maxtr), info(maxtr),&
    &ibind(nbind)
  !  ..scalar arguments..
  integer :: i, j, k
  !  ..
  k = merk
  j = nbind
  do i = 1, nbind
    ibind(j) = info(k)
    k = up(k)
    j = j - 1
  end do

20 k = right(merk)
  if (k /= 0) go to 30
  merk = up(merk)
  if (merk <= 1) return
  go to 20
30 merk = k
  k = left(merk)
  if (k /= 0) go to 30

end subroutine fpseno
