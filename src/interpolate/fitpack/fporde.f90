subroutine fporde(x, y, m, kx, ky, tx, nx, ty, ny, nummer, index, nreg)
  !  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m
  !  according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
  !  to. for each panel a stack is constructed  containing the numbers
  !  of data points lying inside; index(j),j=1,2,...,nreg points to the
  !  first data point in the jth panel while nummer(i),i=1,2,...,m gives
  !  the number of the next data point in the panel.
  !  ..
  !  ..scalar arguments..
  integer :: m, kx, ky, nx, ny, nreg
  !  ..array arguments..
  real(8) :: x(m), y(m), tx(nx), ty(ny)
  integer :: nummer(m), index(nreg)
  !  ..local scalars..
  real(8) :: xi, yi
  integer :: im, k, kx1, ky1, k1, l, l1, nk1x, nk1y, num, nyy
  !  ..
  kx1 = kx + 1
  ky1 = ky + 1
  nk1x = nx - kx1
  nk1y = ny - ky1
  nyy = nk1y - ky

  index(:nreg) = 0

  do im = 1, m
    xi = x(im)
    yi = y(im)
    l = kx1
    l1 = l + 1
    do while (xi >= tx(l1) .and. l /= nk1x)
      l = l1
      l1 = l + 1
    end do
    k = ky1
    k1 = k + 1
    do while (yi >= ty(k1) .and. k /= nk1y)
      k = k1
      k1 = k + 1
    end do

    num = (l - kx1) * nyy + k - ky
    nummer(im) = index(num)
    index(num) = im
  end do

end subroutine fporde

