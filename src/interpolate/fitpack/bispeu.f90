subroutine bispeu(tx, nx, ty, ny, c, kx, ky, x, y, z, m, wrk, lwrk, ier)
  !  subroutine bispeu evaluates on a set of points (x(i),y(i)),i=1,...,m
  !  a bivariate spline s(x,y) of degrees kx and ky, given in the
  !  b-spline representation.
  !
  !  calling sequence:
  !     call bispeu(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk,
  !    * iwrk,kwrk,ier)
  !
  !  input parameters:
  !   tx    : real array, length nx, which contains the position of the
  !           knots in the x-direction.
  !   nx    : integer, giving the total number of knots in the x-direction
  !   ty    : real array, length ny, which contains the position of the
  !           knots in the y-direction.
  !   ny    : integer, giving the total number of knots in the y-direction
  !   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
  !           b-spline coefficients.
  !   kx,ky : integer values, giving the degrees of the spline.
  !   x     : real array of dimension (mx).
  !   y     : real array of dimension (my).
  !   m     : on entry m must specify the number points. m >= 1.
  !   wrk   : real array of dimension lwrk. used as workspace.
  !   lwrk  : integer, specifying the dimension of wrk.
  !           lwrk >= kx+ky+2
  !
  !  output parameters:
  !   z     : real array of dimension m.
  !           on successful exit z(i) contains the value of s(x,y)
  !           at the point (x(i),y(i)), i=1,...,m.
  !   ier   : integer error flag
  !    ier=0 : normal return
  !    ier=10: invalid input data (see restrictions)
  !
  !  restrictions:
  !   m >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
  !   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
  !   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
  !
  !  other subroutines required:
  !    fpbisp,fpbspl
  !
  !  ..scalar arguments..
  integer :: nx, ny, kx, ky, m, lwrk, ier
  !  ..array arguments..
  real(8) :: tx(nx), ty(ny), c((nx - kx - 1) * (ny - ky - 1)), x(m), y(m), z(m), wrk(lwrk)
  !  ..local scalars..
  integer :: iwrk(2)
  integer :: i, lwest
  !  ..
  !  before starting computations a data check is made. if the input data
  !  are invalid control is immediately repassed to the calling program.
  ier = 10
  lwest = kx + ky + 2
  if (lwrk < lwest) return
  if (m < 1) return
  ier = 0
  do i = 1, m
    call fpbisp(tx, nx, ty, ny, c, kx, ky, x(i), 1, y(i), 1, z(i), wrk(1), wrk(kx + 2), iwrk(1), iwrk(2))
  end do

end subroutine bispeu
