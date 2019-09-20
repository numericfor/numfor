      subroutine pardeu(tx, nx, ty, ny, c, kx, ky, nux, nuy, x, y, z, m,
      *wrk, lwrk, iwrk, kwrk, ier)
!  subroutine pardeu evaluates on a set of points (x(i),y(i)),i=1,...,m
!  the partial derivative ( order nux,nuy) of a bivariate spline
!  s(x,y) of degrees kx and ky, given in the b-spline representation.
!
!  calling sequence:
!     call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk,
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
!   nux   : integer values, specifying the order of the partial
!   nuy     derivative. 0<=nux<kx, 0<=nuy<ky.
!   kx,ky : integer values, giving the degrees of the spline.
!   x     : real array of dimension (mx).
!   y     : real array of dimension (my).
!   m     : on entry m must specify the number points. m >= 1.
!   wrk   : real array of dimension lwrk. used as workspace.
!   lwrk  : integer, specifying the dimension of wrk.
!           lwrk >= mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1)
!   iwrk  : integer array of dimension kwrk. used as workspace.
!   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
!
!  output parameters:
!   z     : real array of dimension (m).
!           on successful exit z(i) contains the value of the
!           specified partial derivative of s(x,y) at the point
!           (x(i),y(i)),i=1,...,m.
!   ier   : integer error flag
!   ier=0 : normal return
!   ier=10: invalid input data (see restrictions)
!
!  restrictions:
!   lwrk>=m*(kx+1-nux)+m*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1),
!
!  other subroutines required:
!    fpbisp,fpbspl
!
!  references :
!    de boor c : on calculating with b-splines, j. approximation theory
!                6 (1972) 50-62.
!   dierckx p. : curve and surface fitting with splines, monographs on
!                numerical analysis, oxford university press, 1993.
!
!  author :
!    p.dierckx
!    dept. computer science, k.u.leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1989
!
!  ..scalar arguments..
      integer :: nx, ny, kx, ky, m, lwrk, kwrk, ier, nux, nuy
!  ..array arguments..
      integer :: iwrk(kwrk)
      real(8) :: tx(nx), ty(ny), c((nx - kx - 1) * (ny - ky - 1)), x(m), y(m), z(m),
      *wrk(lwrk)
!  ..local scalars..
      integer :: i, iwx, iwy, j, kkx, kky, kx1, ky1, lx, ly, lwest, l1, l2, mm, m0, m1,
      *nc, nkx1, nky1, nxx, nyy
      real(8) :: ak, fac
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
      kx1 = kx + 1
      ky1 = ky + 1
      nkx1 = nx - kx1
      nky1 = ny - ky1
      nc = nkx1 * nky1
      if (nux < 0 .or. nux >= kx) go to 400
      if (nuy < 0 .or. nuy >= ky) go to 400
      lwest = nc + (kx1 - nux) * m + (ky1 - nuy) * m
      if (lwrk < lwest) go to 400
      if (kwrk < (m + m)) go to 400
      if (m < 1) go to 400
      ier = 0
      nxx = nkx1
      nyy = nky1
      kkx = kx
      kky = ky
!  the partial derivative of order (nux,nuy) of a bivariate spline of
!  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
!  we calculate the b-spline coefficients of this spline
      do 70 i = 1, nc
        wrk(i) = c(i)
70      continue
        if (nux == 0) go to 200
        lx = 1
        do 100 j = 1, nux
          ak = kkx
          nxx = nxx - 1
          l1 = lx
          m0 = 1
          do 90 i = 1, nxx
            l1 = l1 + 1
            l2 = l1 + kkx
            fac = tx(l2) - tx(l1)
            if (fac <= 0.) go to 90
            do 80 mm = 1, nyy
              m1 = m0 + nyy
              wrk(m0) = (wrk(m1) - wrk(m0)) * ak / fac
              m0 = m0 + 1
80            continue
90            continue
              lx = lx + 1
              kkx = kkx - 1
100           continue
200           if (nuy == 0) go to 300
              ly = 1
              do 230 j = 1, nuy
                ak = kky
                nyy = nyy - 1
                l1 = ly
                do 220 i = 1, nyy
                  l1 = l1 + 1
                  l2 = l1 + kky
                  fac = ty(l2) - ty(l1)
                  if (fac <= 0.) go to 220
                  m0 = i
                  do 210 mm = 1, nxx
                    m1 = m0 + 1
                    wrk(m0) = (wrk(m1) - wrk(m0)) * ak / fac
                    m0 = m0 + nky1
210                 continue
220                 continue
                    ly = ly + 1
                    kky = kky - 1
230                 continue
                    m0 = nyy
                    m1 = nky1
                    do 250 mm = 2, nxx
                      do 240 i = 1, nyy
                        m0 = m0 + 1
                        m1 = m1 + 1
                        wrk(m0) = wrk(m1)
240                     continue
                        m1 = m1 + nuy
250                     continue
!  we partition the working space and evaluate the partial derivative
300                     iwx = 1 + nxx * nyy
                        iwy = iwx + m * (kx1 - nux)
                        do 390 i = 1, m
                          call fpbisp(tx(nux + 1), nx - 2 * nux, ty(nuy + 1), ny - 2 * nuy, wrk, kkx, kky,
                          *x(i), 1, y(i), 1, z(i), wrk(iwx), wrk(iwy), iwrk(1), iwrk(2))
390                       continue
400                       return
                        end
