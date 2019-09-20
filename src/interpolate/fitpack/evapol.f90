      real(8) function evapol(tu, nu, tv, nv, c, rad, x, y)
!  function program evacir evaluates the function f(x,y) = s(u,v),
!  defined through the transformation
!      x = u*rad(v)*cos(v)    y = u*rad(v)*sin(v)
!  and where s(u,v) is a bicubic spline ( 0<=u<=1 , -pi<=v<=pi ), given
!  in its standard b-spline representation.
!
!  calling sequence:
!     f = evapol(tu,nu,tv,nv,c,rad,x,y)
!
!  input parameters:
!   tu    : real array, length nu, which contains the position of the
!           knots in the u-direction.
!   nu    : integer, giving the total number of knots in the u-direction
!   tv    : real array, length nv, which contains the position of the
!           knots in the v-direction.
!   nv    : integer, giving the total number of knots in the v-direction
!   c     : real array, length (nu-4)*(nv-4), which contains the
!           b-spline coefficients.
!   rad   : real function subprogram, defining the boundary of the
!           approximation domain. must be declared external in the
!           calling (sub)-program
!   x,y   : real values.
!           before entry x and y must be set to the co-ordinates of
!           the point where f(x,y) must be evaluated.
!
!  output parameter:
!   f     : real
!           on exit f contains the value of f(x,y)
!
!  other subroutines required:
!    bispev,fpbisp,fpbspl
!
!  references :
!    de boor c : on calculating with b-splines, j. approximation theory
!                6 (1972) 50-62.
!    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
!                applics 10 (1972) 134-149.
!    dierckx p. : curve and surface fitting with splines, monographs on
!                 numerical analysis, oxford university press, 1993.
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
        integer :: nu, nv
        real(8) :: x, y
!  ..array arguments..
        real(8) :: tu(nu), tv(nv), c((nu - 4) * (nv - 4))
!  ..user specified function
        real(8) :: rad
!  ..local scalars..
        integer :: ier
        real(8) :: u, v, r, f, one, dist
!  ..local arrays
        real(8) :: wrk(8)
        integer :: iwrk(2)
!  ..function references
        real(8) :: atan2, sqrt
!  ..
!  calculate the (u,v)-coordinates of the given point.
        one = 1
        u = 0.
        v = 0.
        dist = x**2 + y**2
        if (dist <= 0.) go to 10
        v = atan2(y, x)
        r = rad(v)
        if (r <= 0.) go to 10
        u = sqrt(dist) / r
        if (u > one) u = one
!  evaluate s(u,v)
10      call bispev(tu, nu, tv, nv, c, 3, 3, u, 1, v, 1, f, wrk, 8, iwrk, 2, ier)
        evapol = f
        return
      end

