      subroutine splev(t, n, c, k, x, y, m, e, ier)
!  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
!  a spline s(x) of degree k, given in its b-spline representation.
!
!  calling sequence:
!     call splev(t,n,c,k,x,y,m,e,ier)
!
!  input parameters:
!    t    : array,length n, which contains the position of the knots.
!    n    : integer, giving the total number of knots of s(x).
!    c    : array,length n, which contains the b-spline coefficients.
!    k    : integer, giving the degree of s(x).
!    x    : array,length m, which contains the points where s(x) must
!           be evaluated.
!    m    : integer, giving the number of points where s(x) must be
!           evaluated.
!    e    : integer, if 0 the spline is extrapolated from the end
!           spans for points not in the support, if 1 the spline
!           evaluates to zero for those points, if 2 ier is set to
!           1 and the subroutine returns, and if 3 the spline evaluates
!           to the value of the nearest boundary point.
!
!  output parameter:
!    y    : array,length m, giving the value of s(x) at the different
!           points.
!    ier  : error flag
!      ier = 0 : normal return
!      ier = 1 : argument out of bounds and e == 2
!      ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!    m >= 1
!--    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
!
!  other subroutines required: fpbspl.
!
!  references :
!    de boor c  : on calculating with b-splines, j. approximation theory
!                 6 (1972) 50-62.
!    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
!                 applics 10 (1972) 134-149.
!    dierckx p. : curve and surface fitting with splines, monographs on
!                 numerical analysis, oxford university press, 1993.
!
!  author :
!    p.dierckx
!    dept. computer science, k.u.leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
!++ pearu: 11 aug 2003
!++   - disabled cliping x values to interval [min(t),max(t)]
!++   - removed the restriction of the orderness of x values
!++   - fixed initialization of sp to double precision value
!
!  ..scalar arguments..
        integer :: n, k, m, e, ier
!  ..array arguments..
        real(8) :: t(n), c(n), x(m), y(m)
!  ..local scalars..
        integer :: i, j, k1, l, ll, l1, nk1
!++..
        integer :: k2
!..++
        real(8) :: arg, sp, tb, te
!  ..local array..
        real(8) :: h(20)
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
        ier = 10
!--      if(m-1) 100,30,10
!++..
        if (m < 1) go to 100
!..++
!--  10  do 20 i=2,m
!--        if(x(i) < x(i-1)) go to 100
!--  20  continue
        ier = 0
!  fetch tb and te, the boundaries of the approximation interval.
        k1 = k + 1
!++..
        k2 = k1 + 1
!..++
        nk1 = n - k1
        tb = t(k1)
        te = t(nk1 + 1)
        l = k1
        l1 = l + 1
!  main loop for the different points.
        do 80 i = 1, m
!  fetch a new x-value arg.
          arg = x(i)
!  check if arg is in the support
          if (arg < tb .or. arg > te) then
            if (e == 0) then
              goto 35
            else if (e == 1) then
              y(i) = 0
              goto 80
            else if (e == 2) then
              ier = 1
              goto 100
            else if (e == 3) then
              if (arg < tb) then
                arg = tb
              else
                arg = te
              endif
            endif
          endif
!  search for knot interval t(l) <= arg < t(l+1)
!++..
35        if (arg >= t(l) .or. l1 == k2) go to 40
          l1 = l
          l = l - 1
          go to 35
!..++
40        if (arg < t(l1) .or. l == nk1) go to 50
          l = l1
          l1 = l + 1
          go to 40
!  evaluate the non-zero b-splines at arg.
50        call fpbspl(t, n, k, arg, l, h)
!  find the value of s(x) at x=arg.
          sp = 0.0d0
          ll = l - k1
          do 60 j = 1, k1
            ll = ll + 1
            sp = sp + c(ll) * h(j)
60          continue
            y(i) = sp
80          continue
100         return
          end
