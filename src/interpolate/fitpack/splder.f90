subroutine splder(t, n, c, k, nu, x, y, m, e, wrk, ier)
  !  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
  !  the derivative of order nu of a spline s(x) of degree k,given in
  !  its b-spline representation.
  !
  !  calling sequence:
  !     call splder(t,n,c,k,nu,x,y,m,e,wrk,ier)
  !
  !  input parameters:
  !    t    : array,length n, which contains the position of the knots.
  !    n    : integer, giving the total number of knots of s(x).
  !    c    : array,length n, which contains the b-spline coefficients.
  !    k    : integer, giving the degree of s(x).
  !    nu   : integer, specifying the order of the derivative. 0<=nu<=k
  !    x    : array,length m, which contains the points where the deriv-
  !           ative of s(x) must be evaluated.
  !    m    : integer, giving the number of points where the derivative
  !           of s(x) must be evaluated
  !    e    : integer, if 0 the spline is extrapolated from the end
  !           spans for points not in the support, if 1 the spline
  !           evaluates to zero for those points, if 2 ier is set to
  !           1 and the subroutine returns, and if 3 the spline evaluates
  !           to the value of the nearest boundary point.
  !    wrk  : real array of dimension n. used as working space.
  !
  !  output parameters:
  !    y    : array,length m, giving the value of the derivative of s(x)
  !           at the different points.
  !    ier  : error flag
  !      ier = 0 : normal return
  !      ier = 1 : argument out of bounds and e == 2
  !      ier =10 : invalid input data (see restrictions)
  !
  !  restrictions:
  !    0 <= nu <= k
  !    m >= 1
  !    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
  !
  !  other subroutines required: fpbspl
  !
  !  references :
  !    de boor c : on calculating with b-splines, j. approximation theory
  !                6 (1972) 50-62.
  !    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
  !                applics 10 (1972) 134-149.
  !   dierckx p. : curve and surface fitting with splines, monographs on
  !                numerical analysis, oxford university press, 1993.
  !
  !  author :
  !    p.dierckx
  !    dept. computer science, k.u.leuven
  !    celestijnenlaan 200a, b-3001 heverlee, belgium.
  !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
  !
  !  latest update : march 1987
  !
  !++ pearu: 13 aug 2003
  !++   - disabled cliping x values to interval [min(t),max(t)]
  !++   - removed the restriction of the orderness of x values
  !++   - fixed initialization of sp to double precision value
  !
  ! JF: 2019. Converted to f20xx sintaxis and modified to use also as
  ! replacement of splev
  !  ..scalar arguments..
  implicit none
  integer, intent(IN) :: n !<
  integer, intent(IN) :: k !<
  integer, intent(IN) :: nu !<
  integer, intent(IN) :: m !<
  integer, intent(IN) :: e !<
  integer, intent(OUT) :: ier !<
  !  ..array arguments..
  real(8), dimension(n), intent(IN) :: t !<
  real(8), dimension(n), intent(IN) :: c !<
  real(8), dimension(m), intent(IN) :: x !<
  real(8), dimension(m), intent(OUT) :: y !<
  real(8), dimension(n), intent(INOUT) :: wrk !<

  !  ..local scalars..
  integer :: i, j, kk, k1, k2, l, l1, l2, nk1, nk2, nn
  real(8) :: ak, arg, fac, tb, te
  !++..
  integer :: k3
  !..++
  !  ..local arrays ..
  real(8) :: h(6)
  !  before starting computations a data check is made. if the input data
  !  are invalid control is immediately repassed to the calling program.
  ier = 10
  if (nu < 0 .or. nu > k) return
  !++..
  if (m < 1) return
  !..++
  !--  10  do i=2,m
  !--        if(x(i) < x(i-1)) return
  !--      end do
  ier = 0
  !  fetch tb and te, the boundaries of the approximation interval.
  k1 = k + 1
  k3 = k1 + 1
  nk1 = n - k1
  tb = t(k1)
  te = t(nk1 + 1)
  !  the derivative of order nu of a spline of degree k is a spline of
  !  degree k-nu,the b-spline coefficients wrk(i) of which can be found
  !  using the recurrence scheme of de boor.
  l = 1
  kk = k
  nn = n

  wrk(:nk1) = c(:nk1)

  if (nu > 0) then
    nk2 = nk1
    do j = 1, nu
      ak = kk
      nk2 = nk2 - 1
      l1 = l
      do i = 1, nk2
        l1 = l1 + 1
        l2 = l1 + kk
        fac = t(l2) - t(l1)
        IF (fac > 0._8) wrk(i) = ak * (wrk(i + 1) - wrk(i)) / fac
      end do

      l = l + 1
      kk = kk - 1
    end do
  end if

  if (kk == 0) then !  if nu=k the derivative is a piecewise constant function
    j = 1
  else
    l = k1
    k2 = k1 - nu
  end if

  !  main loop for the different points.
  do i = 1, m
    !  fetch a new x-value arg.
    arg = x(i)
    !  check if arg is in the support
    if (arg < tb .or. arg > te) then
      if (e == 3) then
        if (arg < tb) then
          arg = tb
        else
          arg = te
        endif
      else if (e == 2) then
        ier = 1
        return
      else if (e == 1) then
        y(i) = 0
        cycle
      endif
    endif
    !  search for knot interval t(l) <= arg < t(l+1)
    !JF: This should be quite efficient if x(i) are ordered.
    do while (arg < t(l) .and. l /= k1)
      l = l - 1
      j = j - 1
    end do
    do while (arg >= t(l + 1) .and. l /= nk1)
      l = l + 1
      j = j + 1
    end do
    !  evaluate the non-zero b-splines of degree k-nu at arg.
    if (kk == 0) then        !  nu=k: the derivative is a piecewise function
      y(i) = wrk(j)
    else
      call fpbspl(t, n, kk, arg, l, h)
      y(i) = sum(wrk(l - k1 + 1:l - nu) * h(:k2))
    end if
  end do
end subroutine splder
