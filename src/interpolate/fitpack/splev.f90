! JF: Renamed splev => dsplev to avoid conflict with wrapper
! Decided keep splev name in wrapper for compatibility with scipy
! IMPORTANT: NOT LONGER IN USE. Calling splder with nu=0
subroutine dsplev(t, n, c, k, x, y, m, e, ier)
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
  implicit none
  integer, intent(IN) :: n !<
  integer, intent(IN) :: k !<
  integer, intent(IN) :: m !<
  integer, intent(IN) :: e !<
  integer, intent(OUT) :: ier !<

  !  ..array arguments..
  real(8), dimension(n), intent(IN) :: t !<
  real(8), dimension(n), intent(IN) :: c !<
  real(8), dimension(m), intent(IN) :: x !<
  real(8), dimension(m), intent(OUT) :: y !<

  !  ..local scalars..
  integer :: i, k1, l, nk1
  !..++
  real(8) :: arg, tb, te
  !  ..local array..
  real(8) :: h(20)
  !  ..
  !  before starting computations a data check is made. if the input data
  !  are invalid control is immediately repassed to the calling program.
  ier = 10
  !++..
  if (m < 1) return
  !..++
  !--  10  do i=2,m
  !--        if(x(i) < x(i-1)) return
  !--      end do
  ier = 0
  !  fetch tb and te, the boundaries of the approximation interval.
  k1 = k + 1
  !
  nk1 = n - k1
  tb = t(k1)
  te = t(nk1 + 1)
  l = k1
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
    !++..
    !JF: This should be quite efficient if x(i) are ordered.
    do while (arg < t(l) .and. l /= k1)
      l = l - 1
    end do
    do while (arg >= t(l + 1) .and. l /= nk1)
      l = l + 1
    end do

    !  evaluate the non-zero b-splines at arg.
    call fpbspl(t, n, k, arg, l, h)
    !  find the value of s(x) at x=arg.

    y(i) = sum(c(l - k1 + 1:l) * h(:k1))
  end do

end subroutine dsplev

