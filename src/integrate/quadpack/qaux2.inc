!> @file qaux2.inc provides auxiliary routines for quadpack
!! @date "2019-11-19 16:32:13"
!!
!! Routines from QUADPACK, modified by J. Fiol to work with different kind of functions

!! QCHEB computes the Chebyshev series expansion.
!!
!!    This routine computes the Chebyshev series expansion
!!    of degrees 12 and 24 of a function using a fast Fourier transform method
!!
!!      f(x) = sum(k=1, ...,13) (cheb12(k)*t(k-1,x)),
!!      f(x) = sum(k=1, ...,25) (cheb24(k)*t(k-1,x)),
!!
!!    where T(K,X) is the Chebyshev polynomial of degree K.
!!
!!
!!  Author:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner
!!
!!  Reference:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
!!    QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983.
!!    (Section 2.2.3.2)
!!
PRNM(qcheb) (x, fval, cheb12, cheb24)
!*****************************************************************************80
!  Parameters:
!
!    Input, real(8) :: X(11), contains the values of COS(K*PI/24), for K = 1 to 11.
!
!    Input/output, real(8) :: FVAL(25), the function values at the points
!    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24, where (a,b) is the
!    approximation interval.  FVAL(1) and FVAL(25) are divided by two
!    These values are destroyed at output.
!
!    Output, real(8) :: CHEB12(13), the Chebyshev coefficients for degree 12.
!
!    Output, real(8) :: CHEB24(25), the Chebyshev coefficients for degree 24.
implicit none

real(8), intent(IN) :: x(11)
NUMFOR_KINDR, intent(OUT) :: cheb12(13)
NUMFOR_KINDR, intent(OUT) :: cheb24(25)
NUMFOR_KINDR, intent(INOUT) :: fval(25)
integer i
integer j
NUMFOR_KINDR :: alam
NUMFOR_KINDR :: alam1
NUMFOR_KINDR :: alam2
NUMFOR_KINDR :: part1
NUMFOR_KINDR :: part2
NUMFOR_KINDR :: part3
NUMFOR_KINDR :: v(12)

do i = 1, 12
  j = 26 - i
  v(i) = fval(i) - fval(j)
  fval(i) = fval(i) + fval(j)
end do

alam1 = v(1) - v(9)
alam2 = x(6) * (v(3) - v(7) - v(11))
cheb12(4) = alam1 + alam2
cheb12(10) = alam1 - alam2
alam1 = v(2) - v(8) - v(10)
alam2 = v(4) - v(6) - v(12)
alam = x(3) * alam1 + x(9) * alam2
cheb24(4) = cheb12(4) + alam
cheb24(22) = cheb12(4) - alam
alam = x(9) * alam1 - x(3) * alam2
cheb24(10) = cheb12(10) + alam
cheb24(16) = cheb12(10) - alam
part1 = x(4) * v(5)
part2 = x(8) * v(9)
part3 = x(6) * v(7)
alam1 = v(1) + part1 + part2
alam2 = x(2) * v(3) + part3 + x(10) * v(11)
cheb12(2) = alam1 + alam2
cheb12(12) = alam1 - alam2
alam = x(1) * v(2) + x(3) * v(4) + x(5) * v(6) + x(7) * v(8) &
       + x(9) * v(10) + x(11) * v(12)
cheb24(2) = cheb12(2) + alam
cheb24(24) = cheb12(2) - alam
alam = x(11) * v(2) - x(9) * v(4) + x(7) * v(6) - x(5) * v(8) &
       + x(3) * v(10) - x(1) * v(12)
cheb24(12) = cheb12(12) + alam
cheb24(14) = cheb12(12) - alam
alam1 = v(1) - part1 + part2
alam2 = x(10) * v(3) - part3 + x(2) * v(11)
cheb12(6) = alam1 + alam2
cheb12(8) = alam1 - alam2
alam = x(5) * v(2) - x(9) * v(4) - x(1) * v(6) &
       - x(11) * v(8) + x(3) * v(10) + x(7) * v(12)
cheb24(6) = cheb12(6) + alam
cheb24(20) = cheb12(6) - alam
alam = x(7) * v(2) - x(3) * v(4) - x(11) * v(6) + x(1) * v(8) &
       - x(9) * v(10) - x(5) * v(12)
cheb24(8) = cheb12(8) + alam
cheb24(18) = cheb12(8) - alam

do i = 1, 6
  j = 14 - i
  v(i) = fval(i) - fval(j)
  fval(i) = fval(i) + fval(j)
end do

alam1 = v(1) + x(8) * v(5)
alam2 = x(4) * v(3)
cheb12(3) = alam1 + alam2
cheb12(11) = alam1 - alam2
cheb12(7) = v(1) - v(5)
alam = x(2) * v(2) + x(6) * v(4) + x(10) * v(6)
cheb24(3) = cheb12(3) + alam
cheb24(23) = cheb12(3) - alam
alam = x(6) * (v(2) - v(4) - v(6))
cheb24(7) = cheb12(7) + alam
cheb24(19) = cheb12(7) - alam
alam = x(10) * v(2) - x(6) * v(4) + x(2) * v(6)
cheb24(11) = cheb12(11) + alam
cheb24(15) = cheb12(11) - alam

do i = 1, 3
  j = 8 - i
  v(i) = fval(i) - fval(j)
  fval(i) = fval(i) + fval(j)
end do

cheb12(5) = v(1) + x(8) * v(3)
cheb12(9) = fval(1) - x(8) * fval(3)
alam = x(4) * v(2)
cheb24(5) = cheb12(5) + alam
cheb24(21) = cheb12(5) - alam
alam = x(8) * fval(2) - fval(4)
cheb24(9) = cheb12(9) + alam
cheb24(17) = cheb12(9) - alam
cheb12(1) = fval(1) + fval(3)
alam = fval(2) + fval(4)
cheb24(1) = cheb12(1) + alam
cheb24(25) = cheb12(1) - alam
cheb12(13) = v(1) - v(3)
cheb24(13) = cheb12(13)
alam = 1.0_8 / 6.0_8

cheb12(2:12) = cheb12(2:12) * alam

alam = 0.5_8 * alam
cheb12(1) = cheb12(1) * alam
cheb12(13) = cheb12(13) * alam

cheb24 = cheb24 * alam
cheb24(1) = cheb24(1) * 0.5_8
cheb24(25) = cheb24(25) * 0.5_8

end PRNM(qcheb)

!! QEXTR carries out the Epsilon extrapolation algorithm.
!!
!!
!! The routine determines the limit of a given sequence of approximations, by means of the epsilon
!! algorithm of P. Wynn.  An estimate of the absolute error is given.  The condensed epsilon table
!! is computed.  Only those elements needed for the computation of the next diagonal are preserved.
!!
!!  Author:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner
!!
!!  Reference:
!!
!!    Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
!!    QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
!!
PRNM(qextr) (n, epstab, result, abserr, res3la, nres)
!*****************************************************************************80
!  Parameters:
!
!    Input, integer N, indicates the entry of EPSTAB which contains
!    the new element in the first column of the epsilon table.
!
!    Input/output, real(8) :: EPSTAB(52), the two lower diagonals of the triangular
!    epsilon table.  The elements are numbered starting at the right-hand
!    corner of the triangle.
!
!    Output, real(8) :: RESULT, the estimated value of the integral.
!
!    Output, real(8) :: ABSERR, estimate of the absolute error computed from
!    RESULT and the 3 previous results.
!
!    Output, real(8) :: RES3LA(3), the last 3 results.
!
!    Input/output, integer NRES, the number of calls to the routine.  This
!    should be zero on the first call, and is automatically updated
!    before return.
!
!  Local Parameters:
!
!           e0     - the 4 elements on which the
!           e1       computation of a new element in
!           e2       the epsilon table is based
!           e3                 e0
!                        e3    e1    new
!                              e2
!           newelm - number of elements to be computed in the new
!                    diagonal
!           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
!           result - the element in the new diagonal with least value
!                    of error
!           limexp is the maximum number of elements the epsilon table
!           can contain. if this number is reached, the upper diagonal
!           of the epsilon table is deleted.
!
implicit none

NUMFOR_KINDR, intent(INOUT) :: epstab(:)
NUMFOR_KINDR, intent(OUT) :: result
NUMFOR_KINDR, intent(OUT) :: res3la(3)
real(8) :: abserr
NUMFOR_KINDR :: delta1
NUMFOR_KINDR :: delta2
NUMFOR_KINDR :: delta3
real(8) :: epsinf
real(8) :: error
real(8) :: err1
real(8) :: err2
real(8) :: err3
NUMFOR_KINDR :: e0
NUMFOR_KINDR :: e1
NUMFOR_KINDR :: e2
NUMFOR_KINDR :: e3
real(8) :: e1abs
integer :: i
integer :: ib
integer :: ib2
integer :: ie
integer :: indx
integer :: k1
integer :: k2
integer :: k3
integer :: limexp
integer :: n
integer :: newelm
integer :: nres
integer :: num
NUMFOR_KINDR :: res

NUMFOR_KINDR :: ss
real(8) :: tol1
real(8) :: tol2
real(8) :: tol3

nres = nres + 1
abserr = huge(abserr)
result = epstab(n)

if (n < 3) then
  abserr = max(abserr, 0.5_8 * epsilon(abs(result)) * abs(result))
  return
end if

limexp = size(epstab) - 2
epstab(n + 2) = epstab(n)
newelm = (n - 1) / 2
epstab(n) = huge(abs(epstab(n)))
num = n
k1 = n

do i = 1, newelm

  k2 = k1 - 1
  k3 = k1 - 2
  res = epstab(k1 + 2)
  e0 = epstab(k3)
  e1 = epstab(k2)
  e2 = res
  e1abs = abs(e1)
  delta2 = e2 - e1
  err2 = abs(delta2)
  tol2 = max(abs(e2), e1abs) * epsilon(abs(e2))
  delta3 = e1 - e0
  err3 = abs(delta3)
  tol3 = max(e1abs, abs(e0)) * epsilon(abs(e0))
  !
  !  If e0, e1 and e2 are equal to within machine accuracy, convergence
  !  is assumed.
  !
  if (err2 <= tol2 .and. err3 <= tol3) then
    result = res
    abserr = err2 + err3
    abserr = max(abserr, 0.5_8 * epsilon(abs(result)) * abs(result))
    return
  end if

  e3 = epstab(k1)
  epstab(k1) = e1
  delta1 = e1 - e3
  err1 = abs(delta1)
  tol1 = max(e1abs, abs(e3)) * epsilon(abs(e3))
  !
  !  If two elements are very close to each other, omit a part
  !  of the table by adjusting the value of N.
  !
  if (err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3) go to 20

  ss = 1.0_8 / delta1 + 1.0_8 / delta2 - 1.0_8 / delta3
  epsinf = abs(ss * e1)
  !
  !  Test to detect irregular behavior in the table, and
  !  eventually omit a part of the table adjusting the value of N.
  !
  if (epsinf > 1.0e-4_8) go to 30

20 continue

  n = i + i - 1
  exit
  !
  !  Compute a new element and eventually adjust the value of RESULT.
  !
30 continue

  res = e1 + 1.0_8 / ss
  epstab(k1) = res
  k1 = k1 - 2
  error = err2 + abs(res - e2) + err3

  if (error <= abserr) then
    abserr = error
    result = res
  end if

end do
!
!  Shift the table.
!
if (n == limexp) then
  n = 2 * (limexp / 2) - 1
end if

if ((num / 2) * 2 == num) then
  ib = 2
else
  ib = 1
end if

ie = newelm + 1

do i = 1, ie
  ib2 = ib + 2
  epstab(ib) = epstab(ib2)
  ib = ib2
end do

if (num /= n) then

  indx = num - n + 1

  do i = 1, n
    epstab(i) = epstab(indx)
    indx = indx + 1
  end do

end if

if (nres < 4) then
  res3la(nres) = result
  abserr = huge(abserr)
else
  abserr = abs(result - res3la(3)) + abs(result - res3la(2)) &
           + abs(result - res3la(1))
  res3la(1) = res3la(2)
  res3la(2) = res3la(3)
  res3la(3) = result
end if

abserr = max(abserr, 0.5_8 * epsilon(abs(result)) * abs(result))

end PRNM(qextr)

