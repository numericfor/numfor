!! QMOMO computes modified Chebyshev moments.
!!
!!  This routine computes modified Chebyshev moments. The K-th modified Chebyshev moment is defined
!!  as the integral over (-1,1) of W(X)*T(K,X), where T(K,X) is the Chebyshev polynomial of degree K.
!!  It evaluates: rw(k) = integrate(w(x) * T(k-1,x), -1,1) for k=1,25:
!!  \f\[\int_{-1}^{1} w(x)\, T(k-1, x) \, dx \f\]
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
subroutine qmomo(alfa, beta, ri, rj, rg, rh, flgm)
  !*****************************************************************************80
  !  Parameters:
  !
  !    Input, real(8) :: ALFA, a parameter in the weight function w(x), ALFA > -1.
  !
  !    Input, real(8) :: BETA, a parameter in the weight function w(x), BETA > -1.
  !
  !           ri     - real
  !                    vector of dimension 25
  !                    ri(k) is the integral over (-1,1) of
  !                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
  !
  !           rj     - real
  !                    vector of dimension 25
  !                    rj(k) is the integral over (-1,1) of
  !                    (1-x)**beta*t(k-1,x), k = 1, ..., 25.
  !
  !           rg     - real
  !                    vector of dimension 25
  !                    rg(k) is the integral over (-1,1) of
  !                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ...,25.
  !
  !           rh     - real
  !                    vector of dimension 25
  !                    rh(k) is the integral over (-1,1) of
  !                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
  !
  !           flgm   - integer
  !                    input parameter indicating the modified moments
  !                    to be computed
  !                    flgm = 1 compute ri, rj
  !                         = 2 compute ri, rj, rg
  !                         = 3 compute ri, rj, rh
  !                         = 4 compute ri, rj, rg, rh
  !
  implicit none

  real(8), intent(IN) :: alfa !< parameter (exponent power) in the weight function w(x)
  real(8), intent(IN) :: beta !< parameter (exponent power) in the weight function w(x)
  real(8), dimension(25), intent(OUT) :: ri !< momenta for w(x) = (1+x)**alfa
  real(8), dimension(25), intent(OUT) :: rj !< momenta for w(x) = (1-x)**beta
  real(8), dimension(25), intent(OUT) :: rg !< momenta for w(x) = (1+x)**alfa*log((1+x)/2)
  real(8), dimension(25), intent(OUT) :: rh !< momenta for w(x) = (1-x)**alfa*log((1-x)/2)
  integer flgm !< Flag indicating the modified momenta to be computed:
  !!  1: compute ri, rj
  !!  2: compute ri, rj, rg
  !!  3: compute ri, rj, rh
  !!  4: compute ri, rj, rg, rh

  real(8) :: alfp1
  real(8) :: alfp2
  real(8) :: an
  real(8) :: anm1
  real(8) :: betp1
  real(8) :: betp2
  integer i
  integer im1
  real(8) :: ralf
  real(8) :: rbet
  !
  alfp1 = alfa + 1.0_8
  betp1 = beta + 1.0_8
  alfp2 = alfa + 2.0_8
  betp2 = beta + 2.0_8
  ralf = 2.0_8**alfp1
  rbet = 2.0_8**betp1
  !
  !  Compute RI, RJ using a forward recurrence relation.
  !
  ri(1) = ralf / alfp1
  rj(1) = rbet / betp1
  ri(2) = ri(1) * alfa / alfp2
  rj(2) = rj(1) * beta / betp2
  an = 2.0_8
  anm1 = 1.0_8

  do i = 3, 25
    ri(i) = -(ralf + an * (an - alfp2) * ri(i - 1)) / (anm1 * (an + alfp1))
    rj(i) = -(rbet + an * (an - betp2) * rj(i - 1)) / (anm1 * (an + betp1))
    anm1 = an
    an = an + 1.0_8
  end do

  if ((flgm == 2) .or. (flgm == 4)) then    !  Compute RG using a forward recurrence relation.
    rg(1) = -ri(1) / alfp1
    rg(2) = -(ralf + ralf) / (alfp2 * alfp2) - rg(1)
    an = 2.0_8
    anm1 = 1.0_8
    im1 = 2

    do i = 3, 25
      rg(i) = -(an * (an - alfp2) * rg(im1) - an * ri(im1) + anm1 * ri(i)) / &
              (anm1 * (an + alfp1))
      anm1 = an
      an = an + 1.0_8
      im1 = i
    end do
  end if

  if ((flgm == 3) .or. (flgm == 4)) then ! Compute RH using a forward recurrence relation.
    rh(1) = -rj(1) / betp1
    rh(2) = -(rbet + rbet) / (betp2 * betp2) - rh(1)
    an = 2.0_8
    anm1 = 1.0_8
    im1 = 2

    do i = 3, 25
      rh(i) = -(an * (an - betp2) * rh(im1) - an * rj(im1) + &
                anm1 * rj(i)) / (anm1 * (an + betp1))
      anm1 = an
      an = an + 1.0_8
      im1 = i
    end do
    rh(2:25:2) = -rh(2:25:2)
  end if

  rj(2:25:2) = -rj(2:25:2)

end subroutine qmomo

!>  qsort maintains the order of a list of local error estimates.
!!    This routine maintains the descending ordering in the list of the
!!    local error estimates resulting from the interval subdivision process.
!!    At each call two error estimates are inserted using the sequential
!!    search top-down for the largest error estimate and bottom-up for the
!!    smallest error estimate.
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
subroutine qsort(limit, last, maxerr, ermax, elist, iord, nrmax)
  implicit none
  integer, intent(IN) :: limit !<
  integer, intent(IN) :: last !<
  integer, intent(INOUT) :: maxerr !< The index in the list of the NRMAX-th largest error.
  real(dp), intent(OUT) :: ermax !< The NRMAX-th largest error = ELIST(MAXERR).
  real(dp), dimension(last), intent(IN) :: elist !< contains the error estimates
  integer, dimension(last), intent(INOUT) :: iord !< The first k elements contain
  !!    pointers to the error estimates such that elist(iord(1)) through elist(iord(k))
  !!    form a decreasing sequence, with:
  !!      - K = LAST, if LAST <= (LIMIT/2+2),
  !!      - K = LIMIT+1-LAST,  otherwise
  integer, intent(INOUT) :: nrmax !<
  !! Examples:
  !!
  real(dp) :: errmin, errmax
  integer i, j, k
  integer ibeg
  integer isucc
  integer jbnd
  integer jupbn
  !
  !  Check whether the list contains more than two error estimates.
  !
  if (last <= 2) then
    iord(1) = 1
    iord(2) = 2
    go to 90
  end if
  !
  !  This part of the routine is only executed if, due to a
  !  difficult integrand, subdivision increased the error
  !  estimate. in the normal case the insert procedure should
  !  start after the nrmax-th largest error estimate.
  !
  errmax = elist(maxerr)

  do i = 1, nrmax - 1

    isucc = iord(nrmax - 1)

    if (errmax <= elist(isucc)) then
      exit
    end if

    iord(nrmax) = isucc
    nrmax = nrmax - 1

  end do
  !
  !  Compute the number of elements in the list to be maintained
  !  in descending order.  This number depends on the number of
  !  subdivisions still allowed.
  !
  jupbn = last

  if ((limit / 2 + 2) < last) then
    jupbn = limit + 3 - last
  end if

  errmin = elist(last)
  !
  !  Insert errmax by traversing the list top-down, starting
  !  comparison from the element elist(iord(nrmax+1)).
  !
  jbnd = jupbn - 1
  ibeg = nrmax + 1

  do i = ibeg, jbnd
    isucc = iord(i)
    if (elist(isucc) <= errmax) then
      go to 60
    end if
    iord(i - 1) = isucc
  end do

  iord(jbnd) = maxerr
  iord(jupbn) = last
  go to 90
  !
  !  Insert errmin by traversing the list bottom-up.
  !
60 continue

  iord(i - 1) = maxerr
  k = jbnd

  do j = i, jbnd
    isucc = iord(k)
    if (errmin < elist(isucc)) then
      go to 80
    end if
    iord(k + 1) = isucc
    k = k - 1
  end do

  iord(i) = last
  go to 90

80 continue

  iord(k + 1) = last
  !
  !  Set maxerr and ermax.
90 continue

  maxerr = iord(nrmax)
  ermax = elist(maxerr)

end subroutine qsort

!> qwgtc computes the weight function  1 / ( x - c ) used by qc25c
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
function qwgtc(x, c, flgw) result(y)
  implicit none
  real(dp) :: y !< Value of weight function at x
  real(dp), intent(IN) :: x !<
  real(dp), dimension(:), intent(IN) :: c !< First component is c
  integer, intent(IN) :: flgw             !< Not used here
  y = flgw                     ! Dummy statement
  y = 1._8 / (x - c(1))
end function qwgtc

!> qwgtc computes the weight function used by qc25o
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
function qwgto(x, c, flgw) result(y)
  implicit none
  real(dp) :: y !< Value of weight function at x
  real(dp), intent(IN) :: x !< The point at which the weight function is evaluated.
  real(dp), dimension(:), intent(IN) :: c !< First component is omega (omega = c(1))
  integer, intent(IN) :: flgw             !< Specifies which weight function is used:
  !! 1: -> w(x) = cos(omega x)
  !! 2: -> w(x) = sin(omega x)
  if (flgw == 1) then
    y = cos(c(1) * x)
  else
    y = sin(c(1) * x)
  end if
end function qwgto

!> qwgts computes the weight function used by qc25s
!!
!!  The weight function depends on the value of `flgw`:
!!  flgw =
!!        1: => w(x) = (x-a)**alfa * (b-x)**beta
!!        2: => w(x) = (x-a)**alfa * (b-x)**beta * log (x-a)
!!        3: => w(x) = (x-a)**alfa * (b-x)**beta * log (b-x)
!!        4: => w(x) = (x-a)**alfa * (b-x)**beta * log (x-a) * log(b-x)
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
function qwgts(x, c, flgw) result(y)
  implicit none
  real(dp) :: y !< Value of weight function at x
  real(dp), intent(IN) :: x !< The point at which the weight function is evaluated.
  real(dp), dimension(:), intent(IN) :: c !< Gives the values for the parameters to use:
  !! c(1)= a    (lower limit of integration)
  !! c(2)= b    (upper limit of integration)
  !! c(3)= alfa (exponent)
  !! c(4)= beta (exponent)
  integer, intent(IN) :: flgw   !< Specifies which weight function is used

  associate (a=>c(1), b=>c(2), alfa=>c(3), beta=>c(4))
    if ((x == a) .or. (x == b)) then
      y = 0._dp
      return
    end if
    if (alfa == 0._8) then
      y = 1._dp
    else
      y = (x - a)**alfa
    end if
    IF (beta /= 0._8) y = y * (b - x)**beta
    select case (flgw)
    case (2)
      y = y * log(x - a)
    case (3)
      y = y * log(b - x)
    case (4)
      y = y * log(x - a) * log(b - x)
    end select
  end associate

end function qwgts
