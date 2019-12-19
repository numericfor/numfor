!! @date "2019-12-18 23:09:03 fiol"

!> Integration by Adaptive Simpson method of a function
!! \f[ \int_{a}^{b} f(t)\, dt \f]
#ifdef NF_FARG
PRNM(iads) (ff, a, b, args, IntVal, epsabs, epsrel, abserr, neval, ier)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(iads) (f, a, b, Intval, epsabs, epsrel, abserr, neval, ier)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
NUMFOR_KINDR, intent(OUT)  :: IntVal    !< Result
real(dp), intent(IN) :: a      !< Lower limit of integration
real(dp), intent(IN) :: b      !< Upper limit of integration
real(dp), optional, intent(IN) :: epsabs
real(dp), optional, intent(IN) :: epsrel !<  Bound to estimation error desired
real(dp), optional, intent(OUT) :: abserr !< Estimated error
integer, optional, intent(OUT) :: neval !< Number of function evaluations
integer, optional, intent(OUT) :: ier !< Output (error) code

integer, parameter :: maxK = 36
integer :: k, j
real(dp) :: xM, Error, eps
NUMFOR_KINDR :: CR, CRR
real(dp), dimension(maxK + 1) :: xL, xU
NUMFOR_KINDR, dimension(maxK + 1) :: fL, fM, fU
NUMFOR_KINDR, dimension(maxK + 1) :: Int
! Used for swap regions
real(dp) :: swx
NUMFOR_KINDR :: swf
NUMFOR_KINDR :: correction
real(dp):: epsabs_, epsrel_, abserr_
logical :: ready
integer :: neval_

#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif

! Optional input arguments
epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
epsrel_ = DEFAULT_EPSREL; IF (Present(epsrel)) epsrel_ = epsrel
abserr_ = 0.0_8

!     ---- Initial values ----
j = 1; IntVal = Zero; abserr_ = Zero; neval_ = 0

xL(j) = a; fL(j) = NUMFOR_EVAL_F(f, a)
xU(j) = b; fU(j) = NUMFOR_EVAL_F(f, b)
xM = (b + a) / 2; fM(j) = NUMFOR_EVAL_F(f, xM)
Int(j) = ((xU(j) - xL(j)) / 6) * (fL(j) + fU(j) + 4 * fM(j)) !Simpson

do
  neval_ = neval_ + 1
  k = j
  j = j + 1
  CR = Int(k)

  xU(j) = xU(k); fU(j) = fU(k)
  xU(k) = 0.5d0 * (xL(k) + xU(k)); fU(k) = fM(k)
  xL(j) = xU(k); fL(j) = fM(k)
  xM = 0.5d0 * (xU(k) + xL(k)); fM(k) = NUMFOR_EVAL_F(f, xM)
  Int(k) = ((xU(k) - xL(k)) / 6) * (fL(k) + fU(k) + 4 * fM(k)) !Simpson

  xM = 0.5d0 * (xU(j) + xL(j)); fM(j) = NUMFOR_EVAL_F(f, xM)
  Int(j) = ((xU(j) - xL(j)) / 6) * (fL(j) + fU(j) + 4 * fM(j)) !Simpson

  ! New value of integral on the interval is the sum of the value for the two subintervals
  CRR = Int(k) + Int(j)
  ! The correction to the integral is described in:
  ! "Numerical Analysis by Kincaid and Cheney, 3rd edition 2017. Section 7.5"
  correction = CRR + (CRR - CR) / 15._dp
  Error = abs(correction - CR)

  ! Next time bisect first the interval with higher contribution to the integral
  if (abs(Int(k)) > abs(Int(j))) then ! SWAP k <--> j
    swx = xU(j); xU(j) = xU(k); xU(k) = swx
    swx = xL(j); xL(j) = xL(k); xL(k) = swx
    swf = fU(j); fU(j) = fU(k); fU(k) = swf
    swf = fL(j); fL(j) = fL(k); fL(k) = swf
    swf = fM(j); fM(j) = fM(k); fM(k) = swf
    swf = Int(j); Int(j) = Int(k); Int(k) = swf
  end if

  ! Evaluate the exit condition
  eps = max(epsabs_, epsrel_ * abs(IntVal + correction))
  ready = (Error <= eps) .or. & ! Converged within tolerance
    & ((100 * abs(correction) < abs(IntVal)) .and. Error <= 10 * eps) ! The contribution is not significant

  if (ready) then ! End of iteration in k for a subinterval
    IntVal = IntVal + correction
    j = j - 2
    abserr_ = abserr_ + Error
    if (j <= 0) then        ! All intervals converged
      neval_ = 3 + 2 * neval_
      IF (Present(ier)) ier = 0
      IF (Present(abserr)) abserr = abserr_
      IF (Present(neval)) neval = neval_
      return
    end if
  else
    if (j + 1 > maxK) then
      IntVal = IntVal + correction
      neval_ = -(3 + 2 * neval_) ! A negative value means that it did not converge.
      IF (Present(ier)) ier = 1
      IF (Present(neval)) neval = neval_
      call print_msg("Maximum number of subdivisions reached for interval ["//&
        & str(xL(j))//','//str(xU(j))//"]", 'iads', errcode=1)
      return
    end if
  end if
enddo
end PRNM(iads)

!> iadsi Integrate(f,x,a,inf) by means of adaptive simpson method
!! It is assumed than the main contribution comes from a subinterval dx
!!
#ifdef NF_FARG
PRNM(iadsi) (ff, a, brkpts, args, IntVal, epsabs, epsrel, abserr, neval, ier)
implicit none
NUMFOR_KINDF :: ff         !< function to integrate
real(dp), dimension(:), intent(IN) :: args
NF_FARG :: f
#else
PRNM(iadsi) (f, a, brkpts, Intval, epsabs, epsrel, abserr, neval, ier)
implicit none
NUMFOR_KINDF :: f         !< function to integrate
#endif
NUMFOR_KINDR, intent(OUT) :: IntVal    !< Result
real(dp), intent(IN) :: a      !< Lower limit of integration
real(dp), dimension(:), intent(IN) :: brkpts      !< Integration domains will be [xL + (j-1)*dx, xL + j*dx]
real(dp), optional, intent(IN) :: epsabs
real(dp), optional, intent(IN) :: epsrel !<  Bound to estimation error desired
real(dp), optional, intent(OUT) :: abserr !< Estimated error
integer, optional, intent(OUT) :: neval !< Number of function evaluations
integer, optional, intent(OUT) :: ier

!! Examples:
!!
real(dp):: epsabs_, epsrel_, abserr_
real(dp) :: xu
real(dp) :: xl
real(dp) :: ERi, Ratio, Overall_error, ER
NUMFOR_KINDR :: R_i         ! Value of integral for each interval

integer :: i
integer :: ncalls
integer, parameter :: maxdX = 24
real(dp) :: dx
real(dp) :: lambda
integer :: npts, ier_, neval_

#ifdef NF_FARG
f%f => ff
allocate (f%args(size(args)))
f%args = args
#endif

! Optional input arguments
epsabs_ = DEFAULT_EPSABS; IF (Present(epsabs)) epsabs_ = epsabs
epsrel_ = DEFAULT_EPSREL; IF (Present(epsrel)) epsrel_ = epsrel
abserr_ = 0.0_8

npts = size(brkpts)
if (npts < 1) then
  call print_msg("At leas one break point must be given")
  return
end if

xl = a
IntVal = Zero; Overall_error = Zero
neval_ = 0

! First: A loop over critical regions
do i = 1, npts
  xu = xl + brkpts(i)        ! Primera integracion : [a, a + dx]
  call iads(f, xl, xu, R_i, epsabs_, epsrel_, ER, ncalls, ier_)
  IF (Present(ier)) ier = ier_
  IntVal = IntVal + R_i
  Ratio = abs(R_i) / (abs(IntVal) + Small)
  Overall_error = Overall_error + ER * Ratio
  neval_ = neval_ + ncalls
  xl = xu
  if (ier_ /= 0) then
    IF (Present(neval)) neval = neval_
    return
  end if
end do
ER = Overall_error          ! Error of the supposedly important part
dx = brkpts(npts) - a       ! Size of step to add
lambda = 1.0_dp             ! Scale parameter

! Loop over additional regions
! For each step the requirement on the accuracy is relaxed
ERi = max(epsrel_, 0.01_dp)
do i = npts, maxdX
  xl = xu; xu = xl + lambda * dx

  ! R_i = iads(f, xl, xu, ERi, abserr, ncalls, ier_)
  call iads(f, xl, xu, R_i, epsabs_, ERi, abserr, ncalls, ier_)
  IF (Present(ier)) ier = ier_

  if (ier_ /= 0) then
    call print_msg("Error in interval ["//str(xl)//","//str(xu)//"]", "iadsi")
    IF (Present(neval)) neval = ncalls ! return with error code message
    return
  end if

  neval_ = neval_ + ncalls
  IntVal = IntVal + R_i
  Ratio = abs(R_i) / (abs(IntVal) + Small)
  Overall_error = Overall_error + abserr * Ratio
  if (Ratio < epsrel_) then
    IF (Present(neval)) neval = neval_
    return
  end if

  if (Ratio < 0.1_dp) then
    ERi = ER / Ratio ! Relax the error condition
  else
    lambda = max(1.5 * lambda, 10._dp) ! Extends the new region
  end if
  ERi = min(ERi, 0.1_dp)
enddo
!  If reached here, no convergence was obtained
neval_ = -neval_
IF (Present(ier)) ier = 1
IF (Present(neval)) neval = neval_

end PRNM(iadsi)
