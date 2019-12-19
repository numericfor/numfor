  !> @file qnthsh.inc Implementation of tanh-sinh method
  !! @date "2019-12-10 11:55:21"

  !> qnthsh Integration by tanh-sinh method
  !!
  !! This routine is a simple non-adaptive automatic integrator using trapezoidal rule
  !! to integrate a the transformed integrand using the double exponential (also called tanh-sinh) method.
  !! The tahn-sinh scheme is based on the observation that the trapezoid rule converges
  !! very rapidly for functions on the entire real line that go to zero like exp( - exp(t) ).
  !! The change of variables x = tanh( pi sinh(t) /2) transforms an integral over [-1, 1]
  !! into an integral with integrand suited to the double exponential rule.
  !!
  !! The transformed integral is infinite, but we truncate the domain of integration to [-a, a],
  !! where `a <= 3`.
  !! Following J. D. Cook, we chose the limit '3' because for t = 3, the transformed x values
  !! are nearly equal to 1 for 12 or more significant figures.  Accordingly, the
  !! weights are 12 orders of magnitude smaller than the largest weights.
  !!
  !! The integration first applies the trapezoid rule to [-a, a].
  !! Then it subsequently cuts the step size in half each time, comparing the results.
  !! Integration stops when subsequent iterations are close enough together or the maximum
  !! integration points have been used.
  !! By cutting h in half, the previous integral can be reused; we only need evaluate the
  !! integrand at the newly added points.
  !! References:
  !! - https://en.wikipedia.org/wiki/Tanh-sinh_quadrature
  !! - Bailey, David H, "Tanh-Sinh High-Precision Quadrature", 2006 (http://crd.lbl.gov/~dhbailey/dhbpapers/dhb-tanh-sinh.pdf)
  !! - Cook, John D, ["Double Exponential Integration"](https://www.johndcook.com/blog/double_exponential_integration/)
  !!
#ifdef NF_FARG
  PRNM(qnthsh) (ff, a, b, args, IntVal, epsabs, epsrel, abserr, neval, ier)
  implicit none
  NUMFOR_KINDF :: ff         !< function to integrate
  real(dp), dimension(:), intent(IN) :: args !< Optional arguments
  NF_FARG :: f
#else
  PRNM(qnthsh) (f, a, b, IntVal, epsabs, epsrel, abserr, neval, ier)
  implicit none
  NUMFOR_KINDF :: f         !< function to integrate
#endif
  NUMFOR_KINDR :: IntVal    !< Result
  real(dp), intent(IN) :: a      !< Lower limit of integration
  real(dp), intent(IN) :: b      !< Upper limit of integration
  real(dp), optional, intent(IN) :: epsabs
  real(dp), optional, intent(IN) :: epsrel
  real(dp), optional, intent(OUT) :: abserr
  integer, optional, intent(OUT) :: neval
  integer, optional, intent(OUT) :: ier

  real(dp) :: yj, xj, xn, w, sh, chsh, h
  real(dp) :: t
  NUMFOR_KINDR :: contrib, newterm
  real(dp) :: alfa, beta
  real(dp) :: r
  integer :: k, j, n
  integer, parameter :: Nmax = 12
  real(dp):: epsabs_, epsrel_, abserr_

  ! We may use them for extrapolation (later)
  NUMFOR_KINDR, dimension(0:Nmax) :: Ival
  real(dp), dimension(0:Nmax) :: DIval

#ifdef NF_FARG
  f%f => ff
  allocate (f%args(size(args)))
  f%args = args
#endif

  ! Optional input arguments
  epsabs_ = 1.e-7_dp; IF (Present(epsabs)) epsabs_ = epsabs
  epsrel_ = 1.e-5_dp; IF (Present(epsrel)) epsrel_ = epsrel
  abserr_ = Zero

  ! Conversion from (a,b) to (-1,1)
  alfa = (b - a) / 2
  beta = (b + a) / 2

  n = 1
  ! One-point approximation
  w = M_PI_2
  j = 0; xj = beta
  IntVal = w * NUMFOR_EVAL_F(f, xj)
  h = max(4._8, min(6._8, 2.5 * log(-log(max(epsabs_, epsrel_)))))

  Ival(0) = IntVal * h * alfa
  DIval(0) = abs(Ival(0))

  do k = 1, Nmax
    newterm = Zero
    h = h / 2
    do j = 1, n, 2
      t = h * j
      sh = M_PI_2 * sinh(t); chsh = cosh(sh)
      yj = 1 / (exp(sh) * chsh)   ! yj = 1 - x_j
      xj = alfa * (1 - yj) + beta !
      xn = alfa * (yj - 1) + beta !
      w = M_PI_2 * cosh(t) / chsh**2
      contrib = w * (NUMFOR_EVAL_F(f, xj) + NUMFOR_EVAL_F(f, xn))
      newterm = newterm + contrib
    end do
    n = n * 2
    IntVal = IntVal + newterm
    Ival(k) = IntVal * h * alfa
    DIval(k) = abs(Ival(k) - Ival(k - 1))
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check condition and output

    r = log(DIval(k)) / log(DIval(k - 1))

    if (10 * abs(r - 2) < 1) then
      abserr_ = DIval(k)**2
    else
      abserr_ = DIval(k)
    end if

    if (abserr_ <= max(epsabs_, epsrel_ * abs(Ival(k)))) then
      IntVal = IntVal * alfa * h
      IF (Present(neval)) neval = 2 * (n / 2) + 1
      IF (Present(abserr)) abserr = abserr_
      IF (Present(ier)) ier = 0
      return
    end if
  end do
  call print_msg("Maximum number of subdivisions (Nmax = "//str(Nmax)//") reached.")
  IF (Present(ier)) ier = 1       ! Maximum number of subdivisions reached
  end PRNM(qnthsh)

! #ifdef NF_FARG
!   PRNM(qnthsh2) (ff, a, b, args, IntVal, epsabs, epsrel, abserr, neval, ier)
!   implicit none
!   NUMFOR_KINDF :: ff         !< function to integrate
!   real(dp), dimension(:), intent(IN) :: args
!   NF_FARG :: f
! #else
!   PRNM(qnthsh2) (f, a, b, IntVal, epsabs, epsrel, abserr, neval, ier)
!   implicit none
!   NUMFOR_KINDF :: f         !< function to integrate
! #endif

!   NUMFOR_KINDR :: IntVal    !< Result
!   real(dp), intent(IN) :: a      !< Lower limit of integration
!   real(dp), intent(IN) :: b      !< Upper limit of integration
!   real(dp), optional, intent(IN) :: epsabs
!   real(dp), optional, intent(IN) :: epsrel
!   real(dp), optional, intent(OUT) :: abserr
!   integer, optional, intent(OUT) :: neval
!   integer, optional, intent(OUT) :: ier

!   real(dp) :: yk, w, sh, chsh, h
!   real(dp) :: t, xM, xU!, xL
!   real(dp) :: alfa, beta
!   real(dp) :: r
!   integer :: l, k, s
!   integer :: n, lpts, upts
!   integer, parameter :: Nmax = 12
!   real(dp):: epsabs_, epsrel_, abserr_, err_

!   ! integer :: npts = 2
!   real(dp), dimension(0:2) :: brkpts
!   NUMFOR_KINDR, dimension(0:ubound(brkpts, dim=1) - 1) :: J
!   NUMFOR_KINDR, dimension(0:ubound(brkpts, dim=1)) :: fval

! ! ! We may use them for extrapolation (later)
!   NUMFOR_KINDR, dimension(0:Nmax) :: Ival
!   real(dp), dimension(0:Nmax) :: DIval

! #ifdef NF_FARG
!   f%f => ff
!   allocate (f%args(size(args)))
!   f%args = args
! #endif

! ! Optional input arguments
!   epsabs_ = 1.e-7_dp; IF (Present(epsabs)) epsabs_ = epsabs
!   epsrel_ = 1.e-5_dp; IF (Present(epsrel)) epsrel_ = epsrel
!   abserr_ = 0.0_dp

! ! Conversion from (a,b) to (-1,1)
!   alfa = (b - a) / 2
!   beta = (b + a) / 2

!   ! Upper limit of integration (approximation to infinity)
!   xU = max(4._8, min(6._8, 2.2 * log(-log(max(epsabs_, epsrel_))))) / 2
!   ! Break points
!   xM = 5 * xU / 8
!   ! brkpts = [0._dp, xU]
!   brkpts = [0._dp, xm, xU]

!   ! print *, brkpts

!   lpts = lbound(brkpts, dim=1)
!   upts = ubound(brkpts, dim=1)
! ! Values on Break-points
!   w = M_PI_2
!   fval(0) = w * NUMFOR_EVAL_F(f, beta) ! In zero we evaluate only once
!   do k = 1, upts
!     t = brkpts(k)
!     sh = M_PI_2 * sinh(t); chsh = cosh(sh)
!     yk = 1 / (exp(sh) * chsh)   ! yk = 1 - x_k
!     w = M_PI_2 * cosh(t) / chsh**2
!     xU = alfa * (1 - yk) + beta
!     xM = alfa * (yk - 1) + beta
!     fval(k) = w * (NUMFOR_EVAL_F(f, xM) + NUMFOR_EVAL_F(f, xU)) / 2
!     J(k - 1) = (fval(k - 1) + fval(k)) ! J_j
!   end do

!   IntVal = 0
!   neval = 1 + 2 * upts

!   do s = lpts, upts - 1
!     h = brkpts(s + 1) - brkpts(s)
!     Ival(0) = J(s) * h * alfa
!     DIval(0) = abs(Ival(0))
!     n = 2
!     do l = 1, Nmax
!       h = h / 2
!       call integr_level(brkpts(s), J(s))
!       if (err_ <= max(epsabs_, epsrel_ * abs(IntVal + Ival(l)))) then
!         IntVal = IntVal + Ival(l)
!         abserr_ = err_ + abserr_

!         IF (Present(neval)) neval = neval + 2 * (n / 2) + 1
!         ! print "(A)", "l="//str(l)//', neval='//str(neval)
!         ! print "(A)", "I(l)="//str(Ival(l))//', Integral='//str(IntVal)
!         IF (Present(abserr)) abserr = abserr_
!         exit
!       end if
!     end do
!     if (l >= Nmax) then
!       call print_msg("Maximum number of subdivisions (Nmax = "//str(Nmax)//") reached.")
!       IF (Present(ier)) ier = 1       ! Maximum number of subdivisions reached
!     end if
!   end do

! contains
!   subroutine integr_level(x0, J)
!     implicit none
!     real(dp), intent(IN) :: x0
!     NUMFOR_KINDR, intent(INOUT) :: J
!     integer :: k
!     real(dp) :: t, sh, yk, xp, xn
!     NUMFOR_KINDR :: newterm, contrib

!     newterm = 0._8
!     do k = 1, n, 2
!       t = x0 + h * k
!       sh = M_PI_2 * sinh(t); chsh = cosh(sh)
!       yk = 1 / (exp(sh) * chsh)
!       xp = alfa * (1 - yk) + beta !
!       xn = alfa * (yk - 1) + beta !
!       w = M_PI_2 * cosh(t) / chsh**2
!       contrib = w * (NUMFOR_EVAL_F(f, xp) + NUMFOR_EVAL_F(f, xn))
!       newterm = newterm + contrib
!       ! print "(A)", "t="//str(t)//", c="//str(contrib)//", n="//str(newterm)
!     end do
!     n = n * 2
!     J = J + newterm

!     Ival(l) = J * h * alfa
!     DIval(l) = abs(Ival(l) - Ival(l - 1))
!     ! print *, "Ival", n, Ival(l), Dival(l)
!     ! print *, ""
!     !
!     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ! Estimate the error
!     r = log(DIval(l)) / log(DIval(l - 1))

!     if (10 * abs(r - 2) < 1) then
!       err_ = DIval(l)**2
!     else
!       err_ = DIval(l)
!     end if
!     ! print *, 'err', 10 * abs(r - 2), abserr_
!   end subroutine integr_level

!   end PRNM(qnthsh2)
