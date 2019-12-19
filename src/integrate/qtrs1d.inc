PRNM(trapz_func) (f, a, b, N) result(int)
NUMFOR_KINDF :: f
NUMFOR_KINDR :: int
real(dp), optional, intent(IN)   :: a, b
integer, optional, intent(IN)   :: N
real(dp)                         :: dh, a_, b_
integer                         :: N_, i
real(dp) :: xx

N_ = 100; a_ = 0._8; b_ = 1._8; 
IF (Present(N)) N_ = N
if ((Present(a) .and. .not. Present(b)) .or. (Present(b) .and. .not. Present(a))) then
  call print_msg("Error: Both integration limits must be present", 'trapz')
  return
end if
if (Present(a)) then
  a_ = a; b_ = b
end if

int = (NUMFOR_EVAL_F(F, a_) + NUMFOR_EVAL_F(F, b_)) / 2._8
if (N_ == 1) then
  int = int * (b_ - a_)
  return
end if

! Version with low-usage of memory
dh = (b_ - a_) / (N_ - 1)
do i = 1, N_ - 2
  xx = a_ + i * dh
  int = int + NUMFOR_EVAL_F(F, xx)
end do
int = int * dh

end PRNM(trapz_func)
!
! ------------------------------------------------------------------------------
! Trapezoid nonlinearly spaced points
PRNM(trapz_nonlin_func) (f, x) result(int)
implicit none
NUMFOR_KINDF :: f
NUMFOR_KINDR :: int
real(8), dimension(:) :: x
real(8)              :: dh
integer              :: L, i
L = size(x)
!
int = NUMFOR_EVAL_F(f, x(1)) * (x(2) - x(1)) + NUMFOR_EVAL_F(f, x(L)) * (x(L) - x(L - 1))
do i = 2, L - 1
  dh = (x(i + 1) - x(i - 1))
  int = int + NUMFOR_EVAL_F(f, x(i)) * dh
enddo
int = int / 2._8
end PRNM(trapz_nonlin_func)

! ------------------------------------------------------------------------------
! Trapezoid Sampled points
PRNM(trapz_sample) (y, a, b) result(int)
implicit none
NUMFOR_KINDR, dimension(:) :: y
NUMFOR_KINDR :: int
real(dp) :: a, b, dh
integer :: L
L = size(y)
dh = (b - a) / (L - 1)
int = (dh * (sum(y(2:L - 1)) + (y(1) + y(L)) / 2_dp))
end PRNM(trapz_sample)

PRNM(trapz_dh_sample) (y, dh) result(int)
implicit none
NUMFOR_KINDR, dimension(:) :: y
NUMFOR_KINDR :: int
real(dp) :: dh
integer :: L
L = size(y)
int = (dh * (sum(y(2:L - 1)) + (y(1) + y(L)) / 2._dp))
end PRNM(trapz_dh_sample)

PRNM(trapz_nonlin_sample) (y, x) result(int)
implicit none
NUMFOR_KINDR, dimension(:) :: y
NUMFOR_KINDR :: int
real(dp), intent(IN), dimension(size(y)) :: x
integer :: L
L = size(y)
int = y(1) * (x(2) - x(1)) + y(L) * (x(L) - x(L - 1)) / 2.0_dp
int = int + sum(y(2:L - 1) * (x(2:L) - x(1:L - 1)))
end PRNM(trapz_nonlin_sample)

!+-----------------------------------------------------------------+
! Simpson's rule for 1d integration
!+-----------------------------------------------------------------+
! function version
PRNM(simps_func) (f, a, b, N) result(int)
implicit none
NUMFOR_KINDF :: f
NUMFOR_KINDR, dimension(:), allocatable :: y
NUMFOR_KINDR :: int
real(dp), optional, intent(IN) :: a, b
integer, optional, intent(IN)  :: N
real(dp) :: dh, a_, b_
integer  :: N_, i
real(dp), dimension(:), allocatable :: xx, wt

N_ = 100; IF (Present(N)) N_ = N
a_ = 0._dp; IF (Present(a)) a_ = a
b_ = 1._dp; IF (Present(b)) b_ = b
!
allocate (xx(N_), wt(N_), y(N_))
xx = linspace(a_, b_, N_, retstep=dh)

call get_simpson_weights(wt)

do i = 1, N_
  y(i) = NUMFOR_EVAL_F(f, xx(i))
end do

int = sum(y * wt) * dh
deallocate (xx, y, wt)
end PRNM(simps_func)

!+-----------------------------------------------------------------+
PRNM(simpson_dx_sample) (y, dx) result(int)
NUMFOR_KINDR, dimension(:), intent(IN) :: y
real(dp), dimension(:), intent(IN) :: dx
NUMFOR_KINDR :: int
integer :: L, i
L = size(y)
int = 0
i = 0
do while (i < L)
  if ((dx(i) == dx(i + 1)) .AND. (dx(i) == dx(i + 2))) then !Simpson's 3/8 rule
    int = int + &
      & (3._dp * dx(i) * (y(i) + 3._dp * (y(i + 1) + y(i + 2)) + y(i + 3))) / 8._dp
    i = i + 3
  elseif (dx(i) == dx(i + 1)) then      !Simpson's 1/3 rule

    int = int + &
      &(2._dp * dx(i) * (y(i) + 4._dp * y(i + 1) + y(i + 2))) / 6._dp
    i = i + 2
  elseif (dx(i) /= dx(i + 1)) then      !trapezoidal rule
    int = int + dx(i) * (y(i) + y(i + 1)) / 2._dp
    i = i + 1
  endif
enddo
end PRNM(simpson_dx_sample)

PRNM(simps_nonlin_func) (f, x) result(int)
NUMFOR_KINDF :: f
real(dp), intent(IN), dimension(:) :: x
NUMFOR_KINDR :: int
NUMFOR_KINDR, dimension(size(x)) :: y
integer :: L, i
real(dp), dimension(size(x)) :: dx
L = size(y)
dx(1:L - 1) = x(2:L) - x(1:L - 1)
do i = 1, size(x)
  y(i) = NUMFOR_EVAL_F(f, x(i))
end do
! Given arrays y and dx -> make calculation
int = simpson_dx_sample(y, dx)
end PRNM(simps_nonlin_func)

PRNM(simpson_sample) (y, a, b) result(int)
NUMFOR_KINDR, dimension(:), intent(IN) :: y
NUMFOR_KINDR :: int
real(dp), intent(IN) :: a, b
integer :: L
real(dp), dimension(size(y)) :: dx
L = size(y)
dx = (b - a) / real(L - 1, kind=dp)
int = simpson_dx_sample(y, dx)
end PRNM(simpson_sample)

PRNM(simpson_nonlin_sample) (y, x) result(int)
NUMFOR_KINDR, dimension(:), intent(IN) :: y
NUMFOR_KINDR :: int
real(dp), intent(IN), dimension(:) :: x
real(dp), dimension(size(x)) :: dx
integer :: L
L = size(y)
dx(1:L - 1) = x(2:L) - x(1:L - 1)
int = simpson_dx_sample(y, dx)
end PRNM(simpson_nonlin_sample)

PRNM(simpson_dh_sample) (y, dh) result(int)
NUMFOR_KINDR, dimension(:), intent(IN) :: y
NUMFOR_KINDR :: int
real(dp), intent(IN) :: dh
real(dp), dimension(size(y)) :: dx
dx = dh
int = simpson_dx_sample(y, dx)
end PRNM(simpson_dh_sample)
