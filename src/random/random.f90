!> @file random.f90 for random number generator (and distributions)
!> @author Juan Fiol <juanfiol@gmail.com> (modifications, see real authors below)
!! @date "2024-04-08 16:13:26"

!> \defgroup randomnumber Random number generator
!!
!! The random number are generated using the 64-bit version of MT19937,
!! translated from C-program for MT19937-64 (2004/9/29 version)
!! originally coded by Takuji Nishimura and Makoto Matsumoto
!! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
!!
!! Fortran translation by Rémi Piatek
!! see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/fortran.html
module random
  USE utils, only: i8, dp
  use mt19937_64, only: init_genrand64, init_by_array64,&
    & genrand64_real2, genrand64_real3, genrand64_int64
  ! JF: Added range of integer provided
  integer(i8), parameter :: rng_max = 2_i8**62 + 1_i8      ! 2^62+1
  integer(i8), parameter :: rng_min = ishft(-2_i8**62, 1)  ! -2^63

  !> Initialize the random number generator
  interface random_seed
    module procedure init_genrand64
    module procedure init_by_array64
    module procedure init_genrand64_int
    module procedure init_genrand64_empty
  end interface random_seed

  !> @ingroup randomnumber Function returning a real number in the semi-open interval [0,1)
  interface random_real
    module procedure genrand64_real2, rng_real_interval
  end interface random_real

  !> @ingroup randomnumber
  !! Function returning a real number in the open interval (0,1)
  interface random_real_pos
    module procedure genrand64_real3
  end interface random_real_pos

  !> @ingroup randomnumber Function returning an integer number in the open interval (0,1)
  interface random_int
    module procedure genrand64_int64, rng_integer_pos, rng_int
  end interface random_int

  private
  public :: rng_min, rng_max
  public :: random_seed, random_real, random_real_pos, random_int
contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!! random floating !!!!!!!!!!!!!!!!!!!!!!!!!!

  !> rng_real_interval Draw samples from a uniform distribution.
  !!
  !! Samples are uniformly distributed over the half-open interval [low, high) (includes low, but excludes high).
  function rng_real_interval(a, b) result(y)
    implicit none
    real(dp) :: y !< random float number
    real(dp), intent(IN) :: a !< lower limit
    real(dp), intent(IN) :: b !< higher limit
    y = (b - a) * random_real() + a
  end function rng_real_interval

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!! random integer !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @ingroup randomnumber rng_integer_pos returns a random integer in the range [0, N), or equivalently [0, N-1]
  !!
  !!  This function returns a random integer from 0 to N-1 inclusive
  !!  by scaling down and/or discarding samples.
  !!  All integers in the range [0,N-1] are produced with equal probability.
  !!  This function is designed for sampling from ranges smaller than the range of the underlying generator.
  !!  To use the full range use random_int() with no argument
  function rng_integer_pos(N) result(y)
    implicit none
    integer(i8), intent(IN) :: N !< Uper limit of the desired range
    integer(i8) :: y             !< integer random number generated
    !! Examples:
    !! ```{.f90}
    !! print *, rng_int()  ! Random integer in the range [-2^63, 2^62+1]
    !! print *, rng_int(6) ! Random integer in the range [0,5]
    !! ```
    integer(i8) :: scale

    ! Check N > 0 and N < range of generator
    y = 0_i8; IF ((N < 1) .OR. (N > rng_max)) return

    ! We divide first to avoid overflow (scale <= rng_max)
    scale = rng_max; IF (N /= 1) scale = (rng_max / N) * 2

    y = N + 1
    do while ((y >= N))
      y = random_int() / scale
      ! Due to integer division, zero gets twice the amount than the rest
      IF (y < 0) y = -y
    end do
  end function rng_integer_pos

  !> rng_int Computes a random integer in the range [low, high] inclusive
  function rng_int(low, high) result(y)
    implicit none
    integer(i8) :: y !< random integer
    integer(i8), intent(IN) :: low !< lower limit
    integer(i8), intent(IN) :: high !< higher limit
    integer(i8) :: N
    N = high - low + 1
    y = rng_integer_pos(N) + low
  end function rng_int
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!! Seed routines !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> init_genrand64_empty
  subroutine init_genrand64_empty()
    implicit none
    integer(i8) :: seed
    seed = generate_random_seed()
    call random_seed(seed)
  end subroutine init_genrand64_empty

  subroutine init_genrand64_int(seed)
    implicit none
    integer, intent(IN) :: seed
    call random_seed(int(seed, kind=i8))
  end subroutine init_genrand64_int

  !> generate_random_seed uses either /dev/random or the clock
  function generate_random_seed() result(y)
    implicit none
    integer(i8) :: y !< Seed generated
    integer :: u
    integer :: i
    integer :: ost, rst         ! status flags
    character(len=12), dimension(2), parameter :: rdev = ['/dev/random ', '/dev/urandom']
    do i = 1, size(rdev)
      open (newunit=u, file=rdev(i), access='stream', action='read', form='unformatted', iostat=ost)
      if (ost == 0) then
        read (u, iostat=rst) y
        IF (rst == 0) return
      end if
    end do
    ! If got here, it could not read a value. Get it from clock
    y = iclock()
  end function generate_random_seed

  !> Set a positive integer from the clock values
  !> @return a non-zero positive integer
  function iclock() result(ic)
    integer(i8) :: ic
    integer, dimension(8) :: t
    call date_and_time(values=t)
    ic = abs(t(3) + 10 * t(4) + 100 * t(5) + 1000 * t(6) + 10000 * t(7) + 100000 * t(8))
  end function iclock

end module random
! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:
