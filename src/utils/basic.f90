!> @file basic.f90 contains module basic
!! @date "2024-02-21 15:35:55"

!> @ingroup utils
!> This module will provide some basic convenience routines
!! [Description](@ref submodule-utils)
!!
!! Eventually it will provide:
!!  - [X] types (real, double precision, infinite, etc)
!!  - [X] Basic error handling (defining output for warning and errors)
!!    To improve maybe look into https://docs.python.org/3/library/warnings.html
!!  - Some basic file handling (some os.path functionality) (Really needed?)
!!  - [X] Some sys functionality (stdin, stdout, stderr via `iso_fortran_env`)
!!  - [X] Timers, timestamp (a simple implementation)
!!
!! For numerical work we need to use specific minimal precision which is
!! machine-and-compiler-independent. The real type definitions emphatize this.
!! The names are `sp` and `dp` are chosen only by tradition/historic reasons.
!! We define types with at least 15 decimal places and exponent range of 307.
!! (see for instance http://fortranwiki.org/fortran/show/Real+precision).
!! We also could use c_double from module iso_c_binding (fortran 2003 and later)
!!
module basic
  USE, intrinsic :: iso_fortran_env, only: stdin => input_unit, &
    stdout => output_unit, stderr => error_unit,&
    &  int64, real32, real64, real128

  implicit none

  ! JF: Solved in ad-hoc manner. More discussion/work is needed
  ! ! JF: An alternative should be provided, may be something like
  ! USE, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_negative,
  !   & ieee_value, &
  !   & ieee_negative_inf, ieee_positive_inf

  ! integer, parameter :: dp = selected_real_kind(15, 307)
  ! integer, parameter :: i8 = selected_int_kind(9, 307)
  ! integer, parameter :: sp = selected_real_kind(6, 37)
  ! integer, parameter :: qp = selected_real_kind(2 * precision(1._dp))

  ! integer, parameter :: i4 = int32
  integer, parameter :: i8 = int64
  integer, parameter :: dp = real64
  integer, parameter :: r8 = real64
  integer, parameter :: sp = real32
  integer, parameter :: qp = real128

  integer, private, dimension(7), parameter :: time_factors = [0, 12, 30, 24, 60, 60, 1000] ! Factors to convert to higher time-unit

  !> Simple timer. Holds start-time, stop-time, and time-difference
  type timer
    integer, dimension(7) :: start_date !< Start date
    integer, dimension(7) :: stop_date !< Stop date
    integer, dimension(7) :: date_elapsed !< Date-difference between start and stop
    real(dp), private :: start_cputime    !< Start cputime
    real(dp), private :: stop_cputime     !< Stop cputime
    real(dp), public :: elapsed !< Elapsed time (in seconds)
  contains

    procedure :: start => start_timer !< start time
    procedure :: stop => stop_timer   !< stop time
    procedure :: show => print_elapsed !< Show time
    procedure, private :: add_timers   !< add times
    procedure, private :: assign_int   !< assign int to time
    procedure, private :: assign       !< assign to time
    procedure, private :: divide_T_int !< divide time in intervals

    generic :: Operator(+) => add_timers   !< add times
    generic :: Operator(/) => divide_T_int         !< divide time in intervals
    generic :: Assignment(=) => assign_int, assign !< Assign a time

  end type timer

  private
  public :: i8, sp, dp, r8, qp, timer
  public :: nf_minf, nf_inf, is_inf
  public :: stdout, stdin, stderr
  private :: start_timer, stop_timer, print_elapsed
  public :: print_msg

  ! Common constants
  real(dp), public, parameter :: Zero = 0._dp          !< Real cero
  real(dp), public, parameter :: Small = epsilon(1._dp) !< Real small number
  complex(dp), public, parameter :: C_Z0 = (0._dp, 0._dp) !< Complex cero
  complex(dp), public, parameter :: C_R1 = (1._dp, 0._dp) !< Complex unit
  complex(dp), public, parameter :: C_I1 = (0._dp, 1._dp) !< Complex imaginary unit

  real(dp), parameter :: nf_minf = -huge(1._dp) !< Minus infinite
  real(dp), parameter :: nf_inf = huge(1._dp)   !< Positive infinite

  ! Basic mathematical constants (mostly from gsl).
  ! Many have more decimal places than needed/used
  real(dp), public, parameter :: M_pi = 3.14159265358979323846264338328_dp !< \f$ \pi \f$
  real(dp), public, parameter :: M_dpi = 6.28318530717958647692528676656_dp !< 2*pi
  real(dp), public, parameter :: M_pi_2 = 1.57079632679489661923132169164_dp !< pi/2
  real(dp), public, parameter :: M_pi_4 = 0.78539816339744830966156608458_dp !< pi/4
  real(dp), public, parameter :: M_sqrtpi = 1.77245385090551602729816748334_dp !< sqrt(pi)
  real(dp), public, parameter :: M_2_sqrtpi = 1.12837916709551257389615890312_dp !< 2/sqrt(pi)
  real(dp), public, parameter :: M_1_pi = 0.31830988618379067153776752675_dp !< 1/pi
  real(dp), public, parameter :: M_2_pi = 0.63661977236758134307553505349_dp !< 2/pi
  real(dp), public, parameter :: M_lnpi = 1.14472988584940017414342735135_dp !< ln(pi)
  real(dp), public, parameter :: M_ln2 = 0.693147180559945309417232121458176568_dp !< ln(2)
  real(dp), public, parameter :: M_e = 2.71828182845904523536028747135_dp !< e
  real(dp), public, parameter :: M_log2e = 1.442695040888963407359924681001892137_dp !< log_2 (e)
  real(dp), public, parameter :: M_log10e = 0.43429448190325182765112891892_dp !< log_10 (e)

  ! Include physical constants
  include "codata.inc"

  real(dp), public, parameter :: deg2rad = 0.017453292519943295_dp !< pi/180
  real(dp), public, parameter :: rad2deg = 57.295779513082320876654618_dp !< 180/pi

contains

  ! JF: Check if this is a good solution or there are better
  !> is_inf Checks if the argument is plus or minus infinite.
  !!
  function is_inf(x) result(y)
    implicit none
    integer :: y !< Returns 0 if is finite, -1 if is minus-infinite and +1 if is plus-infinite
    real(dp), intent(IN) :: x !< Variable to test

    if (abs(x) < HUGE(x)) then
      y = 0
    else if (x < 0) then
      y = -1
    else
      y = 1
    end if
  end function is_inf

  ! ! ----------------------- Timer functions -----------------------
  !> Reset timer
  subroutine reset_timer(T)
    implicit none
    type(timer), intent(OUT) :: T !<
    T%start_date = 0
    T%stop_date = 0
    T%date_elapsed = 0
    T%start_cputime = Zero
    T%stop_cputime = Zero
    T%elapsed = Zero
  end subroutine reset_timer

  !> Assign an integer value to a timer
  !!
  subroutine assign_int(T, val)
    implicit none
    class(timer), intent(OUT) :: T !< Time
    integer, intent(IN) :: val !< integer time
    T%start_date = 0
    T%stop_date = 0
    T%date_elapsed = val
    T%start_cputime = Zero
    T%stop_cputime = Zero
    T%elapsed = real(val, kind=dp)
  end subroutine assign_int

  !> Divide T by an integer value
  function divide_T_int(T1, val) result(T)
    implicit none
    class(timer), intent(IN) :: T1 !< Time
    integer, intent(IN) :: val  !< integer
    type(timer) :: T !< Resulting time
    integer :: i
    integer :: n
    T%start_date = T1%start_date
    T%stop_date = T1%stop_date
    T%date_elapsed = 0
    T%start_cputime = T1%start_cputime
    T%stop_cputime = T1%stop_cputime
    T%elapsed = T1%elapsed / real(val, kind=dp)
    do i = 1, 6
      n = mod(T1%date_elapsed(i), val)
      T%date_elapsed(i + 1) = (T1%date_elapsed(i + 1) + n * time_factors(i + 1)) / val
    end do

  end function divide_T_int

  !> Assign values from a timer to another
  !!
  subroutine assign(T, T1)
    implicit none
    class(timer), intent(OUT) :: T !< New time
    type(timer), intent(IN) :: T1 !< Original time
    T%start_date = T1%start_date
    T%stop_date = T1%stop_date
    T%date_elapsed = T1%date_elapsed
    T%start_cputime = T1%start_cputime
    T%stop_cputime = T1%stop_cputime
    T%elapsed = T1%elapsed
  end subroutine assign

  !> Record date. Notice that Time difference with UTC is the last
  function record_date() result(y)
    implicit none
    integer, dimension(7) :: y !< Date in format
    integer, dimension(8) :: tmp !< Date in format
    tmp = 0                      ! JF: It should not be necessary
    call date_and_time(values=tmp)
    y(1:3) = tmp(1:3)
    y(4:7) = tmp(5:8)
  end function record_date

  !< Start the timer
  !!
  subroutine start_timer(T)
    implicit none
    class(timer), intent(out) :: T !< Timer
    ! Set stop time and elapsed time to zero
    T%elapsed = Zero
    T%stop_cputime = Zero
    T%stop_date = 0
    T%date_elapsed = 0
    ! Set the start time
    T%start_date = record_date()
    call cpu_time(T%start_cputime)
  end subroutine start_timer

  !> stop the timer and calculate time-differences
  subroutine stop_timer(T)
    implicit none
    class(timer), intent(inout) :: T !< Timer

    call cpu_time(T%stop_cputime)
    T%stop_date = record_date()
    call diff_timer(T)

  end subroutine stop_timer

  !> add_timers Adds two time intervals
  !!
  function add_timers(T1, T2) result(T)
    implicit none
    class(timer), intent(IN) :: T1 !< first time interval
    class(timer), intent(IN) :: T2 !< second time interval
    type(timer) :: T !< Resulting time
    integer :: i

    T%elapsed = T1%elapsed + T2%elapsed
    T%date_elapsed = 0
    T%date_elapsed(4:) = T1%date_elapsed(4:) + T2%date_elapsed(4:)
    ! Convert to the higher unit (e.g: 61 s -> 1m 1s)
    do i = 7, 4, -1
      if (T%date_elapsed(i) > time_factors(i)) then
        T%date_elapsed(i) = T%date_elapsed(i) - time_factors(i)
        T%date_elapsed(i - 1) = T%date_elapsed(i - 1) + 1
      end if
    end do
  end function add_timers

  !> diff_timer
  !!
  !! Examples:
  !!
  subroutine diff_timer(T)
    implicit none
    type(timer), intent(INOUT) :: T !<
    integer :: i

    T%elapsed = T%stop_cputime - T%start_cputime ! in seconds
    T%date_elapsed = T%stop_date - T%start_date

    do i = 7, 4, -1
      if (T%date_elapsed(i) < 0) then
        T%date_elapsed(i) = T%date_elapsed(i) + time_factors(i)
        T%date_elapsed(i - 1) = T%date_elapsed(i - 1) - 1
      end if
    end do
  end subroutine diff_timer

  !> print_elapsed time to unit or stdout
  !!
  subroutine print_elapsed(T, unit)
    implicit none
    class(timer), intent(in) :: T !< Timer
    integer, optional :: unit     !< unit to output
    integer :: unit_

    unit_ = stdout; IF (present(unit)) unit_ = unit
    write (unit_, '(A,f0.2,A)') 'cpu time:   ', T%elapsed, 's'
    write (unit_, '(A)') 'Total time: '//stamp_date_diff(T%date_elapsed)
  end subroutine print_elapsed

  !> stamp_date Creates an string from a date
  !!
  !! @note that date format differs from the output of intrinsic date_and_time
  function stamp_date_diff(dd) result(output)
    implicit none
    integer, dimension(7) :: dd !> date-format differs from usual
    character(len=:), allocatable :: output
    character(len=4) :: tmp
    !
    integer :: i, n
    character(len=3), parameter, dimension(7) :: timeunit = [' Y,', ' M,', ' d,', ' h:', ' m:', ' s ', 'ms ']

    do n = 1, 7     ! Get first non-zero value
      IF (dd(n) /= 0) exit
    end do
    output = ''
    do i = n, 7
      write (tmp, '(I4)') dd(i)
      output = output//trim(adjustl(tmp))//timeunit(i)
    end do
  end function stamp_date_diff
  ! ------------------End of Timer functions -----------------------

  !> Print a message, and optionally stop the program
  !!
  !! If errcode > 0 stop the program
  !! Examples of use:
  !! ---------------
  !! ```
  !!
  !! call print_msg("Invalid argument!", errcode=1) ! will stop the program with error code 1.
  !! call print_msg("Invalid argument","my_sub",1)  ! will stop the program with error code 1.
  !! call print_msg("Fishy values of argument","my_sub", errcode=3) ! will stop the program with error code 3.
  !! call print_msg("Fishy values of argument", errcode=-3) ! will **NOT** stop the program. Show error code 3.
  !! call print_msg("Fishy values of argument","my_sub", errcode=0) ! will **NOT** stop the program.
  !! call print_msg("Fishy values of argument","my_sub", errcode=0, unit=13) ! will **NOT** stop the program.
  !!
  !! ```
  !! The first three examples print the given message and stop the program.
  !! The next two will keep running after printing the message to stderr.
  !! The last will keep running after printing the message to a previously open file.
  subroutine print_msg(msg, sub, errcode, unit)
    implicit none
    character(len=*), intent(IN) :: msg     !< Message to print on stderr
    character(len=*), optional, intent(IN) :: sub !< Routine name. Default = None
    integer, optional, intent(IN) :: errcode !< Error code. Default = 0
    integer, optional, intent(IN) :: unit !< Unit to write. Default = stderr
    integer :: errcode_, unit_
    character(len=2) :: tmp

    errcode_ = 0; IF (present(errcode)) errcode_ = errcode
    write (tmp, '(I2.2)') abs(errcode_)

    unit_ = stderr; IF (present(unit)) unit_ = unit

    if (present(sub)) then
      write (unit_, '(A)') 'Code: '//tmp//' in sub '//sub//' -> '//msg
    else
      write (unit_, '(A)') 'Code: '//tmp//' -> '//msg
    end if
    if (errcode_ > 0) stop 1
  end subroutine print_msg

end module basic

! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:
