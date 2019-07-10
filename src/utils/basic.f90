!> This module will provide some basic convenience routines
!! It should have:
!!  - types (real, double precision, etc) (?)
!!  - Basic error handling (defining output for warning and errors)
!!    Maybe look into https://docs.python.org/3/library/warnings.html
!!  - Some basic file handling (some os.path functionality) (?)
!!  - Some sys functionality (stdin, stdout, stderr)
!!  - Timers, timestamp (?)
!!

module basic
  use, intrinsic :: iso_fortran_env, only: stdin => input_unit, &
    stdout => output_unit, &
    stderr => error_unit
  use strings, only: str

  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.d0)

  !> Simple timer. Holds start-time, stop-time, and time-difference
  type timer
    integer, dimension(8) :: start_date !< Start date
    integer, dimension(8) :: stop_date !< Stop date
    integer, dimension(8) :: date_elapsed !< Date-difference between start and stop
    real(dp), private :: start_cputime    !< Start cputime
    real(dp), private :: stop_cputime     !< Stop cputime
    real(dp), public :: elapsed !< Elapsed time (in seconds)
  contains

    procedure, public :: start_timer
    procedure, public :: stop_timer
    procedure, public :: print_elapsed

  end type timer

  private
  public :: sp, dp, timer
  public :: stdout, stdin, stderr
  public :: start_timer, stop_timer, print_elapsed
  public :: print_msg

  ! Common constants
  real(dp), public, parameter :: Zero = 0._dp
  complex(dp), public, parameter :: C0 = (0._dp, 0._dp)
  complex(dp), public, parameter :: C1 = (1._dp, 0._dp)
  complex(dp), public, parameter :: CI = (0._dp, 1._dp)

  ! Basic mathematical constants (mostly from gsl)
  real(dp), public, parameter :: M_pi = 3.14159265358979323846264338328_dp      !<\< \f$ \pi \f$
  real(dp), public, parameter :: M_dpi = 6.28318530717958647692528676656_dp      !< 2*pi
  real(dp), public, parameter :: M_pi_2 = 1.57079632679489661923132169164_dp      !< pi/2
  real(dp), public, parameter :: M_pi_4 = 0.78539816339744830966156608458_dp      !< pi/4
  real(dp), public, parameter :: M_sqrtpi = 1.77245385090551602729816748334_dp      !< sqrt(pi)
  real(dp), public, parameter :: M_2_sqrtpi = 1.12837916709551257389615890312_dp      !< 2/sqrt(pi)
  real(dp), public, parameter :: M_1_pi = 0.31830988618379067153776752675_dp      !< 1/pi
  real(dp), public, parameter :: M_2_pi = 0.63661977236758134307553505349_dp      !< 2/pi
  real(dp), public, parameter :: M_lnpi = 1.14472988584940017414342735135_dp      !< ln(pi)
  real(dp), public, parameter :: M_ln2 = 0.693147180559945309417232121458176568_dp !< ln(2)
  real(dp), public, parameter :: M_e = 2.71828182845904523536028747135_dp !< e
  real(dp), public, parameter :: M_log2e = 1.442695040888963407359924681001892137_dp !< log_2 (e)
  real(dp), public, parameter :: M_log10e = 0.43429448190325182765112891892_dp !< log_10 (e)

  ! Include physical constants
  include "codata.inc"

  real(dp), public, parameter :: deg2rad = 0.017453292519943295_dp !< pi/180
  real(dp), public, parameter :: rad2deg = 57.295779513082320876654618_dp !< 180/pi
contains

  ! ----------------------- Timer functions -----------------------
  !> Record date. Notice that Time difference with UTC is the last
  function record_date(T) result(y)
    implicit none
    class(timer), intent(in) :: T !< Timer
    integer, dimension(8) :: y !< Date in format
    integer :: d
    call cpu_time(y)
    d = y(4)
    y(4:7) = y(5:8)
    y(8) = d
  end function record_date

  !< Start the timer
  subroutine start_timer(T)
    implicit none
    class(timer), intent(in) :: T !< Timer
    integer, dimension(8) :: d
    T%elapsed = Zero
    T%stop_cputime = Zero
    T%stop_date = 0
    T%date_elapsed = 0
    T%start_date = record_date()
    call cpu_time(T%start_cputime)
  end subroutine start_timer

  !> stop_timer the timer and calculate time-differences
  subroutine stop_timer(T)
    implicit none
    class(timer), intent(in) :: T !< Timer
    integer, dimension(4), parameter :: factors = [24, 60, 60, 1000] ! Factors to convert to higher time unit

    call cpu_time(T%stop_cputime)
    T%stop_date = record_date()
    T%elapsed = T%stop_cputime - T%start_cputime ! in seconds
    T%date_elapsed = T%start_date - T%stop_date
    do i = 7, 4, -1
      if (date_elapsed(i) < 0) then
        T%date_elapsed(i) = T%date_elapsed(i) + factors(i)
        T%date_elapsed(i - 1) = T%date_elapsed(i - 1) - 1
      end if
    end do
  end subroutine stop_timer

  !> print_elapsed time to unit or stdout
  subroutine print_elapsed(T, unit)
    implicit none
    class(timer), intent(in) :: T !< Timer
    integer, optional :: unit     !< unit to output
    integer :: unit_

    unit_ = stdout; IF (present(unit)) unit_ = unit

    write (unit_, '(A,f14.2,A)') 'cpu time: ', elapsed, 's'
    write (unit_, '(A)') stamp_date_diff(T%date_elapsed)
  end subroutine print_elapsed

  !> stamp_date Creates an string from a date
  !!
  !! @note that date format differs from the output of intrinsic date_and_time
  function stamp_date_diff(dd) result(output)
    implicit none
    integer, dimension(8) :: dd !> date-format differs from usual
    character(len=:), allocatable :: output

    integer :: i, n
    character(len=3), parameter, dimension(7) :: timeunit = [' Y,', ' M,', ' d,', ' h:', ' m:', ' s.', 'ms ']
    character(len=:), allocatable :: stamp

    do n = 1, 7     ! Get first non-zero value
      IF (dd(n) /= 0) exit
    end do
    output = ''
    do i = n, 7
      output = output//dd(i)//timeunit(i)
    end do
  end function stamp_date_diff
  ! ------------------End of Timer functions -----------------------

  !> Print an error, and optionally stop the program
  !!
  !! Examples of use:
  !! ---------------
  !! ```
  !!
  !! call print_msg("Invalid argument!")
  !! call print_msg("Invalid argument","my_sub")
  !! call print_msg("Fishy values of argument","my_sub", errcode=0, unit=13)
  !!
  !! ```
  !! The first two examples print the message and stop the program.
  !! The third will keep running after printing the message to a previously open file
  subroutine print_msg(msg, sub, errcode, unit)
    character(len=*) :: msg !< Message to print on stderr
    character(len=*), optional, intent(in) :: sub !< Routine name. Default = None
    integer, optional, intent(in) :: errcode !< Error code. Default = 1
    integer, optional, intent(in) :: unit    !< Unit to write. Default = stderr
    integer, errcode_, unit_
    errcode_ = 1; 
    IF (present(errcode)) errcode_ = errcode
    unit_ = stderr; IF (present(unit)) unit_ = unit
    write (unit_, '(A)') 'fatal error in sub:'//sub//':'//msg
    if (errcode_ /= 0) stop errcode_
  end subroutine print_msg

end module basic

