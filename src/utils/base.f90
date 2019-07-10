!> This module will provide some basic convenience routines
!! It should have:
!!  - types (real, double precision, etc) (?)
!!  - Basic error handling (defining output for warning and errors)
!!    Maybe look into https://docs.python.org/3/library/warnings.html
!!  - Some basic file handling (some os.path functionality) (?)
!!  - Some sys functionality (stdin, stdout, stderr)
module basic
  use, intrinsic :: iso_fortran_env, only: stdin => input_unit, &
    stdout => output_unit, &
    stderr => error_unit

  implicit none

  public :: stdout, stdin, stderr
contains

end module basic

