!> @file   mtrandom.F90
!> @author Juan Fiol <juanfiol@gmail.com> (modifications, see real authors below)
!> @date   Wed Jul 15 15:36:49 2009
!>
!> @brief
!> @remark
!> @verbatim
!> From the Algorithmic Conjurings of Scott Robert Ladd comes...
!>  An implementation of the Mersenne Twister algorithm for generating
!>  pseudo-random sequences.
!>
!>  History
!>  -------
!>   2.0.0   4 January 2004
!>           Corrected erroneous unsigned bit manipulations
!>           Doubled resolution by using 64-bit math
!>           Added rng_rand64
!>
!>  ORIGINAL ALGORITHM COPYRIGHT
!>  ============================
!>  Copyright (C) 1997,2002 Makoto Matsumoto and Takuji Nishimura.
!>  Any feedback is very welcome. For any question, comments, see
!>  http://www.math.keio.ac.jp/matumoto/emt.html or email
!>  matumoto@math.keio.ac.jp
!>---------------------------------------------------------------------
!>  COPYRIGHT NOTICE, DISCLAIMER, and LICENSE:
!>
!>  This notice applies *only* to this specific expression of this
!>  algorithm, and does not imply ownership or invention of the
!>  implemented algorithm.
!>
!>  If you modify this file, you may insert additional notices
!>  immediately following this sentence.
!>
!>  Copyright 2001, 2002, 2004 Scott Robert Ladd.
!>  All rights reserved, except as noted herein.
!>
!>  This computer program source file is supplied "AS IS". Scott Robert
!>  Ladd (hereinafter referred to as "Author") disclaims all warranties,
!>  expressed or implied, including, without limitation, the warranties
!>  of merchantability and of fitness for any purpose. The Author
!>  assumes no liability for direct, indirect, incidental, special,
!>  exemplary, or consequential damages, which may result from the use
!>  of this software, even if advised of the possibility of such damage.
!>
!>  The Author hereby grants anyone permission to use, copy, modify, and
!>  distribute this source code, or portions hereof, for any purpose,
!>  without fee, subject to the following restrictions:
!>
!>      1. The origin of this source code must not be misrepresented.
!>
!>      2. Altered versions must be plainly marked as such and must not
!>         be misrepresented as being the original source.
!>
!>      3. This Copyright notice may not be removed or altered from any
!>         source or altered source distribution.
!>
!>  The Author specifically permits (without fee) and encourages the use
!>  of this source code for entertainment, education, or decoration. If
!>  you use this source code in a product, acknowledgment is not required
!>  but would be appreciated.
!>
!>  Acknowledgement:
!>      This license is based on the wonderful simple license that
!>      accompanies libpng.
!>-----------------------------------------------------------------------
!>  For more information on this software package, please visit
!>      http://www.coyotegulch.com
!>@endverbatim
#include "base.h"

!> @brief Random generator for uniform distributions
!> @remarks they try to have only one (common) interface
!>
!>@details
!> - @b random returns random numbers
!> @code
!> random(state, r)
!> @endcode
!>   - @a state is a @c random_state type variable
!>   - @a @b r is either a real scalar or a 1D, 2D or 3D array
!> - @b seed Seed the generator
!> @code
!> seed(key, state)
!> @endcode
!>   - @a key is, if present, either an integer scalar or 1D array
!>   - @a state is a @c random_state type variable
#ifndef ACCESS
#define ACCESS 'stream'
#endif
#ifndef FORM
#define FORM 'unformatted'
#endif
module mtrandom
  implicit none

  ! Constants
  integer(4), private, parameter :: Nmt = 624_4
  integer(4), private, parameter :: M = 397_4
  ! Convenience factor
  real(dp), private, parameter :: factor = 1._dp / 4294967295._dp

  !> rng_state holds the state of the generator
  type rng_state
    integer(4)                   :: mti = -1
    integer(I64), dimension(0:Nmt - 1) :: mt
  end type rng_state

  !> initializes the generator. Use it as\n
  !> <b> seed(key, state) </b>
  !> - @a key is, if present, either an integer scalar or 1D array
  !> - @a state is a @c random_state type variable
  interface seed
    module procedure rng_init_by_scalar
    module procedure rng_init_by_array
  end interface seed

  !> returns one or several numbers\n
  !> <b> random(state, r) </b>
  !> - @a state is a @a random_state type variable
  !> - @a r is either a real scalar or a 1D, 2D or 3D array
  interface random
    module procedure rng_rand_real_array0d
    module procedure rng_rand_real_array1d
    module procedure rng_rand_real_array2d
    module procedure rng_rand_real_array3d
  end interface random

  private   ! Everything is private unless explicitly made public
  public :: rng_state
  public :: seed, random
  public :: rng_rand64, rng_rand, rng_rand_range
  public :: rng_rand_real1, rng_rand_real2, rng_rand_real3
  public :: test_random

contains

!> Initializes the generator with "semilla"
!! @param semilla
!! @param state n
!!
!! If semilla is zero or is not present, the semilla is read from /dev/urandom or from system
!! clock
  subroutine rng_init_by_scalar(semilla, state)
    integer, optional, intent(in)  :: semilla
    type(rng_state), intent(out) :: state
    integer :: i
    if (.NOT. Present(semilla) .or. semilla == 0) then
      state%mt(0) = read_urandom()
      IF (state%mt(0) == 0) state%mt(0) = iclock()
    else
      state%mt(0) = semilla    ! save semilla
    end if
    ! Set the seed using values suggested by Matsumoto & Nishimura, using a generator by
    !   Knuth. See original source for details.
    do i = 1, Nmt - 1
      state%mt(i) = iand(4294967295_8, 1812433253_8 * ieor(state%mt(i - 1), ishft(state%mt(i - 1), -30_8)) + i)
    end do
    state%mti = Nmt
  end subroutine rng_init_by_scalar

!> Initialize with an array of seeds
!! @param init_key is the seed
!! @param state is the state of the random generator
  subroutine rng_init_by_array(init_key, state)
    type(rng_state), intent(out) :: state
    integer(4), dimension(:), intent(in) :: init_key

    ! working storage
    integer :: key_length
    integer :: i
    integer :: j
    integer :: k

    call rng_init_by_scalar(19650218_4, state)

    i = 1
    j = 0
    key_length = size(init_key)

    do k = max(Nmt, key_length), 0, -1
      state%mt(i) = ieor(state%mt(i), (ieor(state%mt(i - 1), ishft(state%mt(i - 1), -30_8) * 1664525_8))) + init_key(j) + j

      i = i + 1
      j = j + 1

      if (i >= Nmt) then
        state%mt(0) = state%mt(Nmt - 1)
        i = 1
      end if

      if (j >= key_length) j = 0
    end do

    do k = Nmt - 1, 0, -1
      state%mt(i) = ieor(state%mt(i), (ieor(state%mt(i - 1), ishft(state%mt(i - 1), -30_8) * 1566083941_8))) - i

      i = i + 1
      if (i >= Nmt) then
        state%mt(0) = state%mt(Nmt - 1)
        i = 1
      end if
    end do

    state%mt(0) = 1073741824_8 ! 0x40000000, assuring non-zero initial array
  end subroutine rng_init_by_array

  !>  Obtain the next 64-bit (orig said 32??) integer in the pseudo-random sequence
  !> @param state
  !> @return r a 64-bit integer
  !> @remarks variable names match those in original example
  function rng_rand64(state) result(r)
    type(rng_state), intent(inout) :: state
    integer(I64) :: r
    ! internal constants
    integer(I64), dimension(0:1), parameter :: mag01 = (/0_8, -1727483681_8/)
    ! Period parameters
    integer(I64), parameter :: UPPER_MASK = 2147483648_8
    integer(I64), parameter :: LOWER_MASK = 2147483647_8
    ! Tempering parameters
    integer(I64), parameter :: TEMPERING_B = -1658038656_8
    integer(I64), parameter :: TEMPERING_C = -272236544_8

    integer(4) :: kk

    ! Generate N words at a time
    if (state%mti >= Nmt .or. state%mti < 0) then
      ! The value -1 acts as a flag saying that the seed has not been set.
      IF (state%mti == -1) call rng_init_by_scalar(0_4, state)

      ! Fill the mt array
      do kk = 0, Nmt - M - 1
        r = ior(iand(state%mt(kk), UPPER_MASK), iand(state%mt(kk + 1), LOWER_MASK))
        state%mt(kk) = ieor(ieor(state%mt(kk + M), ishft(r, -1_8)), mag01(iand(r, 1_8)))
      end do

      do kk = Nmt - M, Nmt - 2
        r = ior(iand(state%mt(kk), UPPER_MASK), iand(state%mt(kk + 1), LOWER_MASK))
        state%mt(kk) = ieor(ieor(state%mt(kk + (M - Nmt)), ishft(r, -1_8)), mag01(iand(r, 1_8)))
      end do

      r = ior(iand(state%mt(Nmt - 1), UPPER_MASK), iand(state%mt(0), LOWER_MASK))
      state%mt(Nmt - 1) = ieor(ieor(state%mt(M - 1), ishft(r, -1)), mag01(iand(r, 1_8)))

      ! Start using the array from first element
      state%mti = 0
    end if

    ! Here is where we actually calculate the number with a series of transformations
    r = state%mt(state%mti)
    state%mti = state%mti + 1

    r = ieor(r, ishft(r, -11))
    r = iand(4294967295_8, ieor(r, iand(ishft(r, 7), TEMPERING_B)))
    r = iand(4294967295_8, ieor(r, iand(ishft(r, 15), TEMPERING_C)))
    r = ieor(r, ishft(r, -18))
  end function rng_rand64

  !--------------------------------------------------------------------------
  !> Obtain the next 32-bit integer in the pseudo-random sequence
  !> @param state
  !> @return r
  function rng_rand(state) result(r)
    type(rng_state), intent(inout) :: state
    integer(4) :: r
    integer(I64) :: x         !>!<  working storage
    x = rng_rand64(state)
    if (x > 2147483647_8) then
      r = x - 4294967296_8
    else
      r = x
    end if
  end function rng_rand

  !> Obtain a pseudorandom integer in the range [lo,hi]
  !> @param state
  !> @param lo minimum value possible on output
  !> @param hi minimum value possible on output
  !>
  !> @return one pseudorandom integer less or equal than @c lo and greater or equal than @c hi
  function rng_rand_range(state, lo, hi) result(r)
    type(rng_state), intent(inout) :: state
    integer, intent(in) :: lo
    integer, intent(in) :: hi
    integer(4) :: r
    ! Use real value to calculate range
    r = lo + floor((hi - lo + 1.0_8) * rng_rand_real2(state))
  end function rng_rand_range

  !> Obtain a pseudorandom real number in the range [0,1],
  !> @param state
  !> @return a number greater than or equal to 0 and less than or equal to 1.
  function rng_rand_real1(state) result(r)
    type(rng_state), intent(inout) :: state
    real(dp) :: r
    !
    r = real(rng_rand64(state), dp) * factor
  end function rng_rand_real1

  !> Obtain a pseudorandom real number in the range [0,1)
  !> @param state
  !>
  !> @return a number greater than 0 and less than or equal to 1.
  function rng_rand_real2(state) result(r)
    type(rng_state), intent(inout) :: state
    real(dp) :: r
    r = real(rng_rand64(state), dp) * factor
  end function rng_rand_real2

  !> Obtain a pseudorandom real number in the range (0,1)
  !> @param  state
  !> @return a number greater than 0 and less than 1.
  function rng_rand_real3(state) result(r)
    type(rng_state), intent(inout) :: state
    real(dp) :: r
    r = (real(rng_rand64(state), dp) + 0.5_dp) * factor
  end function rng_rand_real3

  !> Obtain a pseudorandom real number in the range [0,1],
  !> \param state
  !> \param r is the scalar output
  !> \return a number greater than or equal to 0 and less than or equal to 1.
  subroutine rng_rand_real_array0d(state, r)
    type(rng_state), intent(INOUT) :: state
    real(dp), intent(INOUT) :: r
    !
    r = rng_rand_real1(state)
  end subroutine rng_rand_real_array0d

  !> Obtain an array (1,...,dim) of pseudorandom real numbers in the range [0,1]
  !> \param state
  !> \param r is the 1D-array output
  !> \return an array of numbers greater than or equal to 0 and less than or equal to 1.
  subroutine rng_rand_real_array1d(state, r)
    type(rng_state), intent(INOUT) :: state
    real(dp), dimension(:), intent(INOUT) :: r
    integer :: i1
    ! Local constant; precalculated to avoid division below
    do i1 = 1, size(r)
      r(i1) = real(rng_rand64(state), dp) * factor
    end do
  end subroutine rng_rand_real_array1d

  !> Obtain a 2D-array (1,...,dim1)(1,...,dim2) of pseudorandom real numbers in the range
  !> [0,1]
  !> \param state
  !> \param r is the 2d-array output
  !> \return a 2d-array of numbers greater than or equal to 0 and less than or equal to 1
  subroutine rng_rand_real_array2d(state, r)
    type(rng_state), intent(inout) :: state
    real(dp), dimension(:, :), intent(inout) :: r
    integer :: i1
    do i1 = 1, size(r, 2)
      call rng_rand_real_array1d(state, r(:, i1))
    end do
  end subroutine rng_rand_real_array2d

  !> Obtain a 3D-array of pseudorandom real numbers in the range [0,1]
  !> \param state
  !> \param r is the 3D-array output
  !> \return a 3d-array of numbers greater than or equal to 0 and less than or equal to 1
  subroutine rng_rand_real_array3d(state, r)
    type(rng_state), intent(inout) :: state
    real(dp), dimension(:, :, :), intent(inout) :: r
    integer :: i1
    do i1 = 1, size(r, 3)
      call rng_rand_real_array2d(state, r(:, :, i1))
    end do
  end subroutine rng_rand_real_array3d

  !> Set a positive integer from the clock values
  !> @return a non-zero positive integer
  integer function iclock()
    integer, dimension(dp) :: tiempo
    call date_and_time(values=tiempo)
    iclock = abs(tiempo(3) + 10 * tiempo(4) + 100 * tiempo(5) + 1000 * tiempo(6) + 10000 * tiempo(7) + 100000 * tiempo(dp))
  end function iclock

  !> gets a random number from /dev/urandom or /dev/random
  !> @remarks Works with ifort, to use with gfortran, first line should read
  !> @code open(89,file='/dev/urandom',access='stream', form='unformatted') @endcode
  function read_urandom()
    integer(4) :: read_urandom
    open (89, file='/dev/random', action='read', access=ACCESS, form=FORM, err=1)
    read (89) read_urandom
    close (89)
    return
1   read_urandom = 0
  end function read_urandom

end module mtrandom
