module histograms

  USE basic, only: dp, Zero, Tiny, print_msg
  USE grids, only: linspace, searchsorted, mean, std

  !> type histogram holds the data from an histogram
  type, public :: histog
    integer :: n = 0            !< Total number of values parsed
    real(dp), dimension(:), pointer :: hist !< Values of the histogram
    real(dp), dimension(:), pointer :: bin_edges !< Array (2xlen(hist)) the edges
  end type histog

  public :: histogram

contains
  !> Calculates an histogram of an array of data
  !!
  function histogram(a, Nbins, bins, range, density) result(h)
    USE basic, only: timer
    type(histog) :: h !< The histogram construct
    real(dp), dimension(:), target, intent(IN) :: a !< Input data
    integer, optional, intent(IN) :: Nbins !< Number of equal-width bins to use
    real(dp), dimension(:), optional, intent(IN) :: bins  !<
    real(dp), dimension(2), optional, intent(IN) :: range !<
    logical, optional, intent(IN) :: density !<
    real(dp), dimension(:), pointer :: p
    real(dp), dimension(:), pointer :: fb, lb
    logical :: d_
    integer :: Nbins_
    integer :: nn
    integer :: i, j, k
    integer :: status
    integer, parameter :: BLOCK = 65536
    ! integer, parameter :: BLOCK = huge(1)
    real(dp), parameter :: weight = 1.0_dp
    type(timer) :: T

    d_ = .FALSE.; IF (Present(density)) d_ = density

    h%n = size(a)
    call get_bins(a, Nbins, bins, range, h)
    Nbins_ = size(h%bin_edges) - 1

    allocate (h%hist(Nbins_), stat=status)
    IF (status /= 0) call print_msg('Error allocating hist', 'histogram')

    fb => h%bin_edges(1:)
    lb => h%bin_edges(2:)
    h%hist = Zero

    call T%start()
    nn = 0
    ! The use of the BLOCK does not appear to be very important for speed
    do concurrent(k=1:size(a):BLOCK)
      p => a(k:)
      N = min(h%n - k + 1, BLOCK) ! Number of elements in block

      do concurrent(j=1:N)
        find_bin: do i = 1, Nbins_
          if ((p(j) >= fb(i)) .and. (p(j) < lb(i))) then
            h%hist(i) = h%hist(i) + weight
            nn = nn + 1
            exit find_bin
          end if
        end do find_bin
      end do

    end do
    h%n = nn

    IF (d_) h%hist = h%hist / (lb - fb) / h%n
    call T%stop()
    call T%show()

  end function histogram

  ! TODO: We need to test which automatic method works best
  subroutine get_bins(a, Nbins, bins, range, h)
    implicit none

    real(dp), dimension(:), intent(IN) :: a !<
    integer, optional, intent(IN) :: Nbins !<
    real(dp), dimension(:), optional, intent(IN) :: bins !<
    real(dp), dimension(2), optional, intent(IN) :: range !<
    real(dp) :: first_edge, last_edge, width
    integer :: Nbins_
    integer :: status
    integer :: i, j
    type(histog), intent(INOUT) :: h

    if (Present(range)) then     ! Clip to range if present
      IF (range(1) > range(2)) call print_msg('Max must be larger than min', 'histogram')
      first_edge = range(1)
      last_edge = range(2)
    else                        ! Otherwise use extreme values
      first_edge = minval(a)
      last_edge = maxval(a)
    end if
    ! finally, correct the limits if they are equal
    if (first_edge == last_edge) then
      first_edge = first_edge - 0.5_dp
      last_edge = last_edge + 0.5_dp
    end if

    if (Present(bins)) then
      ! Check that they are ordered
      IF (any((bins(:size(bins) - 1) > bins(1:)))) &
        & call print_msg('bins must increase monotonically', 'histogram')

      ! We clip bins if needed
      i = searchsorted(bins, first_edge)
      j = searchsorted(bins, last_edge)
      IF (bins(i) - first_edge < -Tiny) i = i + 1

      allocate (h%bin_edges(j - i + 1), stat=status)
      IF (status /= 0) call print_msg('Error allocation bin_edges', 'histogram')

      h%bin_edges = bins(i:j)
      print *, 'first', range(1), bins(i)
      print *, 'last', range(2), bins(j)
      return

    else if (Present(Nbins)) then
      Nbins_ = Nbins + 1
    else
      ! We try to guess the "optimal bin edges"
      ! TODO: Here we should check wich method is best
      ! width = width_sturges(a)
      ! width = width_rice(a)
      width = width_doane(a)
      Nbins_ = ceiling(abs((last_edge - first_edge)) / width)
    end if

    allocate (h%bin_edges(Nbins_))
    h%bin_edges = linspace(first_edge, last_edge, Nbins_, endpoint=.True.)
  end subroutine get_bins

!> width_sturges
!!
  function width_sturges(x) result(w)
    implicit none
    real(8) :: w !<
    real(8), dimension(:), intent(IN) :: x !< Array to analyze
    w = (maxval(x) - minval(x)) * log10(2._dp) / log10(size(x) + 1._dp)
  end function width_sturges!> width_sturges

!>
!!
  function width_rice(x) result(w)
    implicit none
    real(8) :: w !<
    real(8), dimension(:), intent(IN) :: x !< Array to analyze
    w = (maxval(x) - minval(x)) / (2._dp * size(x)**(1._dp / 3._dp))
  end function width_rice

!> Doane's histogram bin estimator.
!!
!! Improved version of Sturges' formula which works better for non-normal data.
!! @see stats.stackexchange.com/questions/55134/doanes-formula-for-histogram-binning
  function width_doane(x) result(w)
    implicit none
    real(dp) :: w !< width
    real(dp), dimension(:), intent(IN) :: x !< Array to analyze
    real(dp), dimension(size(x)) :: temp
    real(dp) :: g1
    real(dp) :: sigma
    real(dp) :: log10_2
    integer :: N
    w = (maxval(x) - minval(x))
    N = size(x)
    IF (N <= 2) return
    log10_2 = log10(2._dp)
    sigma = std(x)
    if (sigma > Zero) then
      temp = x - mean(x)
      temp = (sigma / temp)**3
      g1 = abs(mean(temp)) / sqrt(6._dp * (N - 2) / ((N + 1._dp) * (N + 3._dp)))
      w = w * log10_2 / (log10_2 + log10(real(N, kind=dp)) + log10(1._dp + g1))
    end if
  end function width_doane

end module histograms
