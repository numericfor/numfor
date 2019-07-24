!> The module `histograms` provides infrastructure for the calculation of histograms. i.e: a count of frequency
!!
!!
module histograms

  USE basic, only: dp, Zero, Small, print_msg
  USE grids, only: linspace, searchsorted, mean, std

  !> type histogram holds the data from an histogram
  type, public :: histog
    integer :: n = 0            !< Total number of values parsed
    real(dp), dimension(:), pointer :: hist !< Values of the histogram
    real(dp), dimension(:), pointer :: bin_edges !< Array of size (len(hist)+1)
  end type histog

  public :: histogram

contains
  !> Calculates an histogram of an array of data
  !!
  !! @note It counts the number of occurrences in the interval [bins(i), bins(i+1) )
  !! which is closed on the left and open on the right (does not include the upper limit)
  !! except for the last one which includes both limits
  !!
  !! Example:
  !! ```
  !! h = histogram([1, 2, 1], bins=[0, 1, 2, 3])
  !! h%bin_edges = [0,1,2,3]
  !! h%hist = [1, 2, 1]
  !! ```
  !! Or
  !! ```
  !! hist = histogram(, bins=bins, density=.True.)
  !! hist = histogram(a, bins=bins, density=.False.)

  function histogram(a, Nbins, bins, range, density) result(h)
    type(histog) :: h !< The histogram construct
    real(dp), dimension(:), target, intent(IN) :: a !< Input data
    integer, optional, intent(IN) :: Nbins !< Number of equal-width bins to use
    real(dp), dimension(:), optional, intent(IN) :: bins  !< defines a monotonically increasing array of bin edges,
    !< including the rightmost edge, allowing for non-uniform bin widths.
    real(dp), dimension(2), optional, intent(IN) :: range !< The lower and upper range of the bins. If not provided uses `min` and `max`
    logical, optional, intent(IN) :: density !< If ``True``, the result is the value of the
    !! probability *density* function at the bin, normalized
    !! such that the *integral* over the range is 1.
    real(dp), dimension(:), pointer :: p
    real(dp), dimension(:), pointer :: fb, lb
    logical :: d_
    integer :: Nbins_
    integer :: nn
    integer :: i, j, k
    integer :: status
    integer, parameter :: BLOCK = 65536
    real(dp), parameter :: weight = 1.0_dp
    real(dp) :: first_edge, last_edge

    d_ = .FALSE.; IF (present(density)) d_ = density

    h%n = size(a)
    call get_bins(a, Nbins, bins, range, h)
    Nbins_ = size(h%bin_edges) - 1

    ! Allocate memory for hist component
    allocate (h%hist(Nbins_), stat=status)
    IF (status /= 0) call print_msg('Error allocating hist', 'histogram')

    fb => h%bin_edges(1:Nbins_)
    lb => h%bin_edges(2:)
    first_edge = h%bin_edges(1); last_edge = h%bin_edges(Nbins_ + 1)
    h%hist = Zero

    nn = 0
    ! The use of the BLOCK does not appear to be very important for speed
    ! blocks: do concurrent(k=1:size(a):BLOCK)
    blocks: do k = 1, size(a), BLOCK
      p => a(k:)
      N = min(h%n - k + 1, BLOCK) ! Number of elements in block

      ! For a small number of bins it seems more efficient to just look sequencially
      if (Nbins_ <= 100) then
        ! do concurrent(j=1:N)
        do j = 1, N
          find_bin: do i = 1, Nbins_
            if ((p(j) >= fb(i)) .and. (p(j) < lb(i))) then
              h%hist(i) = h%hist(i) + weight
              nn = nn + 1
              exit find_bin
            end if
          end do find_bin

        end do
      else
        do concurrent(j=1:N)
          i = searchsorted(h%bin_edges, p(j))
          if ((i > 0) .and. (i <= Nbins_)) then
            h%hist(i) = h%hist(i) + weight
            nn = nn + 1
          end if
        end do
      end if

    end do blocks

    h%n = nn

    IF (d_) h%hist = h%hist / (lb - fb) / h%n

  end function histogram

  ! TODO: We need to test which automatic method works best
  ! Discussion. In order to make the last bin a closed-interval on both sides:
  ! [bins(N-2), bins(N-1)] which included both limits
  ! as opposed to [bins(N-2), bins(N-1)], which is closed on the left and open on the right
  ! we add a small quantity to the right edge (2*epsilon)
  ! where epsilon is the smallest difference between two real(dp) numbers
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

    if (present(range)) then     ! Clip to range if present
      IF (range(1) > range(2)) call print_msg('Max must be larger than min', 'histogram')
      first_edge = range(1)
      last_edge = range(2)
    else                        ! Otherwise use extreme values
      first_edge = minval(a)
      last_edge = maxval(a)
    end if

    ! Then, we correct the limits if they are equal
    if (first_edge == last_edge) then
      first_edge = first_edge - 0.5_dp
      last_edge = last_edge + 0.5_dp
    end if

    ! Finally we add a correction in order to compare real numbers
    first_edge = first_edge - Small
    last_edge = last_edge + Small

    if (Present(bins)) then
      ! Check that they are ordered
      IF (any((bins(:size(bins) - 1) > bins(2:)))) &
        & call print_msg('bins must increase monotonically', 'histogram')

      ! We clip bins if needed
      i = searchsorted(bins, first_edge)
      j = searchsorted(bins, last_edge)

      IF (i < 1) i = 1          ! Correct if inferior limit in range smaller than first edge

      allocate (h%bin_edges(i:j), stat=status)
      IF (status /= 0) call print_msg('Error allocation bin_edges', 'histogram')
      h%bin_edges = bins(i:j)
      ! We extend the last one in epsilon to make it inclusive
      h%bin_edges(j) = h%bin_edges(j) + 1.5_dp * epsilon(1._dp)
      return

    else if (present(Nbins)) then
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
    h%bin_edges(1) = h%bin_edges(1) - Small
    h%bin_edges(Nbins_) = h%bin_edges(Nbins_) + Small
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
