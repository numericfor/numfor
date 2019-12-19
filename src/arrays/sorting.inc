! From https://github.com/jacobwilliams/fortran-search-and-sort.git
  recursive subroutine quicksort(ilow, ihigh)

    !! Sort the array (ascending order).

    implicit none

    integer, intent(in) :: ilow
    integer, intent(in) :: ihigh

    integer :: ipivot !! pivot element
    integer :: i      !! counter
    integer :: j      !! counter

    IF (ihigh <= ilow) return   ! Nothing to do
    if (ihigh - ilow <= maxsize_insertion) then
      ! do insertion sort:
      insert: do i = ilow + 1, ihigh
        do j = i, ilow + 1, -1
          if (vec(j) < vec(j - 1)) then
            call swap(vec(j), vec(j - 1))
          else
            exit                ! j-loop
          end if
        end do
      end do insert
    else      ! do normal quicksort
      call partition(ilow, ihigh, ipivot)
      call quicksort(ilow, ipivot - 1)
      call quicksort(ipivot + 1, ihigh)
    end if

  end subroutine quicksort

  subroutine partition(ilow, ihigh, ipivot)

    !! Partition the array

    implicit none

    integer, intent(in)  :: ilow
    integer, intent(in)  :: ihigh
    integer, intent(out) :: ipivot

    integer :: i, ip

    call swap(vec(ilow), vec((ilow + ihigh) / 2))
    ip = ilow
    do i = ilow + 1, ihigh
      if (vec(i) < vec(ilow)) then
        ip = ip + 1
        call swap(vec(ip), vec(i))
      end if
    end do
    call swap(vec(ilow), vec(ip))
    ipivot = ip

  end subroutine partition

