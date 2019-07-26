    integer :: up, lo, mid
    integer :: j

    lo = 1
    up = size(x)
    if (x(lo) - elem > Delta) then ! elem below the array
      n = lo - 1                   ! Outside the array!!
      return
    end if

    if (elem - x(up) > Delta) then ! elem above the array
      n = up
      return
    end if

    ! If there are only a few elements just use a simple lookup and return -----
    if (up < minval_bisection) then
      do j = lo, up
        IF (x(j) - elem > Delta) exit
      end do
      n = j - 1
      return
    end if

    ! Else, we will use bisection ----------------------------------------------
    ! Instead of starting on (1, up//2 , up) we start with a linear guess
    n = int(((elem - x(lo)) / (x(up) - x(lo))) * up)

    IF (n < 1) n = 1

    ! Case where elem is higher or equal to the last element. We just return
    if (n > up) then
      n = up
      IF (elem >= x(n) - Delta) return
    end if

    ! The first two actions cover the case that by chance we found it
    If ((x(n) == elem) .or. ((elem > x(n)) .and. (elem < x(n + 1)))) then
      return                    ! Found the index
    else if (x(n + 1) == elem) then
      n = n + 1
      return
    else if (x(n) < elem) then
      lo = n
    else
      up = n
    end if

    do while ((up - lo > 1))
      mid = (up + lo) / 2 ! integer arithmetic
      if (elem >= x(mid)) then
        lo = mid
      else
        up = mid
      endif
    enddo
    n = lo
