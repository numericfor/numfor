subroutine fprppo(nu, nv, if1, if2, cosi, ratio, c, f, ncoff)
  !  given the coefficients of a constrained bicubic spline, as determined
  !  in subroutine fppola, subroutine fprppo calculates the coefficients
  !  in the standard b-spline representation of bicubic splines.
  !  ..
  !  ..scalar arguments..
  real(8) :: ratio
  integer :: nu, nv, if1, if2, ncoff
  !  ..array arguments
  real(8) :: c(ncoff), f(ncoff), cosi(5, nv)
  !  ..local scalars..
  integer :: i, iopt, ii, j, k, l, nu4, nvv
  !  ..
  nu4 = nu - 4
  nvv = nv - 7
  iopt = if1 + 1

  f(:ncoff) = 0._8

  i = 0
  do l = 1, nu4
    ii = i
    if (l > iopt) then
      if (l /= nu4 .or. if2 == 0) then
        do k = 1, nvv
          i = i + 1
          j = j + 1
          f(i) = c(j)
        end do
      end if

    else if (l == 1) then

      do k = 1, nvv
        i = i + 1
        f(i) = c(1)
      end do

      j = 1

    else if (l == 2) then

      do k = 1, nvv
        i = i + 1
        f(i) = c(1) + c(2) * cosi(1, k) + c(3) * cosi(2, k)
      end do

      j = 3

    else if (l == 3) then

      do k = 1, nvv
        i = i + 1
        f(i) = c(1) + ratio * (c(2) * cosi(1, k) + c(3) * cosi(2, k)) +&
          &c(4) * cosi(3, k) + c(5) * cosi(4, k) + c(6) * cosi(5, k)
      end do

      j = 6

    end if

    do k = 1, 3
      ii = ii + 1
      i = i + 1
      f(i) = f(ii)
    end do

  end do

  c(:ncoff) = f(:ncoff)

end subroutine fprppo

