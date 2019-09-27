subroutine fpbspl(t, n, k, x, l, h)
  !  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
  !  degree k at t(l) <= x < t(l+1) using the stable recurrence
  !  relation of de boor and cox.
  !  Travis Oliphant  2007
  !    changed so that weighting of 0 is used when knots with
  !      multiplicity are present.
  !    Also, notice that l+k <= n and 1 <= l+1-k
  !      or else the routine will be accessing memory outside t
  !      Thus it is imperative that that k <= l <= n-k but this
  !      is not checked.
  !  ..
  !  ..scalar arguments..
  real(8), intent(IN) :: x
  integer, intent(IN) :: n, k, l
  !  ..array arguments..
  real(8), intent(IN), dimension(n):: t
  real(8), intent(OUT), dimension(20) :: h
  !  ..local scalars..
  real(8) :: f
  integer :: i, j, li, lj
  !  ..local arrays..
  real(8) :: hh(19)
  !  ..

  h(1) = 1._8
  do j = 1, k
    hh(:j) = h(:j)
    h(1) = 0.0_8
    do i = 1, j
      li = l + i
      lj = li - j
      if (t(li) == t(lj)) then
        h(i + 1) = 0.0_8
      else
        f = hh(i) / (t(li) - t(lj))
        h(i) = h(i) + f * (t(li) - x)
        h(i + 1) = f * (x - t(lj))
      end if
    end do
  end do

end subroutine fpbspl
