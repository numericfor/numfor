! Template to create routines that have two scalar arguments and one of variable shape
#define RNG JOIN(ROUTINE_NAME,0d)
  !> Fills a scalar
  subroutine RNG(A1, A2, x)
    implicit none
    real(dp), intent(IN) :: A1
    real(dp), intent(IN) :: A2
    real(dp), intent(OUT) :: x
    x = RNG_NAME(A1, A2)
  end subroutine RNG
#undef RNG

#define RNG1 JOIN(ROUTINE_NAME,0d)
#define RNG JOIN(ROUTINE_NAME,1d)
  !> Fills a 1d array
  subroutine RNG(A1, A2, x)
    implicit none
    real(dp), intent(IN) :: A1
    real(dp), intent(IN) :: A2
    real(dp), dimension(:), intent(OUT) :: x !< object to fill
    integer :: i
    do i = 1, size(x, 1)
      call RNG1(A1, A2, x(i))
    end do
  end subroutine RNG
#undef RNG
#undef RNG1

#define RNG1 JOIN(ROUTINE_NAME,1d)
#define RNG JOIN(ROUTINE_NAME,2d)
  !> Fills a 2d array
  subroutine RNG(A1, A2, x)
    implicit none
    real(dp), intent(IN) :: A1
    real(dp), intent(IN) :: A2
    real(dp), dimension(:, :), intent(OUT) :: x !< object to fill
    integer :: i
    do i = 1, size(x, 2)
      call RNG1(A1, A2, x(:, i))
    end do
  end subroutine RNG
#undef RNG
#undef RNG1

#define RNG1 JOIN(ROUTINE_NAME,2d)
#define RNG JOIN(ROUTINE_NAME,3d)
  !> Fills a 3d array
  subroutine RNG(A1, A2, x)
    implicit none
    real(dp), intent(IN) :: A1
    real(dp), intent(IN) :: A2
    real(dp), dimension(:, :, :), intent(OUT) :: x !<object to fill
    integer :: i
    do i = 1, size(x, 3)
      call RNG1(A1, A2, x(:, :, i))
    end do
  end subroutine RNG
#undef RNG
#undef RNG1

#define RNG1 JOIN(ROUTINE_NAME,3d)
#define RNG JOIN(ROUTINE_NAME,4d)
  !> Fills a 4d array
  subroutine RNG(A1, A2, x)
    implicit none
    real(dp), intent(IN) :: A1
    real(dp), intent(IN) :: A2
    real(dp), dimension(:, :, :, :), intent(OUT) :: x !< object to fill
    integer :: i
    do i = 1, size(x, 4)
      call RNG1(A1, A2, x(:, :, :, i))
    end do
  end subroutine RNG
#undef RNG
#undef RNG1

#define RNG1 JOIN(ROUTINE_NAME,4d)
#define RNG JOIN(ROUTINE_NAME,5d)
  !> Fills a 5d array
  subroutine RNG(A1, A2, x)
    implicit none
    real(dp), intent(IN) :: A1
    real(dp), intent(IN) :: A2
    real(dp), dimension(:, :, :, :, :), intent(OUT) :: x !< object to fill
    integer :: i
    do i = 1, size(x, 5)
      call RNG1(A1, A2, x(:, :, :, :, i))
    end do
  end subroutine RNG
#undef RNG
#undef RNG1

! Local variables:
! eval: (add-hook 'before-save-hook 'time-stamp)
! time-stamp-start: "date[ ]+\\\\?[\"]+"
! time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S"
! End:
