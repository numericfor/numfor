PRNM(trapz2d_func) (func, xrange, yrange, Nx, Ny) result(int)
NUMFOR_KINDF :: func
NUMFOR_KINDR :: int
real(8), dimension(2), intent(IN) :: xrange, yrange
integer, intent(IN) :: Nx, Ny
integer :: i, j
real(8) :: xx(Nx), yy(Ny)
real(8) :: hx, hy
hx = (xrange(2) - xrange(1)) / Nx
hy = (yrange(2) - yrange(1)) / Ny
xx = linspace(xrange(1), xrange(2), Nx)
yy = linspace(yrange(1), yrange(2), Ny)

int = &
  NUMFOR_EVAL_F(func, ([xrange(1), yrange(1)])) + &
  NUMFOR_EVAL_F(func, ([xrange(1), yrange(2)])) + &
  NUMFOR_EVAL_F(func, ([xrange(2), yrange(1)])) + &
  NUMFOR_EVAL_F(func, ([xrange(2), yrange(2)]))

do i = 2, Nx
  do j = 2, Ny
    int = int + 4._8 * NUMFOR_EVAL_F(func, ([xx(i), yy(j)]))
  enddo
enddo
do j = 2, Ny
  int = int + 2._8 * (NUMFOR_EVAL_F(func, ([xrange(1), yy(j)])) &
    & + NUMFOR_EVAL_F(func, ([xrange(2), yy(j)])))
enddo
do i = 2, Nx
  int = int + 2._8 * (NUMFOR_EVAL_F(func, ([xx(i), yrange(1)]))&
    & + NUMFOR_EVAL_F(func, ([xx(i), yrange(2)])))
enddo
int = int * hx * hy / 4._8
end PRNM(trapz2d_func)

PRNM(simps2d_func) (func, xrange, yrange, Nx, Ny) result(int)
NUMFOR_KINDF :: func
NUMFOR_KINDR :: int
real(8), dimension(2) :: xrange, yrange
integer :: Nx, Ny, i, j
real(8) :: xx(2 * Nx), yy(2 * Ny)
real(8) :: hx, hy

hx = xrange(2) - xrange(1)
hx = hx / Nx / 2
hy = yrange(2) - yrange(1)
hy = hy / Ny / 2
xx = linspace(xrange(1), xrange(2), 2 * Nx)
yy = linspace(yrange(1), yrange(2), 2 * Ny)
!
int = &
  NUMFOR_EVAL_F(func, ([xrange(1), yrange(1)])) +&
  & NUMFOR_EVAL_F(func, ([xrange(1), yrange(2)])) +&
  & NUMFOR_EVAL_F(func, ([xrange(2), yrange(1)])) +&
  & NUMFOR_EVAL_F(func, ([xrange(2), yrange(2)]))
!

do j = 1, Ny
  int = int + 4._8 * (NUMFOR_EVAL_F(func, ([xrange(1), yy(2 * j - 1)])) + &
    & NUMFOR_EVAL_F(func, ([xrange(2), yy(2 * j - 1)])))
enddo
do j = 1, Ny - 1
  int = int + 2._8 * (NUMFOR_EVAL_F(func, ([xrange(1), yy(2 * j)])) + NUMFOR_EVAL_F(func, ([xrange(2), yy(2 * j)])))
enddo
do i = 1, Nx
  int = int + 4._8 * (NUMFOR_EVAL_F(func, ([xx(2 * i - 1), yrange(1)])) + NUMFOR_EVAL_F(func, ([xx(2 * i - 1), yrange(2)])))
enddo
do i = 1, Nx - 1
  int = int + 2._8 * (NUMFOR_EVAL_F(func, ([xx(2 * i), yrange(1)])) + NUMFOR_EVAL_F(func, ([xx(2 * i), yrange(2)])))
enddo
!
do j = 1, Ny
  do i = 1, Nx
    int = int + 16._8 * NUMFOR_EVAL_F(func, ([xx(2 * i - 1), yy(2 * j - 1)]))
  enddo
enddo
do j = 1, Ny - 1
  do i = 1, Nx
    int = int + 8._8 * NUMFOR_EVAL_F(func, ([xx(2 * i - 1), yy(2 * j)]))
  enddo
enddo
!
do j = 1, Ny
  do i = 1, Nx - 1
    int = int + 8._8 * NUMFOR_EVAL_F(func, ([xx(2 * i), yy(2 * j - 1)]))
  enddo
enddo
do j = 1, Ny - 1
  do i = 1, Nx - 1
    int = int + 4._8 * NUMFOR_EVAL_F(func, ([xx(2 * i), yy(2 * j)]))
  enddo
enddo
int = int * hx * hy / 9.0_8
end PRNM(simps2d_func)

PRNM(trapz2d_func_eps) (func, xrange, yrange, eps, limit) result(int)
implicit none
NUMFOR_KINDF :: func
NUMFOR_KINDR :: int
NUMFOR_KINDR :: int0
#define NUMFOR_INTEG(f,x,y,N,M) trapz(f,x,y,N,M)
#include "integ2deps.inc"
#undef NUMFOR_INTEG
end PRNM(trapz2d_func_eps)

PRNM(simps2d_func_eps) (func, xrange, yrange, eps, limit) result(int)
implicit none
NUMFOR_KINDF :: func
NUMFOR_KINDR :: int !<
NUMFOR_KINDR :: int0
#define NUMFOR_INTEG(f,x,y,N,M) simps(f,x,y,N,M)
#include "integ2deps.inc"
#undef NUMFOR_INTEG
end PRNM(simps2d_func_eps)

