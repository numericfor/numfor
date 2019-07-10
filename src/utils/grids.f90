!> This module provides convenience routines to work with grids and arrays
module grids
  implicit none
  real(8), parameter :: def_base = 10._8

  private
  public :: linspace, logspace, arange

contains

  !> Return evenly spaced numbers over a specified interval
  !!
  !! Returns `num` evenly spaced samples, calculated over the interval [`start`, `stop`].
  !!
  !!  Examples
  !!  ```
  !!  linspace(2.0, 3.0, num=5)
  !!  ! gives:  [ 2.  ,  2.25,  2.5 ,  2.75,  3.  ]
  !!  linspace(2.0, 3.0, num=5, endpoint=False)
  !!  ! gives: [ 2. ,  2.2,  2.4,  2.6,  2.8]
  !!  linspace(2.0, 3.0, num=5, retstep=True)
  !!  ! gives: [ 2.  ,  2.25,  2.5 ,  2.75,  3.  ]
  !!  ! and retstep= 0.25
  !!  !
  !! ```
  !!
  function linspace(start, end, num, endpoint, retstep) result(x)
    implicit none
    real(8), intent(IN) :: start !< The starting value of the sequence.
    real(8), intent(IN) :: end   !< The end value of the sequence,
    integer, intent(IN) :: num !< Number of samples to generate. Must be positive.
    logical, optional, intent(IN) :: endpoint !< If True, `end` is the last sample. Otherwise, it is not included. Default is True
    real(8), optional, intent(OUT) :: retstep !< If present, return the step
    real(8), dimension(num) :: x              !< An array of uniformly spaced numbers
    real(8) :: step
    integer(4) :: i
    x = Zero
    IF (num < 1) return
    x(1) = start

    if (present(endpoint) .and. (.not. endpoint)) then
      step = (end - start) / real(num, kind=8)
    else
      step = (end - start) / (num - 1._8)
    end if

    IF (present(retstep)) retstep = step
    IF (num == 1) return

    ! Version with forall
    forall (i=1, num - 1) x(i + 1) = start + step * i
    ! do i = 1, num - 1
    !   x(i + 1) = start + step * i
    ! end do

  end function linspace

  !> Makes a grid with numbers spaced evenly on a log scale
  !!
  !! In linear space, the sequence starts at ``base**start``
  !! (`base` to the power of `start`) and ends with ``base**end``
  !!
  !! Examples
  !!
  !! ```
  !! logspace(2.0, 3.0, num=4)
  !! ! gives: [  100.        ,   215.443469  ,   464.15888336,  1000.        ]
  !! logspace(2.0, 3.0, num=4, endpoint=False)
  !! ! gives: [ 100.        ,  177.827941  ,  316.22776602,  562.34132519]
  !! logspace(2.0, 3.0, num=4, base=2.0)
  !! ! gives: array([ 4.        ,  5.0396842 ,  6.34960421,  8.        ])
  !! !
  !! ```
  !!
  function logspace(start, end, num, endpoint, base) result(x)
    implicit none
    real(8), intent(IN) :: start !< ``base**start`` is the starting value of the sequence.
    real(8), intent(IN) :: end   !< ``base**end`` is the final value of the sequence.
    integer, optional, intent(IN) :: num !< Number of samples to generate. Must be positive.
    logical, optional, intent(IN) :: endpoint !< If True, `end` is the last sample. Otherwise, it is not included. Default is True
    real(8), optional, intent(IN) :: base     !< The base of the log space. Default is 10.
    real(8), dimension(num) :: x              !< A sequence of numbers spaced evenly on a log scale.

    real(8) :: b_
    real(8) :: step
    integer(4) :: i

    x = Zero
    IF (num < 1) return

    b_ = def_base; IF (present(base)) b_ = base
    x(1) = b_**start
    IF (num == 1) return
    if (present(endpoint) .and. (.not. endpoint)) then
      step = (end - start) / real(num, kind=8)
    else
      step = (end - start) / (num - 1._8)
    end if
    forall (i=1, num - 1) x(i + 1) = b**(start + step * i)
    ! do i = 1, num - 1
    !   x(i + 1) = b**(start + step * i)
    ! end do
  end function logspace

  !> arange: Return evenly spaced integer values within a given interval
  !!
  !! Values are generated within the half-open interval ``[start, end)``
  !! (in other words, the interval including `start` but excluding `end`).
  function arange(start, end, step) result(x)
    implicit none
    integer, intent(IN) :: start !< the starting value of the interval.
    integer, intent(IN) :: end   !< the final value of the interval (not included)
    integer, optional, intent(IN) :: step !< Spacing between values.
    integer, dimension(:), allocatable :: x !< A sequence of numbers spaced evenly
    integer :: num
    !
    IF (step == 0) return
    step_ = 1; IF (present(step)) step_ = step
    num = ceiling((end - start) / step)
    IF (allocated(x) .and. (size(x) /= num)) deallocate (x)
    IF (.not. allocated(x)) allocate (x(num))
    x(1) = start
    do i = 1, num
      x(i + 1) = x(i) + step_
    end do
  end function arange

!> savetxt Guarda un array 2D en un archivo de texto
!!
!! @note
!! Si fname es "stdout" o " ", o no están presente ni fname ni unit, usa stdout
!! Si se da fname el archivo se abre y cierra.
!! Si se da unit, el archivo queda abierto
  subroutine savetxt(a, fmt, fname, unit)
    implicit none
    real(8), dimension(:, :), intent(IN) :: a        !< Array a escribir a archivo de texto
    character(len=*), optional, intent(in) :: fmt    !< formato a usar para los datos. Default 'g0.5'
    character(len=*), optional, intent(in) :: fname  !< Nombre del archivo de salida
    integer, optional, intent(in) :: unit            !< Unidad a escribir si el archivo está abierto

    real(8), dimension(ubound(a, 2), ubound(a, 1)) :: b
    integer, dimension(2) :: sh
    integer :: i
    integer :: u

    character(len=32) :: form = "g0.5" ! Default
    character(len=32) :: formato
    logical :: closef

    ! Si fname está presente => Toma precedencia sobre unit. Si ninguna está usa stdout
    closef = .False.
    u = stdout
    if (Present(fname)) then
      if (trim(fname) /= '' .and. trim(fname) /= 'stdout') then
        open (newunit=u, file=trim(fname))
        closef = .True.
      end if
    else if (Present(unit)) then ! The file was already open before invoking the function
      IF (unit >= 0 .and. unit <= 99) u = unit
    end if

    b = transpose(a)
    sh = shape(b)

    if (Present(fmt) .and. (trim(fmt) /= 'default') .and. (trim(fmt) /= '')) then
      if (index('(', fmt) == 0) then
        write (formato, '(A,I1,A,A,A)') '(', sh(1), '(', trim(fmt), ',1x))'
      else
        formato = fmt
      end if
    else
      write (formato, '(A,I1,A,A,A)') '(', sh(1), '(', trim(form), ',1x))'
    end if
    do i = 1, sh(2)
      write (u, formato) b(:, i)
    end do

    ! write (u, formato) (b(:, i), i=1, sh(2))

    IF (closef) close (u)
  end subroutine savetxt

end module grids
