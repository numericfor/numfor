!> @ingroup integrate
!> Definition of integrable functions
!! Description: @ref docintegrate
module func_integ
  public

  ! These are the type of functions that will be integrated
  abstract interface            ! real f(x)
    function funqr(x) result(f) !< Function to integrate f(x)
      real(8), intent(IN) :: x  !< variable to be integrated
      real(8) :: f
    end function funqr
  end interface
  abstract interface            ! complex f(x)
    function funqc(x) result(f) !< Function to integrate
      real(8), intent(IN) :: x
      complex(8) :: f
    end function funqc
  end interface
  abstract interface              ! real f([x])
    function funqrnd(x) result(f) !< Function to integrate
      real(8), dimension(:), intent(IN) :: x
      real(8) :: f
    end function funqrnd
  end interface
  abstract interface              ! complex f([x])
    function funqcnd(x) result(f) !< Function to integrate
      real(8), dimension(:), intent(IN) :: x
      complex(8) :: f
    end function funqcnd
  end interface

  abstract interface                     ! real f(x, [args])
    function funqrarg(x, args) result(f) !< Function to integrate
      real(8), intent(IN) :: x
      real(8), dimension(:), intent(IN) :: args
      real(8) :: f
    end function funqrarg
  end interface
  abstract interface                     ! complex f(x, [args])
    function funqcarg(x, args) result(f) !< Function to integrate
      real(8), intent(IN) :: x
      real(8), dimension(:), intent(IN) :: args
      complex(8) :: f
    end function funqcarg
  end interface
  abstract interface                       ! real f([x], [args])
    function funqrndarg(x, args) result(f) !< Function to integrate
      real(8), dimension(:), intent(IN) :: x
      real(8), dimension(:), intent(IN) :: args
      real(8) :: f
    end function funqrndarg
  end interface
  abstract interface                       ! complex f([x], [args])
    function funqcndarg(x, args) result(f) !< Function to integrate
      real(8), dimension(:), intent(IN) :: x
      real(8), dimension(:), intent(IN) :: args
      complex(8) :: f
    end function funqcndarg
  end interface

  !> Type to encapsulate real functions and extra arguments
  type nf_rfunction
    procedure(funqrarg), nopass, pointer :: f !< real function f(x,args)
    real(8), dimension(:), allocatable :: args !< extra-arguments
  end type nf_rfunction
  !> Type to encapsulate complex functions and extra arguments
  type nf_cfunction
    procedure(funqcarg), nopass, pointer :: f!< complex function f(x,args)
    real(8), dimension(:), allocatable :: args!< extra-arguments
  end type nf_cfunction
  type nf_rndfunction
    procedure(funqrndarg), nopass, pointer :: f
    real(8), dimension(:), allocatable :: args
  end type nf_rndfunction
  type nf_cndfunction
    procedure(funqcndarg), nopass, pointer :: f
    real(8), dimension(:), allocatable :: args
  end type nf_cndfunction

  abstract interface
    !  Kind of function used as weight functions in quadrature
    function qfwght(x, args, iarg) result(y)
      implicit none
      real(8) :: y
      real(8), intent(IN) :: x !<
      real(8), dimension(:), intent(IN) :: args !<
      integer, intent(IN) :: iarg !<
    end function qfwght
  end interface

end module func_integ
