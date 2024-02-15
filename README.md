
# Numeric Modern Fortran

The goal of this project is to create a Fortran library similar to [Numpy](https://www.numpy.org)/[Scipy](https://www.scipy.org) for scientific computing.

This library consists (or rather will eventually consist) of several modules to ease the scientific work. It aims to provide user-friendly utilities and relatively-efficient routines for scientific computing, numerical and related work.

## Documentation


Documentation may be found in https://numericfor.github.io/numfor/index.html

Further information on the capabilities of the library may be learned by exploring the different [modules](https://numericfor.github.io/numfor/namespaces.html).

It provides:
  + A "submodule" **Utils** with basic, non-specific, functionality used in many (most?) scientific programs.
  + A "submodule" **Array_utils** with basic functionality to work with arrays, including generation of grids, search of elements and sorting.
  + A "submodule" **Interpolate** with routines to perform interpolation, fitting of data, and some work on polynomials.
  + A "submodule" **Integrate** with routines to perform integration of real and complex functions, and of sampled data.


## Installation ##

Firstly, download the library:

  * `git clone https://github.com/numericfor/numfor.git`

The installation follow the usual steps. The simpler procedure is:

  * cd `numfor`
  *  `make`
  Builds the module information file and static library

  * `make tests` (optional, not yet fully implemented)
  Run tests on submodules

  * `make install`
  Installs the library, modules and documentation into `<prefix>` path


## Use ##

Simply compile your program with flags indicating where to find libnumfor.mod and link to libnumfor.a
There are different ways to accomplish this:

  1. Manually adding the information when compiling. For instance, when using `gfortran`
  
  ```bash
  $> gfortran -c -o myprog.o -I <prefix>/include/numfor myprog.f90
  $> gfortran -o myprog -L <prefix>/lib -lnumfor myprog.o
  ```

  2. Using `pkg-config`
  ```bash
  $> gfortran [your FFLAGS] -c -o myprog.o $(pkg-config --cflags numfor) myprog.f90
  $> gfortran -o myprog myprog.o $(pkg-config --libs numfor)
  ``` 
  or, simply in one step:
  
  ```bash
  $> gfortran [your FFLAGS] -o myprog $(pkg-config --cflags numfor) myprog.f90 $(pkg-config --libs numfor)
  ``` 

  **Note** that the library options must be, in both cases, at the end of the command.

  In order to work, first the environment variable PKG_CONFIG_PATH must be set. 
  For instance, in linux you have to add to $HOME/.bashrc (or site /etc/bashrc if you install for all users) the folowing line:
  
  ```bash
  export PKG_CONFIG_PATH="<prefix>/lib/pkgconfig:$PKG_CONFIG_PATH"
  ```
  You will need to open a new terminal or source the `bashrc` file before compiling your program.
  


