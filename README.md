
# Numeric Modern Fortran

The goal of this project is to create a Fortran library similar to [Numpy](https://www.numpy.org)/[Scipy](https://www.scipy.org) for scientific computing.

This library consists (will consist) of several modules to ease the scientific work. It aims to provide user-friendly utilities and relatively-efficient routines for scientific computing, numerical and related work.

## Documentation


Documentation may be found in https://numericfor.github.io/numfor/index.html


Further information on the capabilities of the library may be learned by exploring the different modules.



It provides:
  + A @subpage docutils utils with basic, non-specific, functionality used in many (most?) scientific programming.
  + A @subpage docarrays with basic functionality to work with arrays, including generation of grids, search of elements and sorting.
  + A @subpage docinterpolate with routines to perform interpolation, fitting of data, and some work on polynomials.
  + A @subpage docintegrate with routines to perform integration of real and complex functions, and of sampled data.


You can also visit the [documentation for the code](namespaces.html)

## Installation ##

The installation follow the usual steps. The simpler procedure is:

  * `make`
  Builds the module information file and static library

  * `make tests` (optional, not yet fully implemented)
  Run tests on submodules

  * `make install`
  Installs the library, modules and documentation into `<prefix>` path


## Use ##

Simply compile your program with flags indicating where to find libnumfor.mod and link to libnumfor.a
There are different ways to accomplish this:


  1. Using `pkg-config`
  ```bash
  $> gfortran -c -o myprog.o $(pkg-config --cflags numfor) myprog.f90
  $> gfortran -o myprog $(pkg-config --libs numfor) myprog.o
  ``` 
  In order to work, first the environment variable PKG_CONFIG_PATH must be set. 
  For instance, in linux you have to add to $HOME/.bashrc (or site /etc/bashrc if you install for all users) the folowing line:
  
  ```bash
  export PKG_CONFIG_PATH="<prefix>/lib/pkgconfig:$PKG_CONFIG_PATH"
  ```
  You will have to open a new terminal or source the `bashrc` file before compiling your program.
  
  2. Manually adding the information when compiling. For instance, when using `gfortran`
  
  ```bash
  $> gfortran -c -o myprog.o -I <prefix>/include/numfor myprog.f90
  $> gfortran -o myprog -L <prefix>/lib -lnumfor myprog.o
  ```




## Make the docs ##

To build the html docs just:

```bash
$> make doc
```

Some dependencies are needed:

 - [Doxygen](http://www.doxygen.nl/)
 - A LaTeX distribution
 
