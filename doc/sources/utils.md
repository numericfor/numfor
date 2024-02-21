# Submodule Utils

The module utils defines several objects (types, constants and routines) convenient for common tasks in programming.

It consists of several areas:

  * Basic support (basic types, error handling, timer).
      * Definitions for single precision (`sp`), double precision (`dp`) kinds
      * Definition of positive and negative infinite
      * A timer, with a simple interface

  * Manipulation of strings is provided in two different "flavors". Both modules provide similar functionality, and where posible they are implemented in the same code. It consists of:
      * A framework to work in a "procedural" approach is provided.
      * An "object-oriented" framework with a new type (`fStr()`) to represent strings.

## General support


For numerical work we need to use specific minimal precision which is
machine-and-compiler-independent. The real type definitions emphatize this.
The names are `sp` and `dp`, and are chosen only by tradition/historic reasons.
For double precision we define types with at least 15 decimal places and exponent range of 307.
(see for instance http://fortranwiki.org/fortran/show/Real+precision)
  
It is convenient for several routines to have a symbol indicating infinite values. We define `nf_minf` and `nf_inf` for negative and positive infinite, respectively.



## Constants


There are defined some commonly used constants (where the definition of "commonly" is quite arbitrary).


| Name       | Value                                     |
|------------|-------------------------------------------|
| Zero       | 0._dp                                     |
| nf_minf    | **-infinite**                             |
| nf_inf     | **infinite**                              |
| C0         | (0._dp, 0._dp)                            |
| C1         | (1._dp, 0._dp)                            |
| CI         | (0._dp, 1._dp)                            |
| M_pi       | 3.14159265358979323846264338328_dp        |
| M_dpi      | 6.28318530717958647692528676656_dp        |
| M_pi_2     | 1.57079632679489661923132169164_dp        |
| M_pi_4     | 0.78539816339744830966156608458_dp        |
| M_sqrtpi   | 1.77245385090551602729816748334_dp        |
| M_2_sqrtpi | 1.12837916709551257389615890312_dp        |
| M_1_pi     | 0.31830988618379067153776752675_dp        |
| M_2_pi     | 0.63661977236758134307553505349_dp        |
| M_lnpi     | 1.14472988584940017414342735135_dp        |
| M_ln2      | 0.693147180559945309417232121458176568_dp |
| M_e        | 2.71828182845904523536028747135_dp        |
| M_log2e    | 1.442695040888963407359924681001892137_dp |
| M_log10e   | 0.43429448190325182765112891892_dp        |


There are also defined many @subpage physical-constants, which are mostly those defined in [scipy.constants](https://docs.scipy.org/doc/scipy/reference/constants.html).

## Simple Timing 

We implemented a simple timer that may be use by calling the methods `start()`, `stop()` and `show()`, as shown in the following example:

@include ex_timer.f90

Prints:
```{.shell-session}
cpu time:   7.77s
Total time: 7 s 801ms 
```

## Warning and error messages

There is currently only one routine implemented for warning and error messages, [print_msg()](@ref basic::print_msg()), that may be called to produce informative results to `stdout` or an alternative place (usually an already opened file), and optionally stop the program if the error code is positive.

For instance,

```{.f90}
USE numfor, only: print_msg
! Do something ...
! and after an unexpected result
call print_msg("Result outside of expected range!") 
```
will print the message and keep the program running, while
```{.f90}
USE numfor, only: print_msg
call print_msg("Invalid argument", "my-routine", 1) 
```
will print the message, the routine and stop the program with error code 1.


## Strings ##

The module strings provides routines for string manipulation. It provides some functionality similar to that available to Python strings:

  + [str()](@ref strings::str()): Cast numbers to strings
  
There are several functions that act similarly to Python string methods:
  
  + [upper()](@ref strings::upper()): Converts a string to uppercase
  + [lower()](@ref strings::lower()): Converts a string to lowercase
  + [swapcase()](@ref strings::swapcase()): Swap cases: (lower -> upper) and (upper -> lower)
  + [reverse()](@ref strings::reverse()): Reverses a string "abc" -> "cba"
  + [endswith()](@ref strings::endswith()): Test wheter a string ends with a given suffix
  + [startswith()](@ref strings::startswith()): Test wheter a string starts with a given prefix
  + [lstrip()](@ref strings::lstrip()): Removes leading characters of a string
  + [rstrip()](@ref strings::rstrip()): Removes trailing characters of a string
  + [strip()](@ref strings::strip()): Removes leading and trailing characters of a string
  + [count_sub()](@ref strings::count_sub()) Counts how many times a substring is in a string
  + [find()](@ref strings::find()): Finds a substring in a string
  + [issub()](@ref strings::issub()): Test whether a substring is in a string
  + [count_sub()](@ref strings::count_sub()): Counts the number of occurrences of a substring in a string
  + [lstrip()](@ref strings::lstrip()): Removes leading chars in a string
  + [rstrip()](@ref strings::rstrip()): Removes trailing chars in a string
  + [strip()](@ref strings::strip()): Removes leading and trailing chars in a string
  + [ljust()](@ref strings::ljust()): Left-padds a string with a char
  + [rjust()](@ref strings::rjust()): Right-padds a string with a char
  + [zfill()](@ref strings::zfill()): Padds with zeroes a string
  + [center()](@ref strings::center()): Centers a string in a line
  + [find()](@ref strings::find()): Finds a substring in a string
  + [replace()](@ref strings::replace()): Replaces a substring with other in a string

### @anchor examplestrings Examples of use  ###


#### @anchor ex_strings1 A simple example 

This is a very simple example, of how to create (output) filenames from the parameters of our problem, showcasing only a few functions.

@include ex_strings1.f90

Note here the use of [strip()](@ref strings::strip()) refers to the above function.


After compiling and executing:
Output for a very long directory name:
```{.shell-session}
 Enter directory name

VeryLongDirectoryNamewithCamelCaseName

========================= These will the files created:========================
           verylongdirectorynamewithcamelcasename/out1.3cm42eV_1.dat
           verylongdirectorynamewithcamelcasename/out1.3cm42eV_2.dat
           verylongdirectorynamewithcamelcasename/out1.3cm42eV_3.dat
           verylongdirectorynamewithcamelcasename/out1.3cm42eV_4.dat
           verylongdirectorynamewithcamelcasename/out1.3cm42eV_5.dat
```

while for a short directory name the output is:

```{.shell-session}
 Enter directory name
dd

======= These will the files created:======
           dd/out1.3cm42eV_1.dat
           dd/out1.3cm42eV_2.dat
           dd/out1.3cm42eV_3.dat
           dd/out1.3cm42eV_4.dat
           dd/out1.3cm42eV_5.dat
```

## Object oriented strings fstring

This module provides an object-oriented version of the module `strings`, with functionality similar to that available to Python strings:
  + `fStr()`: Create type and also cast numbers to strings (of type `fStr`)
  + Strings may be read and write with default format
  + Multiplication by integers produce repetition
  + They can be compared for equality (or not)


The use of this module is rather simple: we define a variable of the type `fStr`, and then use its methods. 


### @anchor examplefstring Examples of use

A minimal example would be simply:

```{.f90}
USE numfor, only: fStr
type(fStr) :: s1
type(fStr) :: s2
s1 = '     Variable s1 =    '
s2 = fStr(1.2123)
print *, 'The '//s1%strip(' =')//' value is ' // s2
```

which outputs:

```{.bash}
"The Variable s1 value is 1.21229994"
```

Note here that the method `strip()` removes all leading and trailing spaces **and** occurrences of the equal symbol `"="`.

#### @anchor ex_fstring1 A simple example 

The Object Oriented version of examples presented in `strings` (@ref ex_strings1) has a very similar form:

@include ex_fstring1.f90

which, after compiling and executing, gives exactly the same result that \ref ex_strings1 from module `strings`

