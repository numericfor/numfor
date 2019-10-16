# Submodule Utils {#docutils}

The module utils defines several objects (types, constants and routines) convenient for common tasks in programming.

It consists of several areas:

  * Basic support (basic types, error handling, timer).

  * Manipulation of strings is provided in two different "flavors". Both modules provide similar functionality, and where posible they are implemented in the same code. It consists of:
      * A framework to work in a "procedural" approach is provided.
      * An "object-oriented" framework with a new type (`fStr()`) to represent strings.

## General support

  * Definitions for single precision (`sp`), double precision (`dp`) kinds
  * A timer, with a simple interface
  
@include ex_timer.f90

Prints:
```{.shell-session}
cpu time:   7.77s
Total time: 7 s 801ms 
```
## Constants


There are defined some commonly used constants (where the definition of "commonly" is quite arbitrary).


| Name       | Value                                     |
|------------|-------------------------------------------|
| Zero       | 0._dp                                     |
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


There are also defined many @subpage physconst, like the constants defined in [scipy.constants](https://docs.scipy.org/doc/scipy/reference/constants.html).

## Strings ##

The module strings provides routines for string manipulation. It provides some functionality similar to that available to Python strings:

  + `strings::str()`: Cast numbers to strings
  
There are several functions that act similarly to Python string methods:
  
  + `strings::upper()`: Converts a string to uppercase
  + `strings::lower()`: Converts a string to lowercase
  + `strings::swapcase()`: Swap cases: (lower -> upper) and (upper -> lower)
  + `strings::reverse()`: Reverses a string "abc" -> "cba"
  + `strings::endswith()`: Test wheter a string ends with a given suffix
  + `strings::startswith()`: Test wheter a string starts with a given prefix
  + `strings::lstrip()`: Removes leading characters of a string
  + `strings::rstrip()`: Removes trailing characters of a string
  + `strings::strip()`: Removes leading and trailing characters of a string
  + `strings::count_sub()` Counts how many times a substring is in a string
  + `strings::find()`: Finds a substring in a string
  + `strings::issub()`: Test whether a substring is in a string
  + `strings::count_sub()`: Counts the number of occurrences of a substring in a string
  + `strings::lstrip()`: Removes leading chars in a string
  + `strings::rstrip()`: Removes trailing chars in a string
  + `strings::strip()`: Removes leading and trailing chars in a string
  + `strings::ljust()`: Left-padds a string with a char
  + `strings::rjust()`: Right-padds a string with a char
  + `strings::zfill()`: Padds with zeroes a string
  + `strings::center()`: Centers a string in a line
  + `strings::find()`: Finds a substring in a string
  + `strings::replace()`: Replaces a substring with other in a string

### Examples of use {#examplestrings} ###


### A simple example {#ex_strings1}

This is a very simple example, of how to create (output) filenames from the parameters of our problem, showcasing only a few functions.

@include ex_strings1.f90

Note here the use of `strings::strip()` refers to the above function.


After compiling and executing:
Output for a very long directory name:
```bash
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

```bash
 Enter directory name
dd

======= These will the files created:======
           dd/out1.3cm42eV_1.dat
           dd/out1.3cm42eV_2.dat
           dd/out1.3cm42eV_3.dat
           dd/out1.3cm42eV_4.dat
           dd/out1.3cm42eV_5.dat
```

## Object oriented strings fstring ##

This module provides an object-oriented version of @ref modstrings, with functionality similar to that available to Python strings:
  + `fStr()`: Cast numbers to strings (of type `fStr`)
  + Strings may be read and write with default format
  + Multiplication by integers produce repetition
  + They can be compared for equality (or not)


The use of this module is rather simple: we define a variable of the type `fStr`, and then use its methods. 


### Examples of use {#examplefstring} ###

A minimal example would be simply:

```{.f90}
USE fString, only: fStr
type(fStr) :: s1
type(fStr) :: s2
s1 = '     Variable s1 =    '
s2 = fStr(1.2123)
print *, 'The '//s1%strip(' =')//' values' // s2
```

which outputs:
```{.bash}
"The Variable s1 values 1.21229994"
```

Note here that the method `fstring::strip()` removes all leading and trailing spaces and occurrences of the `"="` symbol.

#### A simple example {#ex_fstring1} ####

The Object Oriented version of examples presented in `strings` (@ref ex_strings1) has a very similar form:

@include ex_fstring1.f90

which, after compiling and executing, gives exactly the same result that \ref ex_strings1 from module `strings`

