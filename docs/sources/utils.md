# Module utils {#utils}

The module utils defines several objects (types, constants and routines) convenient for common tasks in programming.

It consists of several areas:

  * Basic support (basic types, error handling, timer) is provided by @subpage modbasic.

  * Manipulation of strings is provided in two different "flavors". Both modules provide similar functionality, and where posible they are implemented in the same code. It consists of:
      * A framework to work in a "procedural" approach is provided by the functions and subroutines in the @subpage modstrings.
      * An "object-oriented" framework with a new type (`fStr()`) to represent strings. It is defined in the @subpage modfstring.

