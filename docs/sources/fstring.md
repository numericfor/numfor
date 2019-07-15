# Submodule fString {#modfstring}
The module fString provides the type `fStr()`, suitable for string manipulation.

##  Description 

 It provides an object-oriented version of @ref modstrings, with functionality similar to that available to Python strings:
  + `fStr()`: Cast numbers to strings (of type `fStr`)
  + Strings may be read and write with default format
  + Multiplication by integers produce repetition
  + They can be compared for equality (or not)


The use of this module is rather simple: we define a variable of the type `fStr`, and then use its methods. 


# Examples of use {#examplefstring}

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
```{.f90}
"The Variable s1 values 1.21229994"
```

Note here that the method `fstring::strip()` removes all leading and trailing spaces and occurrences of the `"="` symbol.

## A simple example {#ex_fstring1}

The Object Oriented version of examples presented in `strings` (@ref ex_strings1) has a very similar form:

\include ex_fstring1.f90

which, after compiling and executing, gives exactly the same result that \ref ex_strings1 from module `strings`

