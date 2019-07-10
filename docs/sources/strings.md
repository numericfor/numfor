# Module strings {#modstrings}

##  Description

The module strings provides routines for string manipulation. It provides some functionality similar to that available to Python strings:

  + `str()`: Cast numbers to strings
  + Multiplication by integers produce repetition
      + `left_multiply_int()`: Multiply an integer by a string 
      + `right_multiply_int()`: Multiply a string by an integer

  
There are functions that act similarly to Python string methods:
  
  + `upper()`: Converts a string to uppercase
  + `lower()`: Converts a string to lowercase
  + `swapcase()`: Swap cases: (lower -> upper) and (upper -> lower)
  + `reverse()`: Reverses a string "abc" -> "cba"
  + `endswith()`: Test wheter a string ends with a given suffix
  + `startswith()`: Test wheter a string starts with a given prefix
  + `lstrip()`: Removes leading characters of a string
  + `rstrip()`: Removes trailing characters of a string
  + `strip()`: Removes leading and trailing characters of a string
  + `count_sub()` Counts how many times a substring is in a string
  + `find()`: Finds a substring in a string
  + `issub()`: Test whether a substring is in a string
  + `count_sub()`: Counts the number of occurrences of a substring in a string
  + `lstrip()`: Removes leading chars in a string
  + `rstrip()`: Removes trailing chars in a string
  + `strip()`: Removes leading and trailing chars in a string
  + `ljust()`: Left-padds a string with a char
  + `rjust()`: Right-padds a string with a char
  + `zfill()`: Padds with zeroes a string
  + `center()`: Centers a string in a line
  + `find()`: Finds a substring in a string
  + `replace()`: Replaces a substring with other in a string

# Examples of use {#examplestrings}


## A simple example {#ex_strings1}

This is a very simple example, of how to create (output) filenames from the parameters of our problem, showcasing only a few functions.

\include ex_strings1.f90

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
