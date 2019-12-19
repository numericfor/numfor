PRJ=numfor
version=1.1
prefix=/hola
FC=gfortran
rep="sed -e 's/@PRJ@/${PRJ}/' -e 's/@version@/${version}/' -e 's/@prefix@/${prefix}/' -e 's/@FC@/${FC}/'"

echo "$rep"
