from scipy.constants import physical_constants
import re
notice = """
## Notice

Data taken from scipy.constants (see file codata.py)

The values of the constants provided at this site are recommended for
international use by CODATA and are the latest available. Termed the "2014
CODATA recommended values," they are generally recognized worldwide for use in
all fields of science and technology. The values became available on 25 June
2015 and replaced the 2010 CODATA set. They are based on all of the data
available through 31 December 2014. The 2014 adjustment was carried out under
the auspices of the CODATA Task Group on Fundamental Constants. Also available
is an introduction to the constants for non-experts at
https://physics.nist.gov/cuu/Constants/introduction.html

**Reference:**    https://physics.nist.gov/cuu/Constants/

"""

modeline = "! -*- mode: F90 -*-\n"

letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
digits_ = '0123456789' + '_'
valid_fortran_set = letters + digits_


def is_valid_fortran(s):
  isf = all([x in valid_fortran_set for x in s]) and (s[0] not in digits_)
  return isf


# Process for Fortran mode
names = []
descriptions = []
values = []
for k, v in physical_constants.items():
  name = re.sub(r'(\W)(?=\1)', '\1', k)
  name = name.replace('.', '').replace('-', '_').replace(' ', '_')
  # Check if is a valid Fortran name
  if is_valid_fortran(name):
    names.append(name)
    descriptions.append(" {0} ({2}) {1} ".format(*v))
    values.append("{}_dp".format(v[0]))


# Width of columns
wcol1 = len(max(names, key=len))
wcol2 = len(max(descriptions, key=len))

# headers
doc_txt = """## Physical constants {#physconst}

CODATA Recommended Values of the Fundamental Physical Constants 2014.


"""
doc_txt += '|{}|{}|\n'.format('Constant Name'.center(wcol1 + 2),
                              'Constant value'.center(wcol2 + 2))
doc_txt += '|{}|{}|\n'.format((wcol1 + 2) * '-', (wcol2 + 2) * '-')

inc_txt = modeline

# Create the tables/list of constants
for k, v, d in sorted(zip(names, values, descriptions),
                      key=lambda s: s[0].lower()):
  inc_txt += "real(dp), public, parameter :: {} = {}\n".format(k, v)
  doc_txt += "| {} | {} |\n".format(k.ljust(wcol1), d.ljust(wcol2))

# Write to files
fdir = 'src/utils/'
ddir = 'docs/sources/'

forfname = os.path.join(fdir, "codata.inc")
with open(forfname, 'w') as fo:
  fo.write(inc_txt)

mdfname = os.path.join(ddir, "constants.md")
with open(mdfname, 'w') as fo:
  fo.write(doc_txt + notice)
