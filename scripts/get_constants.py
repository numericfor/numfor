from scipy.constants import physical_constants
for k, v in physical_constants.items():
  print("real(dp), public, parameter :: {} = {}_dp".format(
      k.replace('.', '').replace(' ', '_'), v[0]))
