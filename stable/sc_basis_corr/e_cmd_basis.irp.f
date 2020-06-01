 BEGIN_PROVIDER [double precision, e_c_md_basis, (N_states)]
  implicit none
  BEGIN_DOC
  ! Ecmd for the basis set correction
  END_DOC

 BEGIN_SHELL [ /usr/bin/env python3 ]
import os
import glob
from qp_path import QP_SRC
funcdir=QP_SRC+'/../plugins/qp_plugins_lct/stable/func_mu_of_r/'
os.chdir(funcdir)
functionals_tmp = map(lambda x : x.replace(".irp.f","") , glob.glob("func_*.irp.f"))
functionals = map(lambda x : x.replace("func_","") , functionals_tmp)                 
prefix = ""
for f in functionals:
  print("""
  %sif (trim(basis_cor_func) == '%s') then
    e_c_md_basis = e_c_md_basis_%s"""%(prefix, f, f) )
  prefix = "else "

print("""
  else
   print*, 'Basis set Correlation functional required does not exist ...'
   print*,'correlation_functional ',correlation_functional
   stop""")
print("endif")

 END_SHELL

 END_PROVIDER

