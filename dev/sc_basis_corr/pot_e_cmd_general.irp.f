  BEGIN_PROVIDER [double precision, pot_basis_alpha_ao,(ao_num,ao_num,N_states)]
 &BEGIN_PROVIDER [double precision, pot_basis_beta_ao,(ao_num,ao_num,N_states)]
  implicit none
  BEGIN_DOC
 ! general providers for the alpha/beta correlation potentials on the AO basis
  END_DOC

 BEGIN_SHELL [ /usr/bin/env python3 ]
import os
import glob
from qp_path import QP_SRC
funcdir=QP_SRC+'/../plugins/qp_plugins_lct/dev/func_mu_of_r/'
os.chdir(funcdir)
functionals_tmp = map(lambda x : x.replace(".irp.f","") , glob.glob("func_*.irp.f"))
functionals = map(lambda x : x.replace("func_","") , functionals_tmp)                 

prefix = ""
for f in functionals:
  print("""
  %sif (trim(basis_cor_func) == '%s') then
    pot_basis_alpha_ao = pot_basis_alpha_ao_%s
    pot_basis_beta_ao  = pot_basis_beta_ao_%s"""%(prefix, f, f, f) )
  prefix = "else "

print("""
  else
   print*, 'Basis set Correlation functional required does not exist ...'
   stop""" )
print("endif")

 END_SHELL

 END_PROVIDER

  BEGIN_PROVIDER [logical , needs_eff_two_e_ints ]
  implicit none
  BEGIN_DOC
 ! If true then it tells you that you would need to write effective two-e integrals
  END_DOC
    print*,'needs_eff_two_e_ints_su_pbe_ot before = ',needs_eff_two_e_ints_su_pbe_ot

 BEGIN_SHELL [ /usr/bin/env python3 ]
import os
import glob
from qp_path import QP_SRC
funcdir=QP_SRC+'/../plugins/qp_plugins_lct/dev/func_mu_of_r/'
os.chdir(funcdir)
functionals_tmp = map(lambda x : x.replace(".irp.f","") , glob.glob("func_*.irp.f"))
functionals = map(lambda x : x.replace("func_","") , functionals_tmp)                 

prefix = ""
for f in functionals:
  print("""
  %sif (trim(basis_cor_func) == '%s') then
    needs_eff_two_e_ints = needs_eff_two_e_ints_%s
    """%(prefix, f, f) )
  prefix = "else "

print("""
  else
   print*, 'Correlation functional required does not exist ...'
   stop""" )
print("endif")

 END_SHELL


  END_PROVIDER 

  BEGIN_PROVIDER [double precision, d_dn2_e_cmd_basis, (n_points_final_grid,N_states)]
  implicit none
  BEGIN_DOC
 ! general provider for the functional derivative of the on-top in real space
 !
 ! this is needed for the effective two-e integrals 
  END_DOC

 BEGIN_SHELL [ /usr/bin/env python3 ]
import os
import glob
from qp_path import QP_SRC
funcdir=QP_SRC+'/../plugins/qp_plugins_lct/dev/func_mu_of_r/'
os.chdir(funcdir)
functionals_tmp = map(lambda x : x.replace(".irp.f","") , glob.glob("func_*.irp.f"))
functionals = map(lambda x : x.replace("func_","") , functionals_tmp)                 

prefix = ""
for f in functionals:
  print("""
  %sif (trim(basis_cor_func) == '%s') then
    d_dn2_e_cmd_basis = d_dn2_e_cmd_%s"""%(prefix, f, f) )
  prefix = "else "

print("""
  else
   print*, 'Correlation functional required does not exist ...'
   stop""" )
print("endif")

 END_SHELL

 END_PROVIDER



  BEGIN_PROVIDER [double precision, pot_basis_alpha_mo,(mo_num,mo_num,N_states)]
 &BEGIN_PROVIDER [double precision, pot_basis_beta_mo, (mo_num,mo_num,N_states)]
  implicit none
  BEGIN_DOC
 ! general providers for the alpha/beta correlation potentials on the MO basis
  END_DOC
  integer :: istate
  do istate = 1, N_states
     call ao_to_mo(                                                   &
         pot_basis_alpha_ao(1,1,istate),                                 &
         size(pot_basis_alpha_ao,1),                                &
         pot_basis_alpha_mo(1,1,istate),                                 &
         size(pot_basis_alpha_mo,1)                                 &
         )

     call ao_to_mo(                                                   &
         pot_basis_beta_ao(1,1,istate),                                  &
         size(pot_basis_beta_ao,1),                                 &
         pot_basis_beta_mo(1,1,istate),                                  &
         size(pot_basis_beta_mo,1)                                  &
         )
  enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, int_d_dn2_e_cmd_basis, (N_states)]
 implicit none
 integer :: istate,ipoint
 double precision :: weight
 int_d_dn2_e_cmd_basis = 0.d0
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)
   int_d_dn2_e_cmd_basis(istate) += d_dn2_e_cmd_basis(ipoint,istate) * weight
  enddo
 enddo
 END_PROVIDER 
