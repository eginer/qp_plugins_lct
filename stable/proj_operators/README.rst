=========================
projected_operators_utils
=========================

This plugin proposes the various routines/providers to compute the function f_{\Psi^B}(r_1,r_2) of Eq. (22) of J. Chem. Phys.149, 194301 (2018). 

This quantity if fundamental to determine the basis set correction (see same ref for more details). 

This quantity f_{\Psi^B}(r_1,r_2) is computed for TWO TYPES OF \Psi_^B: 
  
   +) HF-like wave function  : a SINGLE Slater determinant built with RESTRICTED ORBITALS contains in "mo_basis"

   +) CAS-like wave function : a LINEAR COMBINATION of Slater determinants all belonging to an ACTIVE SPACE defined by lis_core, list_inact, list_act and so on ...
