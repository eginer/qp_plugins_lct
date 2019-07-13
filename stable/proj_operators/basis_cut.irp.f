
BEGIN_PROVIDER [ integer, n_orb_max_basis  ]
  implicit none
  BEGIN_DOC
! maximum number of MOs to define the basis set
  END_DOC
  logical                        :: has
 if(mu_of_r_potential.EQ."psi_cas_truncated")then
  n_orb_max_basis = n_core_inact_act_orb
 else 
  n_orb_max_basis = mo_num
 endif

END_PROVIDER
