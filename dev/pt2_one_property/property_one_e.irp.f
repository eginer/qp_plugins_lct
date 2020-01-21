BEGIN_PROVIDER [double precision, one_prop_pot_b_provider, (mo_num,mo_num)]
 implicit none
 one_prop_pot_b_provider = -1.d0 * mo_dipole_z
END_PROVIDER 

BEGIN_PROVIDER [double precision, one_prop_pot_a_provider, (mo_num,mo_num)]
 implicit none
 one_prop_pot_a_provider = -1.d0 * mo_dipole_z
END_PROVIDER 
