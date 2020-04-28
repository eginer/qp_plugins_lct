BEGIN_PROVIDER [double precision, psi_energy_wee_sr, (N_states)]
 implicit none
 BEGIN_DOC
! <psi | W_ee^{sr} | psi> = <psi | W_ee | psi> - <psi | W_ee^{lr} |psi
 END_DOC
 psi_energy_wee_sr = psi_energy_two_e - psi_energy_erf
END_PROVIDER 
