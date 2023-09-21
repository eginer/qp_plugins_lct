program save_mu_of_r_ints
 implicit none
 read_wf = .True.
 touch read_wf 
 
 io_mo_two_e_integrals_erf = "None" 
 touch io_mo_two_e_integrals_erf
 io_ao_two_e_integrals_erf = "None" 
 touch io_ao_two_e_integrals_erf
 call routine
end

subroutine routine
 implicit none
 call save_erf_mu_of_r_ao_into_erf_ints_ao
 call save_erf_mu_of_r_mo_into_erf_ints_mo
 call ezfio_set_mu_of_r_mu_of_r_disk(mu_of_r_dft)
 call ezfio_set_dft_keywords_mu_dft_type("Read")
end
