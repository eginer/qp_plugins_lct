program DFT_Utils_ECMD
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
 !print*,'ecmd_pbe_ueg_prov = ',ecmd_pbe_ueg_prov
  read_wf = .True.
  touch read_wf

  io_mo_one_e_integrals = "None"
  touch io_mo_one_e_integrals  
  io_mo_two_e_integrals = "None"
  touch io_mo_two_e_integrals
  io_ao_two_e_integrals = "None"
  touch io_ao_two_e_integrals
  io_mo_two_e_integrals_erf = "None" 
  touch io_mo_two_e_integrals_erf
  io_ao_two_e_integrals_erf = "None" 
  touch io_ao_two_e_integrals_erf
 
  io_mo_integrals_e_n = "None"
  touch io_mo_integrals_e_n
  io_mo_integrals_kinetic = "None"
  touch io_mo_integrals_kinetic 
  io_ao_integrals_e_n = "None"
  touch io_ao_integrals_e_n 
  io_ao_integrals_kinetic = "None"
  touch io_ao_integrals_kinetic 

  call print_ecmd_var_energy
 !call print_ecmd_var_energy_barth
end


