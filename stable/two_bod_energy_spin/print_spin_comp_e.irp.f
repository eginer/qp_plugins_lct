program two_bod_energy_spin
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
 call print_energy_comp

end

subroutine print_energy
 implicit none
 print*,'psi_energy_two_e_ab    = ',psi_energy_two_e_ab
 print*,'psi_energy_two_e_aa    = ',psi_energy_two_e_aa
 print*,'psi_energy_two_e_bb    = ',psi_energy_two_e_bb
 print*,'sum of biele components = ',psi_energy_two_e_ab + psi_energy_two_e_aa + psi_energy_two_e_bb
 print*,'psi_energy_two_e       = ',psi_energy_two_e
 print*,'psi_energy - <h_core>   = ',psi_energy - psi_energy_h_core

end


subroutine print_energy_comp
 implicit none
 print*,'psi_energy_bielec_ab    = ',psi_energy_two_e
 print*,'<h_core>                = ',psi_energy_h_core

end

