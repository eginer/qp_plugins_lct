program two_bod_energy_spin
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
 call print_energy

end

subroutine print_energy
 implicit none
 print*,'psi_energy_bielec_ab    = ',psi_energy_bielec_ab
 print*,'psi_energy_bielec_aa    = ',psi_energy_bielec_aa
 print*,'psi_energy_bielec_bb    = ',psi_energy_bielec_bb
 print*,'sum of biele components = ',psi_energy_bielec_ab + psi_energy_bielec_aa + psi_energy_bielec_bb
 print*,'psi_energy_bielec       = ',psi_energy_bielec
 print*,'psi_energy - <h_core>   = ',psi_energy - psi_energy_h_core

end

