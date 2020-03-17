program slater_rules_erf
  implicit none
  BEGIN_DOC
! TODO : Small example to use the different quantities in this plugin
  END_DOC

  !! You force QP2 to read the wave function in the EZFIO folder
  read_wf = .True.
  touch read_wf 
  call routine_slater_rules_erf
  
end

subroutine routine_slater_rules_erf
 implicit none
 print*,'psi_energy_erf = ',psi_energy_erf
end
