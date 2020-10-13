program print_hole
 implicit none
 BEGIN_DOC
 ! This programs writes the effective RS-DFT Hamiltonian into the EZFIO folder. 
 ! The next programs that will run unto the EZFIO folder will, by default, 
 !
 ! have the one- and two-body integrals loaded from the EZFIO data. 
 END_DOC
 read_wf = .true.
 touch read_wf
 !! total one-e integrals 
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  
 !! Vne integrals on the MO basis 
 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 !! kinetic integrals on the MO basis 
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 !! Vne integrals on the AO basis 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 !! kinetic integrals on the AO basis 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 

 !! regular 1/r12 integrals  on the MO basis
 io_mo_two_e_integrals = "None"
 touch io_mo_two_e_integrals
 !! regular 1/r12 integrals  on the AO basis
 io_ao_two_e_integrals = "None"
 touch io_ao_two_e_integrals

 ! you read the right/left eigenvectors from EZFIO
 read_rl_eigv = .True.
 touch read_rl_eigv
 !! integral of the effective potential 
 call print_psi_exc_psi_trans
end
