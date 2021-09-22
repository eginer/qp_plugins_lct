program write_effective_z_electric_field_hamiltonian
 implicit none
 BEGIN_DOC
 !
 ! The next programs that will run unto the EZFIO folder will, by default, 
 !
 ! have the one- and two-body integrals loaded from the EZFIO data. 
 END_DOC
 ! total one-e integrals 
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  
 ! Vne integrals on the MO basis 
 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 ! kinetic integrals on the MO basis 
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 ! Vne integrals on the AO basis 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 ! kinetic integrals on the AO basis 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 
 io_ao_one_e_integrals = "None"
 touch io_ao_one_e_integrals  

 call save_v_ne_and_z_electric_filed_ao_ints
end
