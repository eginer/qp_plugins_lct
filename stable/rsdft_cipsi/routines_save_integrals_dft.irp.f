
subroutine save_one_e_effective_potential  
 implicit none
 BEGIN_DOC 
! used to save the effective_one_e_potential into the one-body integrals in the ezfio folder
! this effective_one_e_potential is computed with the current density 
! and will couple the WFT with DFT for the next regular WFT calculation
 END_DOC

 call ezfio_set_mo_one_e_ints_mo_integrals_n_e(effective_one_e_potential_without_kin)
 call ezfio_set_mo_one_e_ints_mo_integrals_kinetic(mo_kinetic_integrals)

 call ezfio_set_ao_one_e_ints_ao_integrals_n_e(ao_effective_one_e_potential_without_kin)
 call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(ao_kinetic_integrals)

 call ezfio_set_mo_one_e_ints_io_mo_integrals_n_e("Read")
 call ezfio_set_mo_one_e_ints_io_mo_integrals_kinetic("Read")
 call ezfio_set_ao_one_e_ints_io_ao_integrals_n_e("Read")
 call ezfio_set_ao_one_e_ints_io_ao_integrals_kinetic("Read")

 print *,  'Effective DFT potential is written on disk on the one electron integrals on the AO/MO basis'
                              

end

subroutine write_all_integrals_for_mrdft
 implicit none
 BEGIN_DOC
 ! saves all integrals needed for RS-DFT-MRCI calculation: 
 !
 ! one-body effective potential and two-elec erf integrals
 END_DOC
 PROVIDE mo_two_e_integrals_erf_in_map mo_two_e_integrals_in_map
 call save_one_e_effective_potential
 call save_erf_two_e_ints_mo_into_ints_mo
 call save_erf_two_e_ints_ao_into_ints_ao
end


subroutine write_all_integrals_for_mrdft_read_all
 implicit none
 BEGIN_DOC
 ! saves all integrals needed for RS-DFT-MRCI calculation: 
 !
 ! one-body effective potential and two-elec erf integrals
 END_DOC
 call save_one_e_effective_potential
end

