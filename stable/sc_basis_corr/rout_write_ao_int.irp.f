
subroutine save_one_e_tot_ao_eff_pot_basis
 implicit none
 BEGIN_DOC 
! used to save the effective_one_e_potential into the one-body integrals in the ezfio folder
! this effective_one_e_potential is computed with the current density 
! and will couple the WFT with DFT for the next regular WFT calculation
 END_DOC

 call ezfio_set_ao_one_e_ints_ao_integrals_n_e(ao_tot_eff_pot_basis_no_kin)
 call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(ao_kinetic_integrals)

 call ezfio_set_ao_one_e_ints_io_ao_integrals_n_e("Read")
 call ezfio_set_ao_one_e_ints_io_ao_integrals_kinetic("Read")

 print *,  'Effective DFT potential is written on disk on the one electron integrals on the AO basis'
                              

end

subroutine save_regular_two_e_ints_ao
 implicit none
 integer :: i,j,k,l
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
 call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals('Read')
end


subroutine write_all_ints_basis_ao
 implicit none
 BEGIN_DOC
 ! saves all integrals needed for RS-DFT-MRCI calculation: 
 !
 ! one-body effective potential and two-elec erf integrals
 END_DOC
 call save_one_e_tot_ao_eff_pot_basis
! print*,'needs_eff_two_e_ints = ',needs_eff_two_e_ints
! pause
 if(needs_eff_two_e_ints)then
  call save_eff_basis_ao_two_e_ints
 else 
  call save_regular_two_e_ints_ao
 endif
end


subroutine save_eff_basis_ao_two_e_ints
 implicit none
 integer :: i,j,k,l
 PROVIDE ao_int_mu_of_r_in_map
 print*,'Writting effective integrals in EZFIO in the AO basis'
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_int_mu_of_r_map)
 call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals('Read')
end
