subroutine save_one_e_tot_eff_pot_basis_electric_field
 implicit none
 BEGIN_DOC 
! used to save the effective_one_e_potential into the one-body integrals in the ezfio folder
! this effective_one_e_potential is computed with the current density 
! and will couple the WFT with DFT for the next regular WFT calculation
 END_DOC

 call ezfio_set_mo_one_e_ints_mo_integrals_n_e(mo_one_e_integrals_electric_field) ! mo_integrals_n_e + mo_kinetic_integrals + epsilon_<mu_z>_mo 
 call ezfio_set_mo_one_e_ints_mo_integrals_kinetic(- mo_kinetic_integrals)

 call ezfio_set_ao_one_e_ints_ao_integrals_n_e(ao_one_e_integrals_electric_field) ! ao_integrals_n_e + ao_kinetic_integrals + epsilon_<mu_z>_ao 
 call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(- ao_kinetic_integrals)

 call ezfio_set_mo_one_e_ints_io_mo_integrals_n_e("Read")
 call ezfio_set_mo_one_e_ints_io_mo_integrals_kinetic("Read")
 call ezfio_set_ao_one_e_ints_io_ao_integrals_n_e("Read")
 call ezfio_set_ao_one_e_ints_io_ao_integrals_kinetic("Read")

 print *,  'Effective DFT potential is written on disk on the one electron integrals on the AO/MO basis'
                              

end

!subroutine save_eff_basis_two_e_ints
!implicit none
!integer :: i,j,k,l
!PROVIDE mo_two_e_int_mu_of_r_in_map
!print*,'Writting effective integrals in EZFIO'
!call ezfio_set_work_empty(.False.)
!call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_int_mu_of_r_map)
!call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')
!end

!subroutine save_regular_two_e_ints_mo
!implicit none
!integer :: i,j,k,l
!PROVIDE mo_two_e_integrals_in_map
!call ezfio_set_work_empty(.False.)
!call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_map)
!call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')
!end


subroutine write_all_ints_basis_electric_field
 implicit none
 BEGIN_DOC
 ! saves all integrals needed for RS-DFT-MRCI calculation: 
 !
 ! one-body effective potential and two-elec erf integrals
 END_DOC
 call save_one_e_tot_eff_pot_basis_electric_field
 call save_regular_two_e_ints_mo
end


