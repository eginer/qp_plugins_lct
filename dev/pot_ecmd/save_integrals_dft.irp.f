
!ubroutine save_one_e_effective_potential  
!implicit none
!BEGIN_DOC 
! used to save the effective_one_e_potential into the one-body integrals in the ezfio folder
! this effective_one_e_potential is computed with the current density 
! and will couple the WFT with DFT for the next regular WFT calculation
!END_DOC
!call write_one_e_integrals('mo_ne_integral', effective_one_e_potential_without_kin, &
!     size(effective_one_e_potential_without_kin,1), size(effective_one_e_potential_without_kin,2))
!call write_one_e_integrals('mo_kinetic_integral',mo_kinetic_integral ,&
!     size(mo_kinetic_integral,1), size(mo_kinetic_integral,2))

!print *,  'Effective DFT potential is written on disk on the mo_ne_integral integrals'
!call ezfio_set_integrals_monoelec_disk_access_mo_one_integrals("Read")

!nd

subroutine save_one_e_effective_potential_ecmd_lda
 implicit none
 BEGIN_DOC 
! used to save the effective_one_e_potential into the one-body integrals in the ezfio folder
! this effective_one_e_potential is computed with the current density 
! and will couple the WFT with DFT for the next regular WFT calculation
 END_DOC
 call ezfio_set_mo_one_e_ints_mo_integrals_n_e(effective_one_e_potential_without_kin_ecmd_lda)
 call ezfio_set_mo_one_e_ints_mo_integrals_kinetic(mo_kinetic_integrals)
  
 print *,  'Effective DFT(Ec,md LDA) potential is written on disk on the mo_ne_integral integrals'
 call ezfio_set_mo_one_e_ints_io_mo_integrals_n_e("Read")
 
end


 subroutine save_one_e_effective_potential_ecmd_pbe_ueg
 implicit none
 BEGIN_DOC 
 !used to save the effective_one_e_potential into the one-body integrals in the ezfio folder
 !this effective_one_e_potential is computed with the current density 
 !and will couple the WFT with DFT for the next regular WFT calculation
 END_DOC
 call ezfio_set_mo_one_e_ints_mo_integrals_n_e(effective_one_e_potential_without_kin_ecmd_pbe_ueg)
 call ezfio_set_mo_one_e_ints_mo_integrals_kinetic(mo_kinetic_integrals)
  
 print *,  'Effective DFT(Ec,md PBE UEG) potential is written on disk on the mo_ne_integral integrals'
 call ezfio_set_mo_one_e_ints_io_mo_integrals_n_e("Read")
 
 end


!ubroutine save_erf_bi_elec_integrals_mo
!implicit none
!integer :: i,j,k,l
!PROVIDE mo_bielec_integrals_erf_in_map
!call ezfio_set_work_empty(.False.)
!call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_erf_map)
!call ezfio_set_integrals_bielec_disk_access_mo_integrals("Read")
!nd

!ubroutine save_erf_bi_elec_integrals_ao
!implicit none
!integer :: i,j,k,l
!PROVIDE ao_bielec_integrals_erf_in_map
!call ezfio_set_work_empty(.False.)
!call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_erf_map)
!call ezfio_set_integrals_bielec_disk_access_ao_integrals("Read")
!nd

!ubroutine write_all_integrals_for_mrdft
!implicit none
!call save_one_e_effective_potential
!call save_erf_bi_elec_integrals_ao
!call save_erf_bi_elec_integrals_mo

!nd

subroutine write_all_integrals_for_mrdft_ecmd_lda
 implicit none
 call save_one_e_effective_potential_ecmd_lda

end


subroutine write_all_integrals_for_mrdft_ecmd_pbe_ueg
 implicit none
 call save_one_e_effective_potential_ecmd_pbe_ueg

end

!ubroutine save_erf_mu_of_r_bi_elec_integrals_mo
!implicit none
!integer :: i,j,k,l
!PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map
!call ezfio_set_work_empty(.False.)
!call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_erf_mu_of_r_map)
!call ezfio_set_integrals_bielec_disk_access_mo_integrals("Read")
!nd

!ubroutine save_erf_mu_of_r_bi_elec_integrals_ao
!implicit none
!integer :: i,j,k,l
!PROVIDE ao_bielec_integrals_erf_mu_of_r_in_map
!call ezfio_set_work_empty(.False.)
!call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_erf_mu_of_r_map)
!call ezfio_set_integrals_bielec_disk_access_ao_integrals("Read")
!nd

