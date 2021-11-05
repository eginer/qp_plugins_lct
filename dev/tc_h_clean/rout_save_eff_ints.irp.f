subroutine save_eff_tc_two_e_ints_mo_into_ints_mo
 implicit none
 integer :: i,j,k,l
 PROVIDE mo_two_e_integrals_tc_int_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_tc_int_map)
 call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')
! call ezfio_set_tc_eff_map_zero_tc_eff_map('True')
end

