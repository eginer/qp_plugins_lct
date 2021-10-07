subroutine save_erf_mu_of_r_mo_into_erf_ints_mo
 implicit none
 integer :: i,j,k,l
 PROVIDE mo_int_erf_mu_of_r_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_erf',mo_int_erf_mu_of_r_map)
 call ezfio_set_mo_two_e_erf_ints_io_mo_two_e_integrals_erf('Read')
end


subroutine save_erf_two_e_ints_mo_into_ints_mo
 implicit none
 integer :: i,j,k,l
 PROVIDE mo_int_erf_mu_of_r_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_int_erf_mu_of_r_map)
 call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')
end

subroutine save_erf_mu_of_r_ao_into_erf_ints_ao
 implicit none
 integer :: i,j,k,l
 PROVIDE ao_int_erf_mu_of_r_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_erf',ao_int_erf_mu_of_r_map)
 call ezfio_set_ao_two_e_erf_ints_io_ao_two_e_integrals_erf('Read')
end


subroutine save_erf_two_e_ints_ao_into_ints_ao
 implicit none
 integer :: i,j,k,l
 PROVIDE ao_int_erf_mu_of_r_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_int_erf_mu_of_r_map)
 call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals('Read')
end

