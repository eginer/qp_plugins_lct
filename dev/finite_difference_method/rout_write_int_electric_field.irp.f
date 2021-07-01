subroutine save_v_ne_and_x_electric_filed_ao_ints
 implicit none
 BEGIN_DOC 

 END_DOC

 call ezfio_set_ao_one_e_ints_ao_integrals_n_e(ao_one_e_potential_n_e_and_electric_field_x) ! ao_integrals_n_e - epsilon_<mu_x>_ao 
 call ezfio_set_ao_one_e_ints_io_ao_integrals_n_e("Read")
end

subroutine save_v_ne_and_y_electric_filed_ao_ints
 implicit none
 BEGIN_DOC 

 END_DOC

 call ezfio_set_ao_one_e_ints_ao_integrals_n_e(ao_one_e_potential_n_e_and_electric_field_y) ! ao_integrals_n_e - epsilon_<mu_y>_ao 
 call ezfio_set_ao_one_e_ints_io_ao_integrals_n_e("Read")
end

subroutine save_v_ne_and_z_electric_filed_ao_ints
 implicit none
 BEGIN_DOC 

 END_DOC

 call ezfio_set_ao_one_e_ints_ao_integrals_n_e(ao_one_e_potential_n_e_and_electric_field_z) ! ao_integrals_n_e - epsilon_<mu_z>_ao 
 call ezfio_set_ao_one_e_ints_io_ao_integrals_n_e("Read")
end

 BEGIN_PROVIDER [ double precision, ao_one_e_integrals_tot_electric_field_x,(ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, ao_one_e_integrals_tot_electric_field_y,(ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, ao_one_e_integrals_tot_electric_field_z,(ao_num,ao_num)]
&BEGIN_PROVIDER [ double precision, ao_one_e_potential_n_e_and_electric_field_x,(ao_num,ao_num)]
&BEGIN_PROVIDER [ double precision, ao_one_e_potential_n_e_and_electric_field_y,(ao_num,ao_num)]
&BEGIN_PROVIDER [ double precision, ao_one_e_potential_n_e_and_electric_field_z,(ao_num,ao_num)]
  implicit none
  BEGIN_DOC
 ! One-electron Hamiltonian in the |AO| basis.
  END_DOC

  ao_one_e_integrals_tot_electric_field_x = ao_integrals_n_e + ao_kinetic_integrals - field_strenght*ao_dipole_total_x
  ao_one_e_integrals_tot_electric_field_y = ao_integrals_n_e + ao_kinetic_integrals - field_strenght*ao_dipole_total_y
  ao_one_e_integrals_tot_electric_field_z = ao_integrals_n_e + ao_kinetic_integrals - field_strenght*ao_dipole_total_z
  
  ao_one_e_potential_n_e_and_electric_field_x = ao_integrals_n_e - field_strenght*ao_dipole_total_x
  ao_one_e_potential_n_e_and_electric_field_y = ao_integrals_n_e - field_strenght*ao_dipole_total_y
  ao_one_e_potential_n_e_and_electric_field_z = ao_integrals_n_e - field_strenght*ao_dipole_total_z
END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_dipole_total_x, (ao_num, ao_num)]
&BEGIN_PROVIDER [ double precision, ao_dipole_total_y, (ao_num, ao_num)]
&BEGIN_PROVIDER [ double precision, ao_dipole_total_z, (ao_num, ao_num)]
 implicit none
 integer :: i,j
 double precision :: nuclei_part 
 nuclei_part = 0.d0
 do i = 1,nucl_num
  nuclei_part += nucl_charge(i) * nucl_coord(i,3)
 enddo
 print*,'nuclei_part = ',nuclei_part
 do i = 1, ao_num
   do j = 1, ao_num
!   ao_dipole_total_z(j,i) = ao_overlap(j,i) * nuclei_part - ao_dipole_z(j,i)
    ao_dipole_total_x(j,i) =  - ao_dipole_x(j,i)
    ao_dipole_total_y(j,i) =  - ao_dipole_y(j,i)
    ao_dipole_total_z(j,i) =  - ao_dipole_z(j,i)
   enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [ double precision, mo_one_e_integrals_electric_field_x,(mo_num,mo_num)]
&BEGIN_PROVIDER [ double precision, mo_one_e_integrals_electric_field_y,(mo_num,mo_num)]
&BEGIN_PROVIDER [ double precision, mo_one_e_integrals_electric_field_z,(mo_num,mo_num)]
  implicit none
  integer                        :: i,j,n,l
  BEGIN_DOC
  ! array of the one-electron Hamiltonian on the |MO| basis :
  ! sum of the kinetic and nuclear electronic potentials (and pseudo potential if needed)
  END_DOC
  mo_one_e_integrals_electric_field_x  = mo_integrals_n_e + mo_kinetic_integrals - field_strenght*mo_dipole_x
  mo_one_e_integrals_electric_field_y  = mo_integrals_n_e + mo_kinetic_integrals - field_strenght*mo_dipole_y
  mo_one_e_integrals_electric_field_z  = mo_integrals_n_e + mo_kinetic_integrals - field_strenght*mo_dipole_z
END_PROVIDER
