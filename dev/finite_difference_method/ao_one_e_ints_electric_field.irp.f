 BEGIN_PROVIDER [ double precision, ao_one_e_integrals_electric_field,(ao_num,ao_num)]
&BEGIN_PROVIDER [ double precision, ao_one_e_integrals_diag_electric_field,(ao_num)]
  implicit none
  integer :: i,j,n,l
  BEGIN_DOC
 ! One-electron Hamiltonian in the |AO| basis.
  END_DOC

  IF (read_ao_one_e_integrals) THEN
     call ezfio_get_ao_one_e_ints_ao_one_e_integrals(ao_one_e_integrals_electric_field)
  ELSE
        ao_one_e_integrals_electric_field = ao_integrals_n_e + ao_kinetic_integrals - field_strenght*ao_dipole_z

  ENDIF

  DO j = 1, ao_num
    ao_one_e_integrals_diag_electric_field(j) = ao_one_e_integrals_electric_field(j,j)
  ENDDO

  IF (write_ao_one_e_integrals) THEN
       call ezfio_set_ao_one_e_ints_ao_one_e_integrals(ao_one_e_integrals_electric_field)
       print *,  'AO one-e integrals written to disk'
  ENDIF

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_one_e_integrals_imag_electric_field,(ao_num,ao_num)]
  implicit none
  integer :: i,j,n,l
  BEGIN_DOC
 ! One-electron Hamiltonian in the |AO| basis.
  END_DOC

  IF (read_ao_one_e_integrals) THEN
     call ezfio_get_ao_one_e_ints_ao_one_e_integrals(ao_one_e_integrals_imag_electric_field)
  ELSE
     print *,  irp_here, ': Not yet implemented'
     stop -1
  ENDIF

  IF (write_ao_one_e_integrals) THEN
       call ezfio_set_ao_one_e_ints_ao_one_e_integrals(ao_one_e_integrals_imag_electric_field)
       print *,  'AO one-e integrals written to disk'
  ENDIF

END_PROVIDER

