
 BEGIN_PROVIDER [double precision, mo_bi_orth_bipole_x , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_bi_orth_bipole_y , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_bi_orth_bipole_z , (mo_num,mo_num)]
 BEGIN_DOC
 ! array of the integrals of MO_i * x MO_j
 ! array of the integrals of MO_i * y MO_j
 ! array of the integrals of MO_i * z MO_j
 END_DOC
 implicit none

  call ao_to_mo_bi_ortho(                                                     &
      ao_dipole_x,                                                   &
      size(ao_dipole_x,1),                                           &
      mo_bi_orth_bipole_x,                                                   &
      size(mo_bi_orth_bipole_x,1)                                            &
      )
  call ao_to_mo_bi_ortho(                                                     &
      ao_dipole_y,                                                   &
      size(ao_dipole_y,1),                                           &
      mo_bi_orth_bipole_y,                                                   &
      size(mo_bi_orth_bipole_y,1)                                            &
      )
  call ao_to_mo_bi_ortho(                                                     &
      ao_dipole_z,                                                   &
      size(ao_dipole_z,1),                                           &
      mo_bi_orth_bipole_z,                                                   &
      size(mo_bi_orth_bipole_z,1)                                            &
      )

END_PROVIDER
