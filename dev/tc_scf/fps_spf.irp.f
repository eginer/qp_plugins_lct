BEGIN_PROVIDER [double precision, FPS_SPF_Matrix_AO_tc, (AO_num, AO_num)]
  implicit none
  BEGIN_DOC
  !   Commutator FPS - SPF
  END_DOC
  double precision, allocatable  :: scratch(:,:)
  allocate(                                                          &
      scratch(AO_num, AO_num)                                  &
      )

  ! Compute FP

  call dgemm('N','N',AO_num,AO_num,AO_num,                           &
      1.d0,                                                          &
      Fock_matrix_tc_ao_tot,Size(Fock_matrix_tc_ao_tot,1),                         &
      TCSCF_density_matrix_ao_tot,Size(TCSCF_density_matrix_ao_tot,1),             &
      0.d0,                                                          &
      scratch,Size(scratch,1))

  ! Compute FPS

  call dgemm('N','N',AO_num,AO_num,AO_num,                           &
      1.d0,                                                          &
      scratch,Size(scratch,1),                                       &
      AO_Overlap,Size(AO_Overlap,1),                                 &
      0.d0,                                                          &
      FPS_SPF_Matrix_AO_tc,Size(FPS_SPF_Matrix_AO_tc,1))

  ! Compute SP

  call dgemm('N','N',AO_num,AO_num,AO_num,                           &
      1.d0,                                                          &
      AO_Overlap,Size(AO_Overlap,1),                                 &
      TCSCF_density_matrix_ao_tot,Size(TCSCF_density_matrix_ao_tot,1),             &
      0.d0,                                                          &
      scratch,Size(scratch,1))

  ! Compute FPS - SPF

  call dgemm('N','N',AO_num,AO_num,AO_num,                           &
      -1.d0,                                                         &
      scratch,Size(scratch,1),                                       &
      Fock_matrix_tc_ao_tot,Size(Fock_matrix_tc_ao_tot,1),                         &
      1.d0,                                                          &
      FPS_SPF_Matrix_AO_tc,Size(FPS_SPF_Matrix_AO_tc,1))

END_PROVIDER

BEGIN_PROVIDER [double precision, FPS_SPF_Matrix_MO_tc, (mo_num, mo_num)]
  implicit none
  begin_doc
!   Commutator FPS - SPF in MO basis
  end_doc
  call ao_to_mo_bi_ortho(FPS_SPF_Matrix_AO_tc, size(FPS_SPF_Matrix_AO_tc,1), &
     FPS_SPF_Matrix_MO_tc, size(FPS_SPF_Matrix_MO_tc,1))
END_PROVIDER


