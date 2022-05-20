
BEGIN_PROVIDER [double precision, FQS_SQF_Matrix_AO, (ao_num,ao_num)]

  BEGIN_DOC
  ! the error is defined as: FQS - SQF
  END_DOC

  implicit none
  double precision, allocatable :: scratch(:,:)

  allocate( scratch(ao_num,ao_num) )

  ! compute FQ
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, 1.d0                 &
            , Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot,1)   &
            , Q_ao, size(Q_ao,1), 0.d0, scratch, size(scratch,1) )

  ! compute FQS
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, 1.d0                   &
            , scratch, size(scratch,1), AO_Overlap, size(AO_Overlap,1) &
            , 0.d0, FQS_SQF_Matrix_AO, size(FQS_SQF_Matrix_AO,1) )

  ! compute SQ
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, 1.d0             &
            , AO_Overlap, size(AO_Overlap,1), Q_ao, size(Q_ao,1) &
            , 0.d0, scratch, size(scratch,1) )

  ! compute FQS - SQF
  call dgemm( 'N', 'N', ao_num, ao_num, ao_num, -1.d0              &
            , scratch, size(scratch,1)                             &
            , Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot,1) &
            , 1.d0, FQS_SQF_Matrix_AO, size(FQS_SQF_Matrix_AO,1) )

  deallocate( scratch )

END_PROVIDER



BEGIN_PROVIDER [double precision, FQS_SQF_Matrix_MO, (mo_num, mo_num)]

  BEGIN_DOC
  ! commutator FQS - SQF in MO basis
  END_DOC

  call ao_to_mo( FQS_SQF_Matrix_AO, size(FQS_SQF_Matrix_AO,1) &
               , FQS_SQF_Matrix_MO, size(FQS_SQF_Matrix_MO,1) )

END_PROVIDER


