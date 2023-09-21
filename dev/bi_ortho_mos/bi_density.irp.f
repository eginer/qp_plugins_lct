
! ---

BEGIN_PROVIDER [double precision, TCSCF_bi_ort_dm_ao_alpha, (ao_num, ao_num) ]
  implicit none
  call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0               &
            , mo_l_coef, size(mo_l_coef, 1), mo_r_coef, size(mo_r_coef, 1) &
            , 0.d0, TCSCF_bi_ort_dm_ao_alpha, size(TCSCF_bi_ort_dm_ao_alpha, 1) )
  integer :: i
! print*,'mo_l_coef ==='
! do i = 1, mo_num
!  write(*,'(100(F16.10,X))')mo_l_coef(:,i)
! enddo
 
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, TCSCF_bi_ort_dm_ao_beta, (ao_num, ao_num) ]

  implicit none

  call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0               &
            , mo_l_coef, size(mo_l_coef, 1), mo_r_coef, size(mo_r_coef, 1) &
            , 0.d0, TCSCF_bi_ort_dm_ao_beta, size(TCSCF_bi_ort_dm_ao_beta, 1) )
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, TCSCF_bi_ort_dm_ao, (ao_num, ao_num) ]

  implicit none

  ASSERT ( size(TCSCF_bi_ort_dm_ao, 1) == size(TCSCF_bi_ort_dm_ao_alpha, 1) )

  if( elec_alpha_num==elec_beta_num ) then
    TCSCF_bi_ort_dm_ao = TCSCF_bi_ort_dm_ao_alpha + TCSCF_bi_ort_dm_ao_alpha
  else
    ASSERT ( size(TCSCF_bi_ort_dm_ao, 1) == size(TCSCF_bi_ort_dm_ao_beta, 1))
    TCSCF_bi_ort_dm_ao = TCSCF_bi_ort_dm_ao_alpha + TCSCF_bi_ort_dm_ao_beta
  endif

END_PROVIDER

! ---

