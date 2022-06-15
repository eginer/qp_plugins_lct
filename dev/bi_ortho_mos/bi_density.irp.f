
! ---

BEGIN_PROVIDER [double precision, TCSCF_density_matrix_ao_alpha, (ao_num, ao_num) ]

  implicit none

  call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0               &
            , mo_r_coef, size(mo_r_coef, 1), mo_l_coef, size(mo_l_coef, 1) &
  !call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0               &
  !          , mo_l_coef, size(mo_l_coef, 1), mo_r_coef, size(mo_r_coef, 1) &
            , 0.d0, TCSCF_density_matrix_ao_alpha, size(TCSCF_density_matrix_ao_alpha, 1) )

  !integer          :: i, j
  !double precision :: trace_density
  !trace_density = 0.d0
  !do i = 1, ao_num
  !  do j = 1, ao_num
  !    trace_density = trace_density + TCSCF_density_matrix_ao_alpha(j,i) * ao_overlap(j,i)
  !  enddo
  !enddo
  !print *, ' trace of TCSCF_density_matrix_ao_alpha =', trace_density

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, TCSCF_density_matrix_ao_beta, (ao_num, ao_num) ]

  implicit none

  call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0                &
            , mo_r_coef, size(mo_r_coef, 1), mo_l_coef, size(mo_l_coef, 1) &
  !call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0                &
  !          , mo_l_coef, size(mo_l_coef, 1), mo_r_coef, size(mo_r_coef, 1) &
            , 0.d0, TCSCF_density_matrix_ao_beta, size(TCSCF_density_matrix_ao_beta, 1) )

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, TCSCF_density_matrix_ao, (ao_num, ao_num) ]

  implicit none

  ASSERT ( size(TCSCF_density_matrix_ao, 1) == size(TCSCF_density_matrix_ao_alpha, 1) )

  if( elec_alpha_num==elec_beta_num ) then
    TCSCF_density_matrix_ao = TCSCF_density_matrix_ao_alpha + TCSCF_density_matrix_ao_alpha
  else
    ASSERT ( size(TCSCF_density_matrix_ao, 1) == size(TCSCF_density_matrix_ao_beta, 1))
    TCSCF_density_matrix_ao = TCSCF_density_matrix_ao_alpha + TCSCF_density_matrix_ao_beta
  endif

END_PROVIDER

! ---

