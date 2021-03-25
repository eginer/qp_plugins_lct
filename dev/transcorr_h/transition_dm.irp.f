
 BEGIN_PROVIDER [ double precision, one_e_tm_mo_alpha, (mo_num,mo_num,N_states) ]
&BEGIN_PROVIDER [ double precision, one_e_tm_mo_beta,  (mo_num,mo_num,N_states) ]
  implicit none
  BEGIN_DOC
  ! $\alpha$ and $\beta$ one-body density matrix for each state
  END_DOC

  integer                        :: j,k,l,m1,m2,k_a,k_b
  integer                        :: occ(N_int*bit_kind_size,2)
  double precision               :: ck, cl, ckl
  double precision               :: phase
  integer                        :: h1,h2,p1,p2,s1,s2, degree
  integer(bit_kind)              :: tmp_det(N_int,2), tmp_det2(N_int)
  integer                        :: exc(0:2,2),n_occ(2)
  double precision, allocatable  :: tmp_a(:,:,:), tmp_b(:,:,:)
  integer                        :: krow, kcol, lrow, lcol

  PROVIDE psi_det

  one_e_tm_mo_alpha = 0.d0
  one_e_tm_mo_beta  = 0.d0
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(j,k,k_a,k_b,l,m1,m2,occ,ck, cl, ckl,phase,h1,h2,p1,p2,s1,s2, degree,exc,&
      !$OMP  tmp_a, tmp_b, n_occ, krow, kcol, lrow, lcol, tmp_det, tmp_det2)&
      !$OMP SHARED(psi_det,psi_coef,N_int,N_states,elec_alpha_num,  &
      !$OMP  elec_beta_num,one_e_tm_mo_alpha,one_e_tm_mo_beta,N_det,&
      !$OMP  mo_num,psi_bilinear_matrix_rows,psi_bilinear_matrix_columns,&
      !$OMP  psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns,&
      !$OMP  psi_bilinear_matrix_order_reverse, psi_det_alpha_unique, psi_det_beta_unique,&
      !$OMP  psi_bilinear_matrix_values, psi_bilinear_matrix_transp_values,&
      !$OMP  N_det_alpha_unique,N_det_beta_unique,irp_here)
  allocate(tmp_a(mo_num,mo_num,N_states), tmp_b(mo_num,mo_num,N_states) )
  tmp_a = 0.d0
  !$OMP DO SCHEDULE(dynamic,64)
  do k_a=1,N_det
    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
    tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)

    ! Diagonal part
    ! -------------

    call bitstring_to_list_ab(tmp_det, occ, n_occ, N_int)
!    do m=1,N_states
     m1 = 1
     m2 = 2
      ck = psi_bilinear_matrix_values(k_a,m1)*psi_bilinear_matrix_values(k_a,m2)
      do l=1,elec_alpha_num
        j = occ(l,1)
        tmp_a(j,j,m1) += ck
      enddo
!    enddo

    if (k_a == N_det) cycle
    l = k_a+1
    lrow = psi_bilinear_matrix_rows(l)
    lcol = psi_bilinear_matrix_columns(l)
    ! Fix beta determinant, loop over alphas
    do while ( lcol == kcol )
      tmp_det2(:) = psi_det_alpha_unique(:, lrow)
      call get_excitation_degree_spin(tmp_det(1,1),tmp_det2,degree,N_int)
      if (degree == 1) then
        exc = 0
        call get_single_excitation_spin(tmp_det(1,1),tmp_det2,exc,phase,N_int)
        call decode_exc_spin(exc,h1,p1,h2,p2)
!        do m=1,N_states
        m1 = 1
        m2 = 2
          ckl = psi_bilinear_matrix_values(k_a,m1)*psi_bilinear_matrix_values(l,m2) * phase
          tmp_a(h1,p1,m1) += ckl
          tmp_a(p1,h1,m1) += ckl
!        enddo
      endif
      l = l+1
      if (l>N_det) exit
      lrow = psi_bilinear_matrix_rows(l)
      lcol = psi_bilinear_matrix_columns(l)
    enddo

  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  one_e_tm_mo_alpha(:,:,:) = one_e_tm_mo_alpha(:,:,:) + tmp_a(:,:,:)
  !$OMP END CRITICAL
  deallocate(tmp_a)

  tmp_b = 0.d0
  !$OMP DO SCHEDULE(dynamic,64)
  do k_b=1,N_det
    krow = psi_bilinear_matrix_transp_rows(k_b)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_transp_columns(k_b)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
    tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)

    ! Diagonal part
    ! -------------

    call bitstring_to_list_ab(tmp_det, occ, n_occ, N_int)
!    do m=1,N_states
    m1 = 1
    m2 = 2
      ck = psi_bilinear_matrix_transp_values(k_b,m1)*psi_bilinear_matrix_transp_values(k_b,m2)
      do l=1,elec_beta_num
        j = occ(l,2)
        tmp_b(j,j,m1) += ck
      enddo
!    enddo

    if (k_b == N_det) cycle
    l = k_b+1
    lrow = psi_bilinear_matrix_transp_rows(l)
    lcol = psi_bilinear_matrix_transp_columns(l)
    ! Fix beta determinant, loop over alphas
    do while ( lrow == krow )
      tmp_det2(:) = psi_det_beta_unique(:, lcol)
      call get_excitation_degree_spin(tmp_det(1,2),tmp_det2,degree,N_int)
      if (degree == 1) then
        exc = 0
        call get_single_excitation_spin(tmp_det(1,2),tmp_det2,exc,phase,N_int)
        call decode_exc_spin(exc,h1,p1,h2,p2)
!        do m=1,N_states
        m1 = 1
        m2 = 2
          ckl = psi_bilinear_matrix_transp_values(k_b,m1)*psi_bilinear_matrix_transp_values(l,m2) * phase
          tmp_b(h1,p1,m1) += ckl
          tmp_b(p1,h1,m1) += ckl
!        enddo
      endif
      l = l+1
      if (l>N_det) exit
      lrow = psi_bilinear_matrix_transp_rows(l)
      lcol = psi_bilinear_matrix_transp_columns(l)
    enddo

  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
  one_e_tm_mo_beta(:,:,:)  = one_e_tm_mo_beta(:,:,:)  + tmp_b(:,:,:)
  !$OMP END CRITICAL

  deallocate(tmp_b)
  !$OMP END PARALLEL

END_PROVIDER

 BEGIN_PROVIDER [ double precision, one_e_tm_mo, (mo_num,mo_num,N_states) ]
&BEGIN_PROVIDER [ double precision, one_e_tm_mo_norm ]
 implicit none
 integer :: i
 one_e_tm_mo_norm = 0.d0
 one_e_tm_mo =  one_e_tm_mo_beta + one_e_tm_mo_alpha
 do i = 1, mo_num
  one_e_tm_mo_norm += one_e_tm_mo(i,i,1)
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, left_right_overlap_read]
&BEGIN_PROVIDER [ double precision, left_overlap_read]
&BEGIN_PROVIDER [ double precision, right_overlap_read]
 implicit none
 integer :: i
 left_right_overlap_read = 0.d0
 do i = 1, N_det
  left_right_overlap_read += psi_coef(i,1) * psi_coef(i,2)
  right_overlap_read += psi_coef(i,1)**2.D0
  left_overlap_read += psi_coef(i,2)**2.D0
 enddo

 END_PROVIDER 
