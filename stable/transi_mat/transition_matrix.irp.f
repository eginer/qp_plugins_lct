 BEGIN_PROVIDER [ double precision, transition_matrix_transposed, (N_states,N_states,mo_num, mo_num) ]
  implicit none
  BEGIN_DOC
  ! Transition density matrix 
  END_DOC

  integer                        :: j,k,l,m,k_a,k_b,n
  integer                        :: occ(N_int*bit_kind_size,2)
  double precision               :: ck, cl, ckl
  double precision               :: phase
  integer                        :: h1,h2,p1,p2,s1,s2, degree
  integer(bit_kind)              :: tmp_det(N_int,2), tmp_det2(N_int)
  integer                        :: exc(0:2,2),n_occ(2)
  double precision, allocatable  :: transi_mat_tmp(:,:,:,:)
  integer                        :: krow, kcol, lrow, lcol

  PROVIDE psi_det mo_dipole_x mo_dipole_y mo_dipole_z 

  print*,'providing the transition_matrix_transposed  '
  transition_matrix_transposed = 0.d0
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(j,k,n,k_a,k_b,l,m,occ,ck, cl, ckl,phase,h1,h2,p1,p2,s1,s2, degree,exc,&
      !$OMP   transi_mat_tmp, n_occ, krow, kcol, lrow, lcol, tmp_det, tmp_det2)&
      !$OMP SHARED(psi_det,psi_coef,N_int,N_states,elec_alpha_num,mo_num, &
      !$OMP  elec_beta_num,transition_matrix_transposed,N_det,&
      !$OMP  psi_bilinear_matrix_rows,psi_bilinear_matrix_columns,&
      !$OMP  psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns,&
      !$OMP  psi_bilinear_matrix_order_reverse, psi_det_alpha_unique, psi_det_beta_unique,&
      !$OMP  psi_bilinear_matrix_values, psi_bilinear_matrix_transp_values,&
      !$OMP  N_det_alpha_unique,N_det_beta_unique,irp_here)
  allocate(transi_mat_tmp(N_states,N_states,mo_num, mo_num))
  transi_mat_tmp = 0.d0
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
    do m=1,N_states
      do n = 1,N_states
       ck = psi_bilinear_matrix_values(k_a,m)*psi_bilinear_matrix_values(k_a,n)
       do l=1,elec_alpha_num
         j = occ(l,1)
         transi_mat_tmp(n,m,j,j) += ck 
       enddo
      enddo
    enddo

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
        do m=1,N_states
         do n = 1,N_states
          ckl = psi_bilinear_matrix_values(k_a,m)*psi_bilinear_matrix_values(l,n) * phase
          transi_mat_tmp(n,m,h1,p1) +=  ckl 
          ckl = psi_bilinear_matrix_values(k_a,n)*psi_bilinear_matrix_values(l,m) * phase
          transi_mat_tmp(n,m,h1,p1) +=  ckl 
         enddo
        enddo
      endif
      l = l+1
      if (l>N_det) exit
      lrow = psi_bilinear_matrix_rows(l)
      lcol = psi_bilinear_matrix_columns(l)
    enddo

  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
   transition_matrix_transposed += transi_mat_tmp
  !$OMP END CRITICAL
  deallocate(transi_mat_tmp)
  allocate(transi_mat_tmp(N_states,N_states, mo_num, mo_num))
  transi_mat_tmp = 0.d0

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
    do m=1,N_states
      do n = 1,N_states
       ck = psi_bilinear_matrix_values(k_b,m)*psi_bilinear_matrix_values(k_b,n)
       do l=1,elec_beta_num
         j = occ(l,2)
         transi_mat_tmp(n,m,j,j) += ck 
       enddo
      enddo
    enddo

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
        do m=1,N_states
         do n = 1,N_states
          ckl = psi_bilinear_matrix_transp_values(k_b,m)*psi_bilinear_matrix_transp_values(l,n) * phase
          transi_mat_tmp(n,m,h1,p1) +=  ckl 
          ckl = psi_bilinear_matrix_transp_values(k_b,n)*psi_bilinear_matrix_transp_values(l,m) * phase
          transi_mat_tmp(n,m,h1,p1) +=  ckl 
         enddo
        enddo
      endif
      l = l+1
      if (l>N_det) exit
      lrow = psi_bilinear_matrix_transp_rows(l)
      lcol = psi_bilinear_matrix_transp_columns(l)
    enddo

  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
   transition_matrix_transposed += transi_mat_tmp
  !$OMP END CRITICAL

  deallocate(transi_mat_tmp)
  !$OMP END PARALLEL
  print*,'provided the transition_matrix_transposed '

END_PROVIDER

 BEGIN_PROVIDER [ double precision, transition_matrix, (mo_num, mo_num,N_states,N_states) ]
 implicit none
 integer :: i,j,istate,jstate,m,n,p,h
! do i = 1, mo_num
!  do j = 1, mo_num
!   do istate = 1, N_states
!    do jstate = 1, N_states
!     transition_matrix(j,i,jstate,istate) = transition_matrix_transposed(jstate,istate,j,i)
!    enddo
!   enddo
!  enddo
! enddo
 double precision :: phase
 integer, allocatable           :: occ(:,:)
 integer                        :: n_occ_ab(2),degree,exc(0:2,2,2)
 allocate(occ(N_int*bit_kind_size,2))
 transition_matrix = 0.d0
 do istate = 1, N_states
  do jstate = 1, N_states
   do i = 1, N_det
    do j = 1, N_det
     call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
     if(degree.gt.1)then
      cycle
     else if (degree == 0)then
      call bitstring_to_list_ab(psi_det(1,1,i), occ, n_occ_ab, N_int)
      do p = 1, n_occ_ab(1) ! browsing the alpha electrons
       m = occ(p,1)
       transition_matrix(m,m,istate,jstate)+= psi_coef(i,istate) * psi_coef(j,jstate)
      enddo
      do p = 1, n_occ_ab(2) ! browsing the beta electrons
       m = occ(p,1)
       transition_matrix(m,m,istate,jstate)+= psi_coef(i,istate) * psi_coef(j,jstate)
      enddo
     else
      call get_single_excitation(psi_det(1,1,j),psi_det(1,1,i),exc,phase,N_int)
!      call debug_det(psi_det(1,1,j),N_int)
!      call debug_det(psi_det(1,1,i),N_int)
      if (exc(0,1,1) == 1) then
        ! Single alpha
        h = exc(1,1,1) ! hole in psi_det(1,1,j) 
        p = exc(1,2,1) ! particle in psi_det(1,1,j) 
      else
        ! Single beta
        h = exc(1,1,2) ! hole in psi_det(1,1,j) 
        p = exc(1,2,2) ! particle in psi_det(1,1,j) 
      endif
!      print*,'h,p',h,p
!      pause
      transition_matrix(p,h,istate,jstate)+= phase * psi_coef(i,istate) * psi_coef(j,jstate)
     endif
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER 
