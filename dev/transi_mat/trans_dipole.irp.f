
 BEGIN_PROVIDER [ double precision, trans_dipole_x, (N_states,N_states) ]
&BEGIN_PROVIDER [ double precision, trans_dipole_y, (N_states,N_states) ]
&BEGIN_PROVIDER [ double precision, trans_dipole_z, (N_states,N_states) ]
  implicit none
  BEGIN_DOC
  ! $\alpha$ and $\beta$ one-body density matrix for each state
  END_DOC

  integer                        :: j,k,l,m,k_a,k_b,n
  integer                        :: occ(N_int*bit_kind_size,2)
  double precision               :: ck, cl, ckl
  double precision               :: phase
  integer                        :: h1,h2,p1,p2,s1,s2, degree
  integer(bit_kind)              :: tmp_det(N_int,2), tmp_det2(N_int)
  integer                        :: exc(0:2,2),n_occ(2)
  double precision, allocatable  :: tmp_x(:,:), tmp_y(:,:), tmp_z(:,:)
  integer                        :: krow, kcol, lrow, lcol

  PROVIDE psi_det mo_dipole_x mo_dipole_y mo_dipole_z 

  print*,'providing the trans_dipole  '
  trans_dipole_x = 0.d0
  trans_dipole_y = 0.d0
  trans_dipole_z = 0.d0
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(j,k,n,k_a,k_b,l,m,occ,ck, cl, ckl,phase,h1,h2,p1,p2,s1,s2, degree,exc,&
      !$OMP  tmp_x, tmp_y, tmp_z, n_occ, krow, kcol, lrow, lcol, tmp_det, tmp_det2)&
      !$OMP SHARED(psi_det,psi_coef,N_int,N_states,elec_alpha_num, mo_dipole_x, mo_dipole_y, mo_dipole_z, &
      !$OMP  elec_beta_num,trans_dipole_x,trans_dipole_y,trans_dipole_z,N_det,&
      !$OMP  psi_bilinear_matrix_rows,psi_bilinear_matrix_columns,&
      !$OMP  psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns,&
      !$OMP  psi_bilinear_matrix_order_reverse, psi_det_alpha_unique, psi_det_beta_unique,&
      !$OMP  psi_bilinear_matrix_values, psi_bilinear_matrix_transp_values,&
      !$OMP  N_det_alpha_unique,N_det_beta_unique,irp_here)
  allocate(tmp_x(N_states,N_states), tmp_y(N_states,N_states), tmp_z(N_states,N_states))
  tmp_x = 0.d0
  tmp_y = 0.d0
  tmp_z = 0.d0
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
         tmp_x(n,m) += ck * mo_dipole_x(j,j)
         tmp_y(n,m) += ck * mo_dipole_y(j,j)
         tmp_z(n,m) += ck * mo_dipole_z(j,j)
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
          tmp_x(n,m) +=  ckl * mo_dipole_x(h1,p1)
          tmp_y(n,m) +=  ckl * mo_dipole_y(h1,p1)
          tmp_z(n,m) +=  ckl * mo_dipole_z(h1,p1)
          ckl = psi_bilinear_matrix_values(k_a,n)*psi_bilinear_matrix_values(l,m) * phase
          tmp_x(n,m) +=  ckl * mo_dipole_x(h1,p1)
          tmp_y(n,m) +=  ckl * mo_dipole_y(h1,p1)
          tmp_z(n,m) +=  ckl * mo_dipole_z(h1,p1)
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
   trans_dipole_x += tmp_x 
   trans_dipole_y += tmp_y
   trans_dipole_z += tmp_z
  !$OMP END CRITICAL
  deallocate(tmp_x,tmp_y,tmp_z)
  allocate(tmp_x(N_states,N_states), tmp_y(N_states,N_states), tmp_z(N_states,N_states))
  tmp_x = 0.d0
  tmp_y = 0.d0
  tmp_z = 0.d0

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
         tmp_x(n,m) += ck * mo_dipole_x(j,j)
         tmp_y(n,m) += ck * mo_dipole_y(j,j)
         tmp_z(n,m) += ck * mo_dipole_z(j,j)
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
          tmp_x(n,m) +=  ckl * mo_dipole_x(h1,p1)
          tmp_y(n,m) +=  ckl * mo_dipole_y(h1,p1)
          tmp_z(n,m) +=  ckl * mo_dipole_z(h1,p1)
          ckl = psi_bilinear_matrix_transp_values(k_b,n)*psi_bilinear_matrix_transp_values(l,m) * phase
          tmp_x(n,m) +=  ckl * mo_dipole_x(h1,p1)
          tmp_y(n,m) +=  ckl * mo_dipole_y(h1,p1)
          tmp_z(n,m) +=  ckl * mo_dipole_z(h1,p1)
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
   trans_dipole_x += tmp_x 
   trans_dipole_y += tmp_y
   trans_dipole_z += tmp_z
  !$OMP END CRITICAL

  deallocate(tmp_x,tmp_y,tmp_z)
  !$OMP END PARALLEL
  print*,'provided the trans_dipole '

END_PROVIDER

 BEGIN_PROVIDER [double precision, trans_dipole_x_bourrin, (N_states,N_states)]
&BEGIN_PROVIDER [double precision, trans_dipole_y_bourrin, (N_states,N_states)]
&BEGIN_PROVIDER [double precision, trans_dipole_z_bourrin, (N_states,N_states)]
 implicit none
 integer :: i,j,m,n
 double precision :: xij,yij,zij

 trans_dipole_x_bourrin = 0.d0
 trans_dipole_y_bourrin = 0.d0
 trans_dipole_z_bourrin = 0.d0

  print*,'providing the trans_dipole_bourrin '
 do i = 1, N_det
  do j = 1, N_det
   call i_H_j_eff_pot(psi_det(1,1,i),psi_det(1,1,j),mo_dipole_x,mo_dipole_x,mo_num,N_int,xij)
   call i_H_j_eff_pot(psi_det(1,1,i),psi_det(1,1,j),mo_dipole_y,mo_dipole_y,mo_num,N_int,yij)
   call i_H_j_eff_pot(psi_det(1,1,i),psi_det(1,1,j),mo_dipole_z,mo_dipole_z,mo_num,N_int,zij)
   do m = 1, N_states
    do n = 1, N_states
     trans_dipole_x_bourrin(n,m) += psi_coef(i,n) * psi_coef(j,m) * xij
     trans_dipole_y_bourrin(n,m) += psi_coef(i,n) * psi_coef(j,m) * yij
     trans_dipole_z_bourrin(n,m) += psi_coef(i,n) * psi_coef(j,m) * zij
    enddo
   enddo
  enddo
 enddo
  print*,'provided the trans_dipole_bourrin '


END_PROVIDER 
