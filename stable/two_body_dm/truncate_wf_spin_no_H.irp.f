program truncate
  read_wf = .True.
  SOFT_TOUCH read_wf
  call run
end

subroutine run
  use bitmasks
  implicit none
  integer :: i,j,k, kk, nab, m, l
  double precision :: norm, E, hij, num, ci, cj
  integer, allocatable :: iorder(:)
  double precision , allocatable :: norm_sort(:)
  double precision :: e_0(N_states)

  nab = n_det_alpha_unique+n_det_beta_unique
  allocate ( norm_sort(0:nab), iorder(0:nab) )

  integer(bit_kind), allocatable :: det_i(:,:), det_j(:,:)
  double precision, allocatable  :: u_t(:,:), v_t(:,:), s_t(:,:)
  double precision, allocatable  :: u_0(:,:), v_0(:,:)
  allocate(u_t(N_states,N_det),v_t(N_states,N_det),s_t(N_states,N_det))
  allocate(u_0(N_det,N_states),v_0(N_det,N_states))

  norm_sort(0) = 0.d0
  iorder(0) = 0
  do i=1,n_det_alpha_unique
   norm_sort(i) = det_alpha_norm(i)
   iorder(i) = i
  enddo
  
  do i=1,n_det_beta_unique
   norm_sort(i+n_det_alpha_unique) = det_beta_norm(i)
   iorder(i+n_det_alpha_unique) = -i
  enddo
  
  call dsort(norm_sort(1),iorder(1),nab)


  PROVIDE psi_bilinear_matrix_values psi_bilinear_matrix_rows psi_bilinear_matrix_columns
  PROVIDE nuclear_repulsion 
  print *,  ''
  do j=0,nab
    i = iorder(j)
    if (i<0) then
      !$OMP PARALLEL DO PRIVATE(k)
      do k=1,n_det
        if (psi_bilinear_matrix_columns(k) == -i) then
          do l=1,N_states
            psi_bilinear_matrix_values(k,l) = 0.d0
          enddo
        endif
      enddo
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(k)
      do k=1,n_det
        if (psi_bilinear_matrix_rows(k) ==  i) then
          do l=1,N_states
            psi_bilinear_matrix_values(k,l) = 0.d0
          enddo
        endif
      enddo
      !$OMP END PARALLEL DO
    endif
    if (ci_threshold <= norm_sort(j)) then
      exit
    endif
  enddo

  m = 0
  do k=1,n_det
    if (sum(psi_bilinear_matrix_values(k,1:N_states)) /= 0.d0) then
    m = m+1
    endif
  enddo

  do k=1,N_states
    E = E_0(k) + nuclear_repulsion
  enddo
  print *,  'Number of determinants:', m
  call wf_of_psi_bilinear_matrix(.True.)
  call save_wavefunction

  deallocate (iorder, norm_sort)
end

