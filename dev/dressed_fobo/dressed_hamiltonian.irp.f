BEGIN_PROVIDER [double precision, dressing_ref_fobo_hamiltonian, (n_det_ref_fobo,n_det_ref_fobo,N_states)]
 implicit none
 integer :: i,j,k,l
 integer :: ii,jj,istate
 double precision :: hij,sec_order,H_ref_fobo(N_det_ref_fobo),hik,hkl
 integer          :: idx(0:N_det_ref_fobo)
 double precision :: accu_negative,accu_positive,phase
 integer :: degree_exc_ionic,degree_exc_neutral,exc(0:2,2,2)
 dressing_ref_fobo_hamiltonian = 0.d0
 accu_negative = 0.d0
 accu_positive = 0.d0
 integer :: h1,p1,h2,p2,s1,s2
 print*, 'providing the dressing_ref_fobo_hamiltonian ...'
 ! IDEA : invert the loops over ref and non ref and parallelize the computation of HIJ
 do istate = 1, N_states
   do i = 1, N_det_non_ref_fobo
    call filter_connected_i_H_psi0(psi_ref_fobo,psi_non_ref_fobo(1,1,i),N_int,N_det_ref_fobo,idx)
    H_ref_fobo = 0.d0
    do ii=1,idx(0)
      k = idx(ii)
      !DEC$ FORCEINLINE
      call i_H_j(psi_ref_fobo(1,1,k),psi_non_ref_fobo(1,1,i),N_int,hij)
      H_ref_fobo(k)  = hij
    enddo
    do ii= 1, idx(0)
     k = idx(ii)
     hik = H_ref_fobo(k) * lambda_special(istate,i)
     do jj = 1, idx(0)
      l = idx(jj)
      ! IDEA : change the indices of the Hamiltonian 
      dressing_ref_fobo_hamiltonian(k,l,istate) += hik * H_ref_fobo(l)
     enddo
    enddo
   enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, hamiltonian_total_dressed, (n_det_ref_fobo,n_det_ref_fobo,N_states)]
 implicit none
 integer :: i,j,k
 do k = 1, N_states
  do i = 1, N_det_ref_fobo
   do j = 1, N_det_ref_fobo
    hamiltonian_total_dressed(j,i,k) = dressing_ref_fobo_hamiltonian(j,i,k) + ref_fobo_hamiltonian_matrix(j,i)
   enddo
  enddo
 enddo

END_PROVIDER 



 BEGIN_PROVIDER [ double precision, lambda_special, (N_states, N_det_non_ref_fobo) ]
&BEGIN_PROVIDER [ integer, lambda_special_pt2, (0:psi_det_size) ]
&BEGIN_PROVIDER [ integer, lambda_special_kept, (0:psi_det_size) ]
  implicit none
  BEGIN_DOC
  ! cm/<Psi_0|H|D_m> or perturbative 1/Delta_E(m)
  END_DOC
  integer :: i,k
  double precision               :: ihpsi_current(N_states)
  integer                        :: i_pert_count
  double precision               :: hii, lambda_pert
  integer                        :: N_lambda_special_pt2, N_lambda_special_pt3
  
  i_pert_count = 0
  lambda_special = 0.d0
  N_lambda_special_pt2 = 0
  N_lambda_special_pt3 = 0
  lambda_special_pt2(0) = 0
  lambda_special_kept(0) = 0

  print*, 'providing the amplitudes ...'
  do i=1,N_det_non_ref_fobo
    call i_h_psi(psi_non_ref_fobo(1,1,i), psi_ref_fobo, psi_ref_fobo_coef, N_int, N_det_ref_fobo,&
        size(psi_ref_fobo_coef,1), N_states,ihpsi_current)
    call i_H_j(psi_non_ref_fobo(1,1,i),psi_non_ref_fobo(1,1,i),N_int,hii)
    do k=1,N_states
      if (ihpsi_current(k) == 0.d0) then
        ihpsi_current(k) = 1.d-32
      endif
!      lambda_special(k,i) = psi_non_ref_fobo_coef(i,k)/ihpsi_current(k) 
      lambda_special(k,i) = min(-1.d-32,psi_non_ref_fobo_coef(i,k)/ihpsi_current(k) )
      lambda_pert = 1.d0 / (psi_ref_fobo_energy_diagonalized(k)-hii)
      if (lambda_pert / lambda_special(k,i)  < 0.299250d0) then
        ! Ignore lamdba
        i_pert_count += 1
        lambda_special(k,i) = 0.d0
        if (lambda_special_pt2(N_lambda_special_pt2) /= i) then
          N_lambda_special_pt2 += 1
          lambda_special_pt2(N_lambda_special_pt2) = i
        endif
      else
        ! Keep lamdba
        if (lambda_special_kept(N_lambda_special_pt3) /= i) then
          N_lambda_special_pt3 += 1
          lambda_special_kept(N_lambda_special_pt3) = i
        endif
      endif
    enddo
  enddo
  lambda_special_pt2(0) = N_lambda_special_pt2
  lambda_special_kept(0) = N_lambda_special_pt3
  print*,'N_det_non_ref_fobo = ',N_det_non_ref_fobo
  print*,'psi_coef_ref_fobo_ratio = ',psi_ref_fobo_coef(2,1)/psi_ref_fobo_coef(1,1)
  print*,'lambda max = ',maxval(dabs(lambda_special))
  print*,'Number of ignored determinants = ',i_pert_count  

END_PROVIDER

BEGIN_PROVIDER [double precision, ref_fobo_hamiltonian_matrix, (n_det_ref_fobo,n_det_ref_fobo)]
 BEGIN_DOC
 ! H matrix in the ref_foboerence space
 END_DOC
 implicit none
 integer :: i,j
 double precision :: hij
 do i = 1, N_det_ref_fobo
  do j = 1, N_det_ref_fobo
   call i_H_j(psi_ref_fobo(1,1,i),psi_ref_fobo(1,1,j),N_int,hij)
   ref_fobo_hamiltonian_matrix(i,j) = hij
  enddo
 enddo
END_PROVIDER

 BEGIN_PROVIDER [double precision, psi_ref_fobo_coef_diagonalized, (N_det_ref_fobo,N_states)]
&BEGIN_PROVIDER [double precision, psi_ref_fobo_energy_diagonalized, (N_states)]
 implicit none
 integer :: i,j
  double precision, allocatable  :: eigenvectors(:,:), eigenvalues(:)
  allocate (eigenvectors(size(ref_fobo_hamiltonian_matrix,1),N_det_ref_fobo))
  allocate (eigenvalues(N_det_ref_fobo))
  call lapack_diag(eigenvalues,eigenvectors,                       &
      ref_fobo_hamiltonian_matrix,size(ref_fobo_hamiltonian_matrix,1),N_det_ref_fobo)
  do i = 1, N_states
   psi_ref_fobo_energy_diagonalized(i) = eigenvalues(i)
   do j = 1, N_det_ref_fobo
    psi_ref_fobo_coef_diagonalized(j,i) = eigenvectors(j,i)
   enddo
  enddo
  deallocate (eigenvectors)
  deallocate (eigenvalues)


 END_PROVIDER

