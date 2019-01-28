 BEGIN_PROVIDER [double precision, psi_ref_fobo_coef_dressed, (n_det_ref_fobo,N_states) ]
&BEGIN_PROVIDER [double precision, energies_ref_fobo_dressed, (N_states) ]
 implicit none
 integer :: i,j,k,l,istate,igoodstate
 double precision, allocatable :: H_matrix_tmp(:,:)
 double precision, allocatable :: eigvalues(:),eigvectors(:,:),psi_coef_ref_fobo_tmp(:)
 double precision :: accu, accu1
 allocate(H_matrix_tmp(n_det_ref_fobo,n_det_ref_fobo))
 allocate(eigvalues(n_det_ref_fobo))
 allocate(eigvectors(n_det_ref_fobo,n_det_ref_fobo))
 allocate(psi_coef_ref_fobo_tmp(n_det_ref_fobo))
 do istate = 1, N_states
  accu1 = 0.d0
  do j = 1, n_det_ref_fobo
   accu1 += psi_ref_fobo_coef(j,istate)**2   ! norm of the "istate" eigenvector in the projected in the ref_foboerence space
   do k = 1, n_det_ref_fobo
    H_matrix_tmp(j,k) = hamiltonian_total_dressed(j,k,istate)
   enddo
  enddo
  accu1 = 1.d0/dsqrt(accu1)
  do j = 1, n_det_ref_fobo
   psi_coef_ref_fobo_tmp(j) = psi_ref_fobo_coef(j,istate) * accu1
  enddo
  call lapack_diagd(eigvalues,eigvectors,H_matrix_tmp,n_det_ref_fobo,n_det_ref_fobo)
  do j = 1, n_det_ref_fobo
   accu = 0.d0
   do k = 1, n_det_ref_fobo 
    accu += eigvectors(k,j) * psi_coef_ref_fobo_tmp(k)
   enddo
   print*, 'accu',accu
   if(dabs(accu).gt.0.9d0)then
    igoodstate = j
    exit
   endif
  enddo
  energies_ref_fobo_dressed(istate) = eigvalues(igoodstate)
  do j = 1,n_det_ref_fobo
   psi_ref_fobo_coef_dressed(j,istate) = eigvectors(j,igoodstate)
  enddo
 enddo
 
 print*,''
 print*,'igoodstate = ',igoodstate
 print*,'extra informations on the spectrum of the dressed matrix '
 print*,''

 double precision, allocatable :: s2_vec(:)
 allocate(s2_vec(n_det_ref_fobo))
 s2_vec = 0.d0

 double precision :: s2
 do i = 1, min(3,n_det_ref_fobo)
  print*,'eigvalues = ',i,eigvalues(i) + nuclear_repulsion
  do k = 1, n_det_ref_fobo 
   do l = 1, n_det_ref_fobo
    call get_s2(psi_ref_fobo(1,1,l),psi_ref_fobo(1,1,k),N_int,s2)
    s2_vec(i) += eigvectors(k,i) * eigvectors(l,i) * s2 
   enddo
  enddo
  print*,'< S2 >    = ',s2_vec(i)
 enddo
 print*,''
 do i = 2, min(3,n_det_ref_fobo)
  print*,'DE = ',eigvalues(i) -eigvalues(1)
 enddo
 print*,''

END_PROVIDER 
