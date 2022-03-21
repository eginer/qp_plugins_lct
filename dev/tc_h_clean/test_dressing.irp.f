program test_dressing
 implicit none
 read_wf  = .True.
 touch read_wf 
 call routine

end

subroutine routine
 implicit none
! provide reigvec_tc
 call test
 double precision, allocatable :: u0(:),h_dressed(:,:),eigenvalues(:),eigenvectors(:,:)
 integer :: idress
 allocate(u0(N_det),h_dressed(N_det,N_det),eigenvalues(N_det),eigenvectors(N_det,N_det))
 idress = 1
 u0(:) = psi_coef(:,1)
 call get_dressed_matrix(u0,psi_det,h_dressed,idress)
 call lapack_diag(eigenvalues,eigenvectors,h_dressed,size(h_dressed,1),N_det)
 print *,'eigenvalues(1)',eigenvalues(1)
end

subroutine test
 implicit none
 double precision, allocatable :: psicoef(:,:)
 double precision, allocatable :: H_jj(:), Dress_jj(:), Dressing_vec(:,:), energies(:)
 logical :: dagger,converged
 external H_u_0_nstates_openmp
 integer :: idx_dress,i,istate,j
 allocate(H_jj(N_det), Dress_jj(N_det), Dressing_vec(N_det,N_states))
 allocate(energies(N_states_diag))
 allocate(psicoef(N_det,N_states_diag))
 idx_dress = 1
 dagger = .False.
 !!! Create a Guess 
 psicoef = 0.d0
 do i = 1, N_det
  psicoef(i,1) = psi_coef(i,1)
 enddo
 do istate = N_states + 1, N_states_diag
  psicoef(istate,istate) = 1.d0
 enddo
 do j = 1, N_det
  H_jj(j) = H_matrix_diag_all_dets(j)
  Dress_jj(j) = 0.d0
 enddo
  call get_dressing_tc_for_dav(psicoef(1,1), psi_det, N_det, N_int, N_states, dagger, Dressing_vec)
  !!! Call the Davidson for dressing vector 
  call dav_double_dressed(psicoef,H_jj,Dress_jj,Dressing_vec,idx_dress,energies,N_det,N_states,N_states_diag,converged,H_u_0_nstates_openmp)
  print*, ' energies = ',energies(1)

end
