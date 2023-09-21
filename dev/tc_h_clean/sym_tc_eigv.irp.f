 BEGIN_PROVIDER [double precision, eigval_sym_tc, (N_det)]
&BEGIN_PROVIDER [double precision, eigvec_sym_tc, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, exp_h_tilde_eigvec_sym_tc]
  implicit none
  integer :: i,idx_dress,j
  logical :: converged, dagger
 BEGIN_DOC
 ! eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian 
 END_DOC
! if(full_sym_tc_h_solver.or.N_det.le.N_det_max_full)then
   double precision, allocatable :: H_sym_tc(:,:)
   allocate(H_sym_tc(N_det,N_det))
   H_sym_tc = 0.5d0 * (htilde_matrix_elmt_tranp + htilde_matrix_elmt)
   call lapack_diag(eigval_sym_tc,eigvec_sym_tc,H_sym_tc,size(H_sym_tc,1),N_det)
   exp_h_tilde_eigvec_sym_tc = 0.d0
   do i = 1, N_det
    do j = 1, N_det
     exp_h_tilde_eigvec_sym_tc += eigvec_sym_tc(i,1) * eigvec_sym_tc(j,1) * htilde_matrix_elmt(i,j)
    enddo
   enddo
! else
!   allocate(reigvec_sym_tc_tmp(N_det,N_states_diag))
!   idx_dress = 1
!   dagger = .False. ! to compute the RIGHT eigenvector 
!   reigvec_sym_tc_tmp = 0.d0
!   do i = 1, N_states
!    do j = 1, N_det
!     reigvec_sym_tc_tmp(j,i) = psi_coef(j,i)
!    enddo
!   enddo
!   call iterative_davidson_sym_tc(psi_det,reigvec_sym_tc_tmp,N_det,N_states,N_states_diag,idx_dress,dagger, & 
!                              eigval_right_sym_tc,converged)
!   do i = 1, N_states
!    do j = 1, N_det
!     reigvec_sym_tc(j,i) = reigvec_sym_tc_tmp(j,i) 
!    enddo
!   enddo
!   if(.not.converged)then
!    print*,'Stopping ... Davidson did not converge ...'
!    stop
!   endif
!   if(comp_left_eigv)then
!     idx_dress = 1
!     dagger = .True. ! to compute the LEFT  eigenvector 
!     reigvec_sym_tc_tmp = 0.d0
!     do i = 1, N_states
!      do j = 1, N_det
!       reigvec_sym_tc_tmp(j,i) = psi_coef(j,i)
!      enddo
!     enddo
!     call iterative_davidson_sym_tc(psi_det,reigvec_sym_tc_tmp,N_det,N_states,N_states_diag,idx_dress,dagger, & 
!                                eigval_left_sym_tc,converged)
!     do i = 1, N_states
!      do j = 1, N_det
!       leigvec_sym_tc(j,i) = reigvec_sym_tc_tmp(j,i) 
!      enddo
!     enddo
!   endif
! endif
END_PROVIDER 

