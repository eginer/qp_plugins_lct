 BEGIN_PROVIDER [double precision, htilde_matrix_elmt, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, htilde_matrix_elmt_tranp, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, htilde_matrix_elmt_eff, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, htilde_matrix_elmt_deriv, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, htilde_matrix_elmt_hcore, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, htilde_matrix_elmt_hthree, (N_det,N_det)]
 implicit none
 BEGIN_DOC
! htilde_matrix_elmt(j,i) = <J| H^tilde |I> 
!
! WARNING !!!!!!!!! IT IS NOT HERMITIAN !!!!!!!!!
 END_DOC
 integer :: i,j
 double precision :: hmono,heff,hderiv,hthree,htot
 do i = 1, N_det
  do j = 1, N_det
  ! < J | Htilde | I >
   call htilde_mu_mat(psi_det(1,1,j),psi_det(1,1,i),hmono,heff,hderiv,hthree,htot)
   htilde_matrix_elmt(j,i) = htot
   htilde_matrix_elmt_eff(j,i) = heff
   htilde_matrix_elmt_deriv(j,i) = hderiv
   htilde_matrix_elmt_hcore(j,i) = hmono
   htilde_matrix_elmt_hthree(j,i) = hthree
  enddo
 enddo
 do i = 1, N_det
  do j = 1, N_det
   htilde_matrix_elmt_tranp(j,i) = htilde_matrix_elmt(i,j)
  enddo
 enddo
! do i = 1, N_det
!  write(*,'(1000(F10.5,X))')htilde_matrix_elmt(i,:)
! enddo
! htilde_matrix_elmt = H_matrix_all_dets
END_PROVIDER 


 BEGIN_PROVIDER [double precision, eigval_trans, (N_det)]
&BEGIN_PROVIDER [integer, n_good_trans_eigval]
&BEGIN_PROVIDER [double precision, reigvec_trans, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, leigvec_trans, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, reigvec_trans_norm, (N_det)]
&BEGIN_PROVIDER [double precision, leigvec_trans_norm, (N_det)]
 implicit none
 integer :: i,j
 BEGIN_DOC
! n_good_trans_eigval = number of PURE REAL eigenvalues for H^tilde (should be all if everything is fine)
!
! eigval_trans(i) = ith eigenvalue of H^tilde (sorted by increasing order)
!
! reigvec_trans(j,i) = <J|Psi_i> where |Psi_i> is the RIGHT EIGENVECTOR
!
! Fulfilling H^tilde |Psi_i> = E_i |Psi_i> 
!
! reigvec_trans_norm(i) = \sum_J |<J|Psi_i>|^2 (not normalized in principle)
!
! leigvec_trans(j,i) = <J|Psi_i> where |Psi_i> is the LEFT  EIGENVECTOR
!
! Fulfilling (H^tilde)^transpose |Psi_i> = E_i |Psi_i> 
!
! leigvec_trans_norm(i) = \sum_J |<J|Psi_i>|^2 (not normalized in principle)
!
 END_DOC

 if(read_rl_eigv)then
  reigvec_trans = 0.d0 
  leigvec_trans = 0.d0 
  double precision, allocatable :: reigvec_trans_tmp(:,:),leigvec_trans_tmp(:,:)
  allocate(reigvec_trans_tmp(N_det,N_states),leigvec_trans_tmp(N_det,N_states))
  print*,'Reading the right/left eigenvectors in the EZFIO....'
   call ezfio_get_transcorr_h_reigvec_trans(reigvec_trans_tmp)
   call ezfio_get_transcorr_h_leigvec_trans(leigvec_trans_tmp)
   reigvec_trans_norm = 0.d0
   leigvec_trans_norm = 0.d0
   do i = 1, N_states
    do j = 1, N_det
     reigvec_trans(j,i) = reigvec_trans_tmp(j,i)
     leigvec_trans(j,i) = leigvec_trans_tmp(j,i)
     reigvec_trans_norm(i) += reigvec_trans_tmp(j,i) * reigvec_trans_tmp(j,i)
     leigvec_trans_norm(i) += leigvec_trans_tmp(j,i) * leigvec_trans_tmp(j,i)
    enddo
   enddo
 else if(full_tc_h_solver)then
  call non_hrmt_real_diag(N_det,htilde_matrix_elmt,reigvec_trans,leigvec_trans,n_good_trans_eigval,eigval_trans)
  do i = 1, n_good_trans_eigval
   reigvec_trans_norm(i) = 0.d0
   leigvec_trans_norm(i) = 0.d0
   do j = 1, N_det
    reigvec_trans_norm(i) += reigvec_trans(j,i) * reigvec_trans(j,i)
    leigvec_trans_norm(i) += leigvec_trans(j,i) * leigvec_trans(j,i)
   enddo
  enddo
 else
  eigval_trans = 0.d0
  reigvec_trans = 0.d0
  leigvec_trans = 0.d0
  n_good_trans_eigval = 1
  double precision :: e0
  double precision, allocatable :: u(:), v(:)
  print*,'Non hermitian Davidson to be plugged in here ...'
  stop
  allocate(u(N_det), v(N_det))
  u = 0.d0
  u(1) = 1.d0
!  call project_ground(u,v,htilde_matrix_elmt,e0,1,N_det)
  do j = 1, N_det
   reigvec_trans(j,1) = v(j)  
  enddo
  print*,'e0 from right eigenvector = ',e0
  eigval_trans(1) = e0

!  call project_ground(u,v,htilde_matrix_elmt_tranp,e0,1,N_det)
  do j = 1, N_det
   leigvec_trans(j,1) = v(j)  
  enddo
  print*,'e0 from left  eigenvector = ',e0
  do i = 1, n_good_trans_eigval
   reigvec_trans_norm(i) = 0.d0
   leigvec_trans_norm(i) = 0.d0
   do j = 1, N_det
    reigvec_trans_norm(i) += reigvec_trans(j,i) * reigvec_trans(j,i)
    leigvec_trans_norm(i) += leigvec_trans(j,i) * leigvec_trans(j,i)
   enddo
  enddo
 endif

END_PROVIDER 


 BEGIN_PROVIDER [ double precision, diag_htilde, (N_det)]
 implicit none
 integer :: i
 double precision :: hmono,heff,hderiv,hthree,htot
 do i = 1, N_det
  call htilde_mu_mat(psi_det(1,1,i),psi_det(1,1,i),hmono,heff,hderiv,hthree,htot)
  diag_htilde(i) = htot
 enddo
 END_PROVIDER 

BEGIN_PROVIDER [ integer, n_tc_ovlp_print]
 implicit none
 n_tc_ovlp_print = min(n_good_trans_eigval,10)
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, left_right_overlap, (n_tc_ovlp_print, n_tc_ovlp_print)]
&BEGIN_PROVIDER [ double precision, left_left_overlap, (n_tc_ovlp_print, n_tc_ovlp_print)]
&BEGIN_PROVIDER [ double precision, right_right_overlap, (n_tc_ovlp_print, n_tc_ovlp_print)]
 implicit none
 double precision :: accu_lr, accu_ll, accu_rr
 integer :: i,j,k,l
 print*,''
 print*,''
 print*,'Computing overlap between states'
 print*,''
 print*,''
 do i = 1, n_tc_ovlp_print
  do j = 1, n_tc_ovlp_print
   accu_rr = 0.d0
   accu_lr = 0.d0
   accu_ll = 0.d0
   do l = 1, N_det
    accu_rr += reigvec_trans(l,j) * reigvec_trans(l,i) 
    accu_lr += leigvec_trans(l,j) * reigvec_trans(l,i) 
    accu_ll += leigvec_trans(l,j) * leigvec_trans(l,i) 
   enddo
   right_right_overlap(j,i) = accu_rr
   left_right_overlap(j,i) = accu_lr
   left_left_overlap(j,i) = accu_ll
  enddo
 enddo
END_PROVIDER 


subroutine write_left_right
 implicit none
 double precision, allocatable :: reigvec_trans_tmp(:,:),leigvec_trans_tmp(:,:)
 allocate(leigvec_trans_tmp(N_det,N_states),reigvec_trans_tmp(N_det,N_states))
 integer :: i,j
 do i = 1, N_states
  do j = 1, N_det
   leigvec_trans_tmp(j,i) = leigvec_trans(j,i)
   reigvec_trans_tmp(j,i) = reigvec_trans(j,i)
  enddo
 enddo
 call ezfio_set_transcorr_h_reigvec_trans(reigvec_trans_tmp)
 call ezfio_set_transcorr_h_leigvec_trans(leigvec_trans_tmp)
end

