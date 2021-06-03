 BEGIN_PROVIDER [double precision, htilde_matrix_elmt, (N_det,N_det)]
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
   call htilde_mat(psi_det(1,1,j),psi_det(1,1,i),hmono,heff,hderiv,hthree,htot)
   htilde_matrix_elmt(j,i) = htot
   htilde_matrix_elmt_eff(j,i) = heff
   htilde_matrix_elmt_deriv(j,i) = hderiv
   htilde_matrix_elmt_hcore(j,i) = hmono
   htilde_matrix_elmt_hthree(j,i) = hthree
  enddo
 enddo
 do i = 1, N_det
  write(*,'(1000(F10.5,X))')htilde_matrix_elmt(i,:)
 enddo
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
 else
  call non_hrmt_real_diag(N_det,htilde_matrix_elmt,reigvec_trans,leigvec_trans,n_good_trans_eigval,eigval_trans)
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

subroutine test_left_right_eigenvalues(ith)
 implicit none
 integer, intent(in) :: ith
 BEGIN_DOC
! rights the "ith" eigenvalue in many different ways which must all be the same
 END_DOC
 integer :: i,j,k
 ! test for the right eigenvector
 double precision :: accu1
 print*,'Printing out the Right, Left and usual eigenvectors'
 print*,'All have been normalized'
 print*,'Norm = ',reigvec_trans_norm(ith),leigvec_trans_norm(ith)
 print*,''
 print*,'Right      Left      Usual'
 do i = 1, N_det
  write(*,'(3(F10.7,X))')reigvec_trans(i,ith)/dsqrt(reigvec_trans_norm(ith)) & 
                        ,leigvec_trans(i,ith)/dsqrt(leigvec_trans_norm(ith)),psi_coef(i,ith)
 enddo
 print*,''
 print*,'Eigenvalue of the usual Hamiltonian '
 print*,'',ci_electronic_energy(ith)
 print*,''
 print*,'Checking the eigenvalue as an expectation value'
 print*,'With the RIGHT eigenvector'

 ! Expectation value 
 accu1 = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   ! \sum_IJ c^R_J <J|H^tilde|I> c^R_I = E * <Psi_ith | Psi_ith>
   accu1 += reigvec_trans(j,ith) * htilde_matrix_elmt(j,i) * reigvec_trans(i,ith)
  enddo
 enddo
 print*,'accu1/norm = ',accu1/reigvec_trans_norm(ith) 
 print*,'eigval_trans ',eigval_trans(ith)

 print*,'With the LEFT  eigenvector'
 accu1 = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   ! \sum_IJ c^L_J <J|H^tilde|I> c^L_I = E * <Psi_ith | Psi_ith>
   accu1 += leigvec_trans(j,ith) * htilde_matrix_elmt(j,i) * leigvec_trans(i,ith)
  enddo
 enddo
 print*,'accu1/norm = ',accu1/leigvec_trans_norm(ith) 
 print*,'eigval_trans ',eigval_trans(ith)

 print*,''
 print*,'Testing the projection scheme'
 print*,'with the RIGHT eigenvector'
 ! Projection
 accu1 = 0.d0
 do i = 1, N_det 
  ! \sum_I <1|H^tilde|I> c^R_I = E * c_1^R 
  accu1 += htilde_matrix_elmt(1,i) * reigvec_trans(i,ith)
 enddo
 print*,'accu1/norm = ',accu1/reigvec_trans(1,1)
 print*,'eigval_trans ',eigval_trans(ith)
 print*,''
 print*,'with the LEFT  eigenvector'
 ! Projection
 accu1 = 0.d0
 do i = 1, N_det 
  ! \sum_I  <I|H^tilde|1> c^L_I = E * * c_1^L
  accu1 +=  htilde_matrix_elmt(i,1) * leigvec_trans(i,ith)
 enddo
 print*,'accu1/norm = ',accu1/leigvec_trans(1,1)
 print*,'eigval_trans ',eigval_trans(ith)
 print*,''

end

subroutine test_overlap_matrix
 implicit none
 integer :: i,j,k
 double precision :: accurr,accull,acculr
 double precision, allocatable :: overlap_rr(:,:)
 double precision, allocatable :: overlap_ll(:,:)
 double precision, allocatable :: overlap_lr(:,:)
 allocate(overlap_rr(n_good_trans_eigval,n_good_trans_eigval), &
          overlap_ll(n_good_trans_eigval,n_good_trans_eigval), & 
          overlap_lr(n_good_trans_eigval,n_good_trans_eigval))
 do i = 1, n_good_trans_eigval
  do j = 1, n_good_trans_eigval
   accurr = 0.d0
   accull = 0.d0
   acculr = 0.d0
   do k = 1, N_det
    accurr += reigvec_trans(k,i) * reigvec_trans(k,j)
    accull += leigvec_trans(k,i) * leigvec_trans(k,j)
    acculr += leigvec_trans(k,i) * reigvec_trans(k,j)
   enddo
   overlap_rr(j,i) = accurr
   overlap_ll(j,i) = accull
   overlap_lr(j,i) = acculr
  enddo
 enddo

 print*,''
 print*,'Printing the overlap between RIGHT eigenvectors'
 do i = 1, n_good_trans_eigval
  write(*,'(1000(F16.12,X))')overlap_rr(i,:)
 enddo

 print*,''
 print*,'Printing the overlap between LEFT  eigenvectors'
 do i = 1, n_good_trans_eigval
  write(*,'(1000(F16.12,X))')overlap_ll(i,:)
 enddo

 print*,''
 print*,'Printing the overlap between LEFT and RIGHT  eigenvectors'
 do i = 1, n_good_trans_eigval
  write(*,'(1000(F16.12,X))')overlap_lr(i,:)
 enddo
end

BEGIN_PROVIDER [double precision, overlap_psi_det_r_eigevec, (N_states)]
 implicit none
 double precision :: accu
 BEGIN_DOC
! overlap between psi_det and the right eigenvectors of Htilde
 END_DOC
 integer :: i,j
 do i = 1, N_states
  accu = 0.d0
  do j = 1, N_det
   accu += psi_coef(j,i) * reigvec_trans(j,i)
  enddo
  overlap_psi_det_r_eigevec(i) = accu
 enddo
END_PROVIDER 
