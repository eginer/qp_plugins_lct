 BEGIN_PROVIDER [ integer, index_HF_psi_det]                                                                                                            
 implicit none
 integer :: i,degree
 do i = 1, N_det
   call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
   if(degree == 0)then
    index_HF_psi_det = i
    exit
   endif
 enddo
 END_PROVIDER



 BEGIN_PROVIDER [double precision, eigval_right_tc_bi_orth, (N_states)]
&BEGIN_PROVIDER [double precision, eigval_left_tc_bi_orth, (N_states)]
&BEGIN_PROVIDER [double precision, reigvec_tc_bi_orth, (N_det,N_states)]
&BEGIN_PROVIDER [double precision, leigvec_tc_bi_orth, (N_det,N_states)]

  BEGIN_DOC
  ! eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian 
  END_DOC

  implicit none
  integer                       :: i, idx_dress, j
  logical                       :: converged, dagger
  integer                       :: n_real_tc_bi_orth_eigval_right,igood_r,igood_l
  double precision, allocatable :: reigvec_tc_bi_orth_tmp(:,:),leigvec_tc_bi_orth_tmp(:,:),eigval_right_tmp(:)

  PROVIDE N_det N_int


    allocate(reigvec_tc_bi_orth_tmp(N_det,N_det),leigvec_tc_bi_orth_tmp(N_det,N_det),eigval_right_tmp(N_det))
    call non_hrmt_real_diag(N_det,htilde_matrix_elmt_bi_ortho,& 
         leigvec_tc_bi_orth_tmp,reigvec_tc_bi_orth_tmp,& 
         n_real_tc_bi_orth_eigval_right,eigval_right_tmp)
    double precision, allocatable :: coef_hf_r(:),coef_hf_l(:)
    integer, allocatable :: iorder(:)
    allocate(coef_hf_r(N_det),coef_hf_l(N_det),iorder(N_det))
    do i = 1,N_det
     iorder(i) = i
     coef_hf_r(i) = -dabs(reigvec_tc_bi_orth_tmp(index_HF_psi_det,i))
    enddo
    call dsort(coef_hf_r,iorder,N_det)
    igood_r = iorder(1)
    print*,'igood_r, coef_hf_r = ',igood_r,coef_hf_r(1)
    do i = 1,N_det
     iorder(i) = i
     coef_hf_l(i) = -dabs(leigvec_tc_bi_orth_tmp(index_HF_psi_det,i))
    enddo
    call dsort(coef_hf_l,iorder,N_det)
    igood_l = iorder(1)
    print*,'igood_l, coef_hf_l = ',igood_l,coef_hf_l(1)

    if(igood_r.ne.igood_l.and.igood_r.ne.1)then
     print *,''
     print *,'Warning, the left and right eigenvectors are "not the same" '
     print *,'Warning, the ground state is not dominated by HF...'
     print *,'State with largest RIGHT coefficient of HF ',igood_r
     print *,'coef of HF in RIGHT eigenvector = ',reigvec_tc_bi_orth_tmp(index_HF_psi_det,igood_r)
     print *,'State with largest LEFT  coefficient of HF ',igood_l
     print *,'coef of HF in LEFT  eigenvector = ',leigvec_tc_bi_orth_tmp(index_HF_psi_det,igood_l)
    endif
    if(state_following_tc)then
     print *,'Following the states with the largest coef on HF'
     print *,'igood_r,igood_l',igood_r,igood_l
     i= igood_r
     eigval_right_tc_bi_orth(1) = eigval_right_tmp(i)
     do j = 1, N_det
       reigvec_tc_bi_orth(j,1) = reigvec_tc_bi_orth_tmp(j,i)
     enddo
     i= igood_l
     eigval_left_tc_bi_orth(1)  = eigval_right_tmp(i)
     do j = 1, N_det
       leigvec_tc_bi_orth(j,1) = leigvec_tc_bi_orth_tmp(j,i)
     enddo
    else 
     do i = 1, N_states
       eigval_right_tc_bi_orth(i) = eigval_right_tmp(i)
       eigval_left_tc_bi_orth(i)  = eigval_right_tmp(i)
       do j = 1, N_det
         reigvec_tc_bi_orth(j,i) = reigvec_tc_bi_orth_tmp(j,i)
         leigvec_tc_bi_orth(j,i) = leigvec_tc_bi_orth_tmp(j,i)
       enddo
     enddo
    endif

  !!!! Normalization of the scalar product of the left/right eigenvectors
  double precision  :: accu, tmp 
  do i = 1, N_states
   !!!! Normalization of right eigenvectors |Phi>
   accu = 0.d0
   do j = 1, N_det
    accu += reigvec_tc_bi_orth(j,i) * reigvec_tc_bi_orth(j,i)
   enddo
   accu = 1.d0/dsqrt(accu)
   do j = 1, N_det
    reigvec_tc_bi_orth(j,i) *= accu 
   enddo
   tmp = reigvec_tc_bi_orth(1,i) / dabs(reigvec_tc_bi_orth(1,i))
   do j = 1, N_det
    reigvec_tc_bi_orth(j,i) *= tmp
   enddo
   !!!! Adaptation of the norm of the left eigenvector such that <chi|Phi> = 1
   accu = 0.d0
   do j = 1, N_det
    accu += leigvec_tc_bi_orth(j,i) * reigvec_tc_bi_orth(j,i)
   enddo
   if(accu.gt.0.d0)then
    accu = 1.d0/dsqrt(accu)
   else
    accu = 1.d0/dsqrt(-accu)
   endif
   tmp = (leigvec_tc_bi_orth(1,i) * reigvec_tc_bi_orth(1,i) )/dabs(leigvec_tc_bi_orth(1,i) * reigvec_tc_bi_orth(1,i))
   do j = 1, N_det
    leigvec_tc_bi_orth(j,i) *= accu * tmp
    reigvec_tc_bi_orth(j,i) *= accu 
   enddo
   print*,'leigvec_tc_bi_orth(1,i),reigvec_tc_bi_orth(1,i) = ',leigvec_tc_bi_orth(1,i),reigvec_tc_bi_orth(1,i)
   accu = 0.d0
   do j = 1, N_det
    accu += leigvec_tc_bi_orth(j,i) * reigvec_tc_bi_orth(j,i)
   enddo
   print*,'norm l/r = ',accu
  enddo

END_PROVIDER 
