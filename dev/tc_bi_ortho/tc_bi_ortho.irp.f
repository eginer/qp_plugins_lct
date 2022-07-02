program tc_bi_ortho
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  read_wf = .True.
  touch read_wf
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
  PROVIDE mo_bi_ortho_tc_two_e mo_bi_ortho_tc_one_e
  call routine_diag
! call test
end

subroutine test
 implicit none
 integer :: i,j 
 double precision :: hmono,htwoe,hthree,htot
 use bitmasks

 print*,'test'
! call htilde_mu_mat_bi_ortho(psi_det(1,1,1), psi_det(1,1,2), N_int, hmono, htwoe, hthree, htot)
 call double_htilde_mu_mat_bi_ortho(N_int,psi_det(1,1,1), psi_det(1,1,2), hmono, htwoe, htot)
 print*,hmono, htwoe,  htot

end

subroutine routine_diag
 implicit none
! provide eigval_right_tc_bi_orth
  provide overlap_bi_ortho
  provide htilde_matrix_elmt_bi_ortho
 integer ::i,j
 print*,'eigval_right_tc_bi_orth = ',eigval_right_tc_bi_orth(1)
 print*,'e_tc_left_right         = ',e_tc_left_right
 print*,'e_tilde_bi_orth_00      = ',e_tilde_bi_orth_00
 print*,'e_pt2_tc_bi_orth        = ',e_pt2_tc_bi_orth
 print*,'e_pt2_tc_bi_orth_single = ',e_pt2_tc_bi_orth_single
 print*,'e_pt2_tc_bi_orth_double = ',e_pt2_tc_bi_orth_double
 print*,'***'
 print*,'e_corr_bi_orth          = ',e_corr_bi_orth
 print*,'e_corr_bi_orth_proj     = ',e_corr_bi_orth_proj
 print*,'e_corr_single_bi_orth   = ',e_corr_single_bi_orth
 print*,'e_corr_double_bi_orth   = ',e_corr_double_bi_orth
 print*,'Left/right eigenvectors'
 do i = 1,N_det
  write(*,'(I5,X,(100(F12.7,X)))')i,leigvec_tc_bi_orth(i,1),reigvec_tc_bi_orth(i,1)
 enddo
 do j=1,N_states
   do i=1,N_det
     psi_l_coef_bi_ortho(i,j) = leigvec_tc_bi_orth(i,j)
     psi_r_coef_bi_ortho(i,j) = reigvec_tc_bi_orth(i,j)
   enddo
 enddo
 SOFT_TOUCH  psi_l_coef_bi_ortho psi_r_coef_bi_ortho
 call routine_save_left_right_bi_ortho
end

!subroutine routine_three_e
! implicit none
! double precision :: new,ref, accu_diag,accu_single, accu_double
! double precision :: n_diag, n_single, n_double
! integer :: i,j,degree
! use bitmasks
! accu_diag = 0.d0
! accu_single = 0.d0
! accu_double = 0.d0
! n_diag = 0.d0
! n_single = 0.d0
! n_double = 0.D0
! do i = 1, N_det
!  do j = 1, N_det
!   call get_excitation_degree(psi_det(1,1,i), psi_det(1,1,j), degree, N_int)
!   if(degree==0)then
!    n_diag += 1.d0
!    call diag_htilde_three_body_ints_bi_ort(N_int, psi_det(1,1,i), new)
!    call diag_htilde_mu_mat_three_body(N_int, psi_det(1,1,i), ref)
!    if(dabs(ref - new).gt.1.d-10.and.dabs(ref).gt.1.d-10)then
!     print*,'Diagonal element '
!     print*,'ref,new,de'
!     print*,ref,new,ref-new
!    endif
!    accu_diag += dabs(ref-new)
!   else if(degree == 1)then
!    n_single += 1.d0
!    call single_htilde_three_body_ints_bi_ort(N_int,psi_det(1,1,i), psi_det(1,1,j), new)
!    call single_htilde_mu_mat_three_body(N_int,psi_det(1,1,i), psi_det(1,1,j), ref)
!    if(dabs(ref - new).gt.1.d-10.and.dabs(ref).gt.1.d-10)then
!     print*,'Single excitation element'
!     print*,'ref,new,de'
!     print*,ref,new,ref-new
!    endif
!    accu_single += dabs(ref-new)
!   else if(degree == 2)then
!    n_double += 1.d0
!    call double_htilde_mu_mat_three_body(N_int, psi_det(1,1,i), psi_det(1,1,j), new)
!    call double_htilde_mu_mat_three_body(N_int, psi_det(1,1,i), psi_det(1,1,j), ref)
!    if(dabs(ref - new).gt.1.d-10.and.dabs(ref).gt.1.d-10)then
!     print*,'Double excitation element'
!     print*,'ref,new,de'
!     print*,ref,new,ref-new
!    endif
!    accu_double += dabs(ref-new)
!   endif
!  enddo
! enddo
! 
! print*,'accu_diag = ',accu_diag / n_diag
! print*,'accu_singl= ',accu_single / n_single 
! print*,'accu_doubl= ',accu_double / n_double
!end
!
!subroutine test_whole_hij
! implicit none
! double precision :: new,ref, accu_diag,accu_single, accu_double
! double precision :: n_diag, n_single, n_double
! integer :: i,j,degree
! use bitmasks
! accu_diag = 0.d0
! accu_single = 0.d0
! accu_double = 0.d0
! n_diag = 0.d0
! n_single = 0.d0
! n_double = 0.D0
! do i = 1, N_det
!  do j = 1, N_det
!  do j = i, i
!   call htilde_mu_mat_bi_ortho_tot(psi_det(1,1,i), psi_det(1,1,j), N_int, new)
!   call htilde_mu_mat_tot(psi_det(1,1,i), psi_det(1,1,j), N_int, ref)
!   call get_excitation_degree(psi_det(1,1,i), psi_det(1,1,j), degree, N_int)
!   if(degree==0)then
!    n_diag += 1.d0
!    if(dabs(ref - new).gt.1.d-10.and.dabs(ref).gt.1.d-10)then
!     print*,'Diagonal element '
!     print*,'ref,new,de'
!     print*,ref,new,ref-new
!     print*,i,j
!     stop
!    endif
!    accu_diag += dabs(ref-new)
!   else if(degree == 1)then
!    n_single += 1.d0
!    if(dabs(ref - new).gt.1.d-10.and.dabs(ref).gt.1.d-10)then
!     print*,'Single excitation element'
!     print*,'ref,new,de'
!     print*,ref,new,ref-new
!    endif
!    accu_single += dabs(ref-new)
!   else if(degree == 2)then
!    n_double += 1.d0
!    if(dabs(ref - new).gt.1.d-10.and.dabs(ref).gt.1.d-10)then
!     print*,'Double excitation element'
!     print*,'ref,new,de'
!     print*,ref,new,ref-new
!    endif
!    accu_double += dabs(ref-new)
!   endif
!  enddo
! enddo
! 
! print*,'accu_diag = ',accu_diag / n_diag
! print*,'accu_singl= ',accu_single / n_single 
! print*,'accu_doubl= ',accu_double / n_double
! 
!
!
!end
