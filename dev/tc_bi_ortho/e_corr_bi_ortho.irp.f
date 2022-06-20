  use bitmasks ! you need to include the bitmasks_module.f90 features

 BEGIN_PROVIDER [ double precision, e_pt2_tc_bi_orth]
&BEGIN_PROVIDER [ double precision, e_pt2_tc_bi_orth_single]
&BEGIN_PROVIDER [ double precision, e_pt2_tc_bi_orth_double]
 implicit none 
 integer :: i,degree
 double precision :: hmono,htwoe,hthree,htilde_ij,coef_pt1,e_i0,delta_e
 e_pt2_tc_bi_orth = 0.d0
 e_pt2_tc_bi_orth_single = 0.d0
 e_pt2_tc_bi_orth_double = 0.d0
 do i = 1, N_det
  call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
  if(degree == 1 .or. degree == 2)then
   call htilde_mu_mat_bi_ortho(psi_det(1,1,i),HF_bitmask,N_int,hmono,htwoe,hthree,htilde_ij)
   call htilde_mu_mat_bi_ortho(psi_det(1,1,i),psi_det(1,1,i),N_int,hmono,htwoe,hthree,e_i0)
   delta_e = e_tilde_00 - e_i0
   coef_pt1 = htilde_ij / delta_e
   call htilde_mu_mat_bi_ortho(HF_bitmask,psi_det(1,1,i),N_int,hmono,htwoe,hthree,htilde_ij)
   e_pt2_tc_bi_orth += coef_pt1 * htilde_ij
   if(degree == 1)then
    e_pt2_tc_bi_orth_single += coef_pt1 * htilde_ij
   else 
    e_pt2_tc_bi_orth_double += coef_pt1 * htilde_ij
   endif
  endif
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, e_tilde_bi_orth_00]
 implicit none
 double precision :: hmono,htwoe,hthree,htilde_ij
 call htilde_mu_mat_bi_ortho(HF_bitmask,HF_bitmask,N_int,hmono,htwoe,hthree,e_tilde_bi_orth_00)
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, e_corr_bi_orth ]
&BEGIN_PROVIDER [ double precision, e_corr_bi_orth_proj ]
&BEGIN_PROVIDER [ double precision, e_corr_single_bi_orth ]
&BEGIN_PROVIDER [ double precision, e_corr_double_bi_orth ]
 implicit none 
 integer :: i,degree
 double precision :: hmono,htwoe,hthree,htilde_ij
 
 e_corr_bi_orth = 0.d0
 e_corr_single_bi_orth = 0.d0
 e_corr_double_bi_orth = 0.d0
 do i = 1, N_det
  call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
  call htilde_mu_mat_bi_ortho(HF_bitmask,psi_det(1,1,i),N_int,hmono,htwoe,hthree,htilde_ij)
  print*,reigvec_tc_bi_orth(i,1) , htilde_ij,reigvec_tc_bi_orth(1,1)
  if(degree == 1)then
   e_corr_single_bi_orth += reigvec_tc_bi_orth(i,1) * htilde_ij/reigvec_tc_bi_orth(1,1)
  else if(degree == 2)then
   e_corr_double_bi_orth += reigvec_tc_bi_orth(i,1) * htilde_ij/reigvec_tc_bi_orth(1,1)
  endif
 enddo
 e_corr_bi_orth_proj = e_corr_single_bi_orth + e_corr_double_bi_orth
 e_corr_bi_orth = eigval_right_tc_bi_orth(1) - e_tilde_bi_orth_00
 END_PROVIDER 

