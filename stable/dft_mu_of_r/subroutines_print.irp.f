subroutine print_contribution_dft_mu_of_r
 implicit none

 print*,  '****************************************'
 print*, 'Functional used   = ',mu_of_r_functional
 print*,  '****************************************'
 print*, 'mu_of_r_potential = ',mu_of_r_potential
 print*,  ' MR DFT energy with pure correlation part for the DFT '
 if(mu_of_r_functional.EQ."basis_set_LDA")then
   print*, ''
   write(*, '(A28,X,F16.10)') 'Energy_c_md_mu_of_r_LDA   = ',Energy_c_md_mu_of_r_LDA
   print*, ''
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else if(mu_of_r_functional.EQ."basis_set_PBE")then
   write(*, '(A28,X,F16.10)') 'Energy ECMD PBE        = ',Energy_c_md_PBE_mu_of_r
   print*,''
   write(*, '(A28,X,F16.10)') 'Energy ECMD LDA        = ',Energy_c_md_mu_of_r_LDA
 else if(mu_of_r_functional.EQ."basis_set_on_top_PBE")then
   provide Energy_c_md_n_and_on_top_PBE_mu_of_r Energy_c_md_on_top_PBE_mu_of_r Energy_c_md_PBE_mu_of_r Energy_c_md_mu_of_r_LDA Energy_c_md_n_and_on_top_LYP_mu_of_r
   write(*, '(A28,X,F16.10)') 'Energy ECMD PBE ontop  = ',Energy_c_md_on_top_PBE_mu_of_r
   write(*, '(A28,X,F16.10)') 'Energy ECMD LYP ontop  = ',Energy_c_md_n_and_on_top_LYP_mu_of_r
   write(*, '(A28,X,F16.10)') 'ECMD PBE ontop/total n = ',Energy_c_md_n_and_on_top_PBE_mu_of_r
   write(*, '(A28,X,F16.10)') 'Energy ECMD PBE        = ',Energy_c_md_PBE_mu_of_r
   print*,''
   write(*, '(A28,X,F16.10)') 'Energy ECMD LDA        = ',Energy_c_md_mu_of_r_LDA
 endif
  if(.true.)then
   write(*, '(A28,X,F16.10)') 'mu_average for basis set  = ',mu_average
  endif

end

subroutine print_contribution_dft_mu_of_r_cc_cv_vv
 implicit none

 print*,  '****************************************'
 print*, 'Functional used   = ',mu_of_r_functional
 print*,  '****************************************'
 print*, 'mu_of_r_potential = ',mu_of_r_potential
 print*,  ' MR DFT energy with pure correlation part for the DFT '
!if(mu_of_r_functional.EQ."basis_set_LDA")then
   print*, ''
   write(*, '(A28,X,F16.10)') 'e_cmd_mu_of_r_lda_cc      = ',e_cmd_mu_of_r_lda_cc
   write(*, '(A28,X,F16.10)') 'e_cmd_mu_of_r_lda_vv      = ',e_cmd_mu_of_r_lda_vv
   write(*, '(A28,X,F16.10)') 'e_cmd_mu_of_r_lda_cv      = ',e_cmd_mu_of_r_lda_cv
   write(*, '(A28,X,F16.10)') 'e_cmd_mu_of_r_lda_cc+vv   = ',e_cmd_mu_of_r_lda_cc + e_cmd_mu_of_r_lda_vv
   write(*, '(A28,X,F16.10)') 'e_cmd_mu_of_r_lda_cc+vv+cv= ',e_cmd_mu_of_r_lda_cc + e_cmd_mu_of_r_lda_vv + e_cmd_mu_of_r_lda_cv 
   print*, ''
   write(*, '(A28,X,F16.10)') 'mu_average_hf_coal_cc     = ',mu_average_hf_coal_cc 
   write(*, '(A28,X,F16.10)') 'mu_average_hf_coal_cv     = ',mu_average_hf_coal_cv 
   write(*, '(A28,X,F16.10)') 'mu_average_hf_coal_vv     = ',mu_average_hf_coal_vv 
   print*, ''
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!else if(mu_of_r_functional.EQ."basis_set_PBE")then
   write(*, '(A28,X,F16.10)') 'Energy ECMD PBE        = ',Energy_c_md_PBE_mu_of_r
   print*,''
   write(*, '(A28,X,F16.10)') 'Energy ECMD LDA        = ',Energy_c_md_mu_of_r_LDA
!endif
  if(.true.)then
   write(*, '(A28,X,F16.10)') 'mu_average for basis set  = ',mu_average
  endif

end

