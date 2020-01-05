subroutine print_contribution_dft_mu_of_r
 implicit none

 print*,  '****************************************'
 print*, 'Functional used   = ',mu_of_r_functional
 print*,  '****************************************'
 print*, 'mu_of_r_potential = ',mu_of_r_potential
 print*,  ' MR DFT energy with pure correlation part for the DFT '
 if(mu_of_r_functional.EQ."basis_set_LDA")then
   print*, ''
   write(*, '(A28,X,F16.10)') 'Energy_c_md_LDA_mu_of_r   = ',Energy_c_md_LDA_mu_of_r
   print*, ''
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else if(mu_of_r_functional.EQ."basis_set_PBE")then
   print*,''
   print*,'Corrections using single determinant mu'
   print*,''
   write(*, '(A28,X,F16.10)') 'Energy ECMD LDA        = ',Energy_c_md_LDA_mu_of_r
   print*,''
   write(*, '(A28,X,F16.10)') 'Energy ECMD PBE        = ',Energy_c_md_PBE_mu_of_r
!  print*,''
!  write(*, '(A28,X,F16.10)') 'Energy ECMD SCAN        = ',Energy_c_md_SCAN_mu_of_r               
 else if(mu_of_r_functional.EQ."basis_set_on_top_PBE")then
   provide Energy_c_md_PBE_mu_of_r  Energy_c_md_n_and_PBE_mu_of_r !Energy_c_md_on_top_SCAN_mu_of_r Energy_c_md_n_and_SCAN_mu_of_r
   provide Energy_c_md_on_top_PBE_mu_of_r Energy_c_md_n_and_on_top_PBE_mu_of_r !Energy_c_md_on_top_SCAN_mu_of_r Energy_c_md_n_and_on_top_SCAN_mu_of_r
   print*,''
   print*,'Corrections using Multi determinant mu'
   print*,''
   print*,'Functionals with UEG ontop pair density at large mu'
   print*,''
   write(*, '(A40,X,F16.10)') 'ECMD LDA        regular    spin dens = ',Energy_c_md_LDA_mu_of_r
   write(*, '(A40,X,F16.10)') 'ECMD LDA        effective  spin dens = ',Energy_c_md_n_and_LDA_mu_of_r
   write(*, '(A40,X,F16.10)') 'ECMD PBE        regular    spin dens = ',Energy_c_md_PBE_mu_of_r
   write(*, '(A40,X,F16.10)') 'ECMD PBE        effective  spin dens = ',Energy_c_md_n_and_PBE_mu_of_r
   write(*, '(A40,X,F16.10)') 'ECMD PBE        NO         spin dens = ',Energy_c_md_zero_spin_PBE_UEG_mu_of_r
!  write(*, '(A40,X,F16.10)') 'ECMD SCAN       regular    spin dens = ',Energy_c_md_SCAN_mu_of_r               
!  write(*, '(A40,X,F16.10)') 'ECMD SCAN       effective  spin dens = ',Energy_c_md_n_and_SCAN_mu_of_r         
!  write(*, '(A40,X,F16.10)') 'ECMD HOLPEG     regular    spin dens = ',Energy_c_md_holpeg_mu_of_r
   print*,''
   print*,'Functionals with extrapolated exact ontop based on current wave function '
   print*,''
   write(*, '(A40,X,F16.10)') 'ECMD PBE/ontop  regular    spin dens = ',Energy_c_md_on_top_PBE_mu_of_r
   write(*, '(A40,X,F16.10)') 'ECMD PBE/ontop  effective  spin dens = ',Energy_c_md_n_and_on_top_PBE_mu_of_r
   write(*, '(A40,X,F16.10)') 'ECMD PBE/ontop  NO         spin dens = ',Energy_c_md_no_spin_dens_and_on_top_PBE_mu_of_r
!  write(*, '(A40,X,F16.10)') 'ECMD SCAN/ontop regular    spin dens = ',Energy_c_md_on_top_SCAN_mu_of_r
!  write(*, '(A40,X,F16.10)') 'ECMD SCAN/ontop effective  spin dens = ',Energy_c_md_n_and_on_top_SCAN_mu_of_r
   print*,''
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
   write(*, '(A28,X,F16.10)') 'Energy ECMD LDA        = ',Energy_c_md_LDA_mu_of_r
!endif
  if(.true.)then
   write(*, '(A28,X,F16.10)') 'mu_average for basis set  = ',mu_average
  endif

end

