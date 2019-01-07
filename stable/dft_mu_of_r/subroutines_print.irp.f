subroutine print_contribution_dft_mu_of_r
 implicit none

 print*,  '****************************************'
 print*, 'Functional used   = ',md_correlation_functional
 print*,  '****************************************'
 print*, 'mu_of_r_potential = ',mu_of_r_potential
 print*,  ' MR DFT energy with pure correlation part for the DFT '
 if(md_correlation_functional.EQ."basis_set_LDA")then
   print*, ''
   write(*, '(A28,X,F16.10)') 'Energy_c_md_mu_of_r_LDA   = ',Energy_c_md_mu_of_r_LDA
   print*, ''
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else if(md_correlation_functional.EQ."basis_set_on_top_PBE")then
   write(*, '(A28,X,F16.10)') 'Energy ECMD UEG        = ',Energy_c_md_on_top_PBE_mu_of_r_UEG
   write(*, '(A28,X,F16.10)') 'Energy ECMD NO UEG     = ',Energy_c_md_on_top_PBE_mu_of_r
   write(*, '(A28,X,F16.10)') 'Energy ECMD large mu(r)= ',Energy_c_md_on_top_mu_of_r
   print*,''
   write(*, '(A28,X,F16.10)') 'Energy ECMD LDA        = ',Energy_c_md_mu_of_r_LDA
 endif
  if(.true.)then
   write(*, '(A28,X,F16.10)') 'mu_average for basis set  = ',mu_average
  endif

end

