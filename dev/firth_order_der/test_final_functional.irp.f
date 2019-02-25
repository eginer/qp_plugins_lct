program pouet
 read_wf = .True.
 touch read_wf
 call test_ecmd_alpha 
end

 subroutine test_ecmd_alpha
 implicit none
 provide Energy_c_LDA_mu_of_r 
 provide e_c_md_mur_ab_LDA_a
 provide e_c_md_mur_ab_LDA_b
 provide e_c_md_mur_aa_LDA_a
 provide e_c_md_mur_bb_LDA_b
 print*,'****************************************'
 print*,' ecmd mu(r) ab LDA ab   =', Energy_c_LDA_mu_of_r
 print*,' ecmd mu(r) ab LDA a    =', e_c_md_mur_ab_LDA_a 
 print*,' ecmd mu(r) ab LDA b    =', e_c_md_mur_ab_LDA_b
 print*,' ecmd mu(r) aa LDA a    =', e_c_md_mur_aa_LDA_a
 print*,' ecmd mu(r) bb LDA b    =', e_c_md_mur_bb_LDA_b
 print*,'****************************************'
 print*,' '
 print*,'****************************************'
 print*,' mu average ab       =', mu_average
 print*,' mu average aa       =', mu_average_aa
 print*,' mu average bb       =', mu_average_bb
 print*,'****************************************'
 print*,' '
 print*,'****************************************'
 print*,' ecmd mu(r) ab LDA ab                                             =', Energy_c_LDA_mu_of_r
 print*,' ecmd mu(r) ab LDA ab - ecmd mu(r) ab LDA a - ecmd mu(r) ab LDA b =', Energy_c_LDA_mu_of_r - e_c_md_mur_ab_LDA_a - e_c_md_mur_ab_LDA_b
 print*,' ecmd mu(r) ab LDA ab new                                         =', Energy_c_LDA_mu_of_r - e_c_md_mur_ab_LDA_a - e_c_md_mur_ab_LDA_b + e_c_md_mur_aa_LDA_a + e_c_md_mur_bb_LDA_b

 end


