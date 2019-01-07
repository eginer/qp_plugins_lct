program DFT_Utils_mu_of_r
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf 
  call routine
end

subroutine routine
 implicit none
 print*,'Energy_c_md_mu_of_r_LDA        = ',Energy_c_md_mu_of_r_LDA
 print*,'Energy_c_md_on_top_PBE_mu_of_r = ',Energy_c_md_on_top_PBE_mu_of_r
 print*,'mu_average                     = ',mu_average



end
