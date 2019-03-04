program DFT_Utils_ECMD
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
! call print_ecmd_var_energy
 call routine
end

subroutine routine
 implicit none
 print*,'Ecmd_pbe_n2_hf_aa = ',Ecmd_pbe_n2_hf_aa
 print*,'Ecmd_pbe_n2_hf_bb = ',Ecmd_pbe_n2_hf_bb


end
