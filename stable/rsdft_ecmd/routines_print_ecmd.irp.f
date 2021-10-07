subroutine print_energy_ecmd
 implicit none
 double precision :: var_e
 var_e = psi_energy_with_nucl_rep(1)
 print*,'<Psi | H | Psi> = ',var_e
 print*,'*******'
 print*,'Energy_c_md_LDA = ',Energy_c_md_LDA
 print*,'E tot LDA       = ',Energy_c_md_LDA + var_e
 print*,'*******'
 print*,'ecmd_pbe_on_top = ',ecmd_pbe_on_top_at_mu
 print*,'E tot ON-TOP    = ',ecmd_pbe_on_top_at_mu+ var_e
 print*,'*******'
 print*,'ecmd_pbe_ueg    = ',ecmd_pbe_ueg_prov
 print*,'E tot PBE-UEG   = ',ecmd_pbe_ueg_prov+ var_e

end
