program print_eff_basis_potential
  implicit none
  BEGIN_DOC
! Compute <Psi | V^B |Psi> 
  END_DOC
  double precision :: potential_tot_alpha, potential_tot_beta
  integer :: i, j, ii, jj, nst 

  provide pot_basis_alpha_mo_basis_pbe_ueg
  provide pot_basis_beta_mo_basis_pbe_ueg
  provide data_one_e_dm_alpha_mo
  provide data_one_e_dm_beta_mo
  nst=1 ! state

  potential_tot_alpha = 0.d0
  potential_tot_beta = 0.d0

  do ii=1, n_act_orb
   i=list_act(ii)
   do jj=1, n_act_orb
    j=list_act(jj)
    potential_tot_alpha+= data_one_e_dm_alpha_mo(i,j,nst)*pot_basis_alpha_mo_basis_pbe_ueg(ii,jj,nst)
    potential_tot_beta+= data_one_e_dm_beta_mo(i,j,nst)*pot_basis_beta_mo_basis_pbe_ueg(ii,jj,nst)
   enddo
  enddo

 print*,'POTENTIAL ALPHA =', potential_tot_alpha
 print*,'POTENTIAL BETA =', potential_tot_beta
 print*,'POTENTIAL TOT =', potential_tot_alpha + potential_tot_beta

end
