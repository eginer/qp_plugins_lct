program write_integrals_for_dft_ecmd_lda
 implicit none
 read_wf = .true.
 touch read_wf
 no_core_density = "no_core_dm"
 touch no_core_density
 call truc_a_faire    
end


 subroutine truc_a_faire 
 implicit none
 integer :: i
 double precision :: r(3)
 double precision :: rhoa(n_states),rhob(n_states),mos_array(mo_num)
 double precision :: d_total_deltarho_rhoa,d_total_deltarho_rhob 
 double precision :: accu,istate,mu 
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '


!!VErrrrifffff
!logical :: dospin
!double precision :: ec(N_states),poto(N_states)
!dospin = .true. ! JT dospin have to be set to true for open shell
!poto =0.d0

 double precision :: threshold
 threshold = 1d-15

 do istate = 1, n_states

  accu = 0.d0 
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   mu = mu_of_r_vector(i)

   call dm_dft_alpha_beta_at_r(r,rhoa(istate),rhob(istate))
   call give_all_mos_at_r(r,mos_array)
   if(dabs(rhoa(istate)+rhob(istate)).lt.threshold) cycle 
   accu += (d_total_deltarho_rhoa(rhoa(istate),rhob(istate),mu)+d_total_deltarho_rhob(rhoa(istate),rhob(istate),mu))*mos_array(elec_alpha_num)**2.d0*final_weight_at_r_vector(i)

  !print*,d_total_deltarho_rhoa(rhoa(istate),rhob(istate),mu),mos_array(elec_alpha_num)**2.d0 
  !!VErrrrifffff
  !call ESRC_MD_LDAERF (mu,rhoa(istate),rhob(istate),dospin,ec(istate))
  !if(isnan(ec(istate)))then
  ! print*,'ec is nan'
  ! stop
  !endif
  !poto(istate) += ec(istate)*final_weight_at_r_vector(i)

 
  enddo
  print*,'< HOMO | dE_c/dn | HOMO >    =',accu
 !print*,'Vrai Ecmd LDA                =',Energy_c_md_LDA_mu_of_r(istate)
 !print*,'Ecmd LDA bart                =',poto(istate)
 enddo
end
