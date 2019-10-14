program write_integrals_for_dft_ecmd_lda
 implicit none
 read_wf = .true.
 touch read_wf
 no_core_density = "no_core_dm"
 touch no_core_density
 mu_of_r_potential = "hf_valence_coallescence"
 touch mu_of_r_potential
 call truc_a_faire_lda
!call truc_a_faire_lda_old_school
 call truc_a_faire_pbeueg    
end


 subroutine truc_a_faire_lda_old_school 
 implicit none
 integer :: i,istate
 double precision :: r(3)
 double precision :: rhoa(n_states),rhob(n_states),mos_array(mo_num)
 double precision :: d_total_deltarho_rhoa,d_total_deltarho_rhob 
 double precision :: accu,accu_2,accu_3,mu 
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '


!!VErrrrifffff
!logical :: dospin
 double precision :: e_c,sr_vc_a,sr_vc_b
!dospin = .true. ! JT dospin have to be set to true for open shell
!poto =0.d0

 double precision :: threshold
 double precision :: homo
 threshold = 1d-15

 do istate = 1, n_states

  accu = 0.d0 
  accu_2 = 0.d0 
  accu_3 = 0.d0 
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   mu = mu_of_r_vector(i)

   rhoa(istate) = one_e_dm_alpha_at_r(i,istate)
   rhob(istate) = one_e_dm_beta_at_r(i,istate)
   call ec_lda_sr(mu,rhoa(istate),rhob(istate),e_c,sr_vc_a,sr_vc_b)
   call give_all_mos_at_r(r,mos_array)
   homo = mos_in_r_array_transp(i,elec_alpha_num)
   if(dabs(rhoa(istate)+rhob(istate)).lt.threshold) cycle 
   accu += (d_total_deltarho_rhoa(rhoa(istate),rhob(istate),mu)+d_total_deltarho_rhob(rhoa(istate),rhob(istate),mu)+sr_vc_a+sr_vc_b)*homo**2.d0*final_weight_at_r_vector(i)
   accu_2 += (d_total_deltarho_rhoa(rhoa(istate),rhob(istate),mu)+sr_vc_a)*homo**2.d0*final_weight_at_r_vector(i)
   accu_3 += sr_vc_a*homo**2.d0*final_weight_at_r_vector(i)
 
  enddo


  print*,'< HOMO | dE_c,md/dn_alpha+dE_c/dn_beta | HOMO >          =',accu,'/', potential_deltarho_ecmd_alpha_mo(elec_alpha_num,elec_alpha_num,istate)+ potential_e_c_lda_ecmd_alpha_mo(elec_alpha_num,elec_alpha_num,istate) + potential_deltarho_ecmd_beta_mo(elec_alpha_num,elec_alpha_num,istate) + potential_e_c_lda_ecmd_beta_mo(elec_alpha_num,elec_alpha_num,istate)
  print*,'< HOMO | dE_c,md/dn_alpha | HOMO >                       =',accu_2 ,'/', potential_deltarho_ecmd_alpha_mo(elec_alpha_num,elec_alpha_num,istate)+ potential_e_c_lda_ecmd_alpha_mo(elec_alpha_num,elec_alpha_num,istate)
  print*,'< HOMO | dE_c/dn_alpha | HOMO >                          =',accu_3,'/',potential_e_c_lda_ecmd_alpha_mo(elec_alpha_num,elec_alpha_num,istate)
 enddo
 end


 subroutine truc_a_faire_lda   
 implicit none
 integer :: istate
 double precision :: homo1,homo2,homo3
 double precision :: lumo1,lumo2,lumo3

 print*,' '
 print*,' '
 print*,'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
 print*,'\\\\\\\\\\\\\\\\\WARNING: only for close shell\\\\\\\\\\\\\\\\\'
 print*,'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'

 do istate = 1, n_states
  print*,'************************HOMO*************************************'
  homo1 = potential_deltarho_ecmd_alpha_mo(elec_alpha_num,elec_alpha_num,istate)+ potential_e_c_lda_ecmd_alpha_mo(elec_alpha_num,elec_alpha_num,istate) + potential_deltarho_ecmd_beta_mo(elec_alpha_num,elec_alpha_num,istate) + potential_e_c_lda_ecmd_beta_mo(elec_alpha_num,elec_alpha_num,istate)
  print*,'< HOMO | dE_c,md/dn_alpha+dE_c/dn_beta | HOMO >          =',homo1

  homo2 = potential_deltarho_ecmd_alpha_mo(elec_alpha_num,elec_alpha_num,istate)+ potential_e_c_lda_ecmd_alpha_mo(elec_alpha_num,elec_alpha_num,istate)
  print*,'< HOMO | dE_c,md/dn_alpha | HOMO >                       =',homo2

  homo3 = potential_e_c_lda_ecmd_alpha_mo(elec_alpha_num,elec_alpha_num,istate)
  print*,'< HOMO | dE_c/dn_alpha | HOMO >                          =',homo3
  print*,' '
  print*,' '
  print*,'***********************LUMO**************************************'
  lumo1 = potential_deltarho_ecmd_alpha_mo(elec_alpha_num+1,elec_alpha_num+1,istate)+ potential_e_c_lda_ecmd_alpha_mo(elec_alpha_num+1,elec_alpha_num+1,istate) + potential_deltarho_ecmd_beta_mo(elec_alpha_num+1,elec_alpha_num+1,istate) + potential_e_c_lda_ecmd_beta_mo(elec_alpha_num+1,elec_alpha_num+1,istate)
  print*,'< LUMO | dE_c,md/dn_alpha+dE_c/dn_beta | LUMO >          =',lumo1

  lumo2 = potential_deltarho_ecmd_alpha_mo(elec_alpha_num+1,elec_alpha_num+1,istate)+ potential_e_c_lda_ecmd_alpha_mo(elec_alpha_num+1,elec_alpha_num+1,istate)
  print*,'< LUMO | dE_c,md/dn_alpha | LUMO >                       =',lumo2

  lumo3 = potential_e_c_lda_ecmd_alpha_mo(elec_alpha_num+1,elec_alpha_num+1,istate)
  print*,'< LUMO | dE_c/dn_alpha | LUMO >                          =',lumo3
  print*,' '
  print*,' '
  print*,'LUMO-HOMO'
  print*,'< LUMO | dE_c,md/dn_alpha+dE_c/dn_beta | LUMO > - < HOMO | dE_c,md/dn_alpha+dE_c/dn_beta | HOMO >         =',lumo1 - homo1
  print*,'< LUMO | dE_c,md/dn_alpha | LUMO > - < HOMO | dE_c,md/dn_alpha | HOMO >                                   =',lumo2 - homo2
  print*,'< LUMO | dE_c/dn_alpha | LUMO > - < HOMO | dE_c/dn_alpha | HOMO >                                         =',lumo3 - homo3



 enddo
end


 subroutine truc_a_faire_pbeueg
 implicit none
 integer :: istate
 double precision :: homo1,homo2
 double precision :: lumo1,lumo2
 print*,' '
 print*,' '
 print*,'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
 print*,'\\\\\\\\\\\\\\\\\WARNING: only for close shell\\\\\\\\\\\\\\\\\'
 print*,'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'

 do istate = 1, n_states
  print*,'************************HOMO*************************************'
  homo1 = potential_c_alpha_mo_sr_pbe_ueg(elec_alpha_num,elec_alpha_num,istate) + potential_c_beta_mo_sr_pbe_ueg(elec_alpha_num,elec_alpha_num,istate)
  print*,'< HOMO | dE_c,md/dn_alpha+dE_c/dn_beta PBE UEG | HOMO >          =',homo1

  homo2 = potential_c_alpha_mo_sr_pbe_ueg(elec_alpha_num,elec_alpha_num,istate)
  print*,'< HOMO | dE_c,md/dn_alpha PBE UEG | HOMO >                       =',homo2
  print*,' '
  print*,' '
  print*,'***********************LUMO**************************************'
  lumo1 = potential_c_alpha_mo_sr_pbe_ueg(elec_alpha_num+1,elec_alpha_num+1,istate) + potential_c_beta_mo_sr_pbe_ueg(elec_alpha_num+1,elec_alpha_num+1,istate)
  print*,'< LUMO | dE_c,md/dn_alpha+dE_c/dn_beta PBE UEG | LUMO >          =',lumo1

  lumo2 = potential_c_alpha_mo_sr_pbe_ueg(elec_alpha_num+1,elec_alpha_num+1,istate)
  print*,'< LUMO | dE_c,md/dn_alpha PBE UEG | LUMO >                       =',lumo2
  print*,' '
  print*,' '
  print*,'LUMO-HOMO'
  print*,'< LUMO | dE_c,md/dn_alpha+dE_c/dn_beta PBE UEG | LUMO > - < HOMO | dE_c,md/dn_alpha+dE_c/dn_beta PBE UEG | HOMO >   =',lumo1-homo1
  print*,'< LUMO | dE_c,md/dn_alpha PBE UEG | LUMO > - < HOMO | dE_c,md/dn_alpha PBE UEG | HOMO >                             =',lumo2-homo2

 enddo
end
