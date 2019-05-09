program two_body_dm
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
 read_wf = .True.
 touch read_wf
 call test_on_top
!call comp_test
!touch two_bod_alpha_beta_mo_physicist
 
!call routine_print
! call provide_everything
!call print_gamma
!call zero_gamma
!call comparaison_decomp_tensor
!call print_gamma
end

subroutine comp_test
 implicit none
 integer :: i
 double precision :: accu
 accu = 0.d0
 do i = 1, n_points_final_grid
  accu += final_weight_at_r_vector(i)
 enddo
 print*,'accu = ',accu
!provide two_bod_alpha_beta_mo_physicist 
!call comp_test2
 
end

subroutine comp_test2
 implicit none
 provide int_on_top_of_r_approx_svd

end

 subroutine zero_gamma 
 implicit none
 integer :: istate,i,j,k,l
 do istate = 1, N_states
!!aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa 
  do i = elec_alpha_num+1, mo_num
   do j = elec_beta_num+1, mo_num
    do k = elec_alpha_num+1, mo_num
     do l = 1, elec_beta_num 
      !                               1 2 1 2                                 1
      !                               1 2 2 
      two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) = 0d0 
     enddo
    enddo
   enddo
  enddo
 
 
  do i = elec_alpha_num+1, mo_num
   do j = elec_beta_num+1, mo_num
    do k = 1,elec_alpha_num
     do l = elec_beta_num+1, mo_num 
      !                               1 2 1 2                                 1
      !                               1 2 2 
      two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) = 0d0 
     enddo
    enddo
   enddo
  enddo
 
 
  do i = elec_alpha_num+1, mo_num
   do j = 1, elec_beta_num 
    do k = elec_alpha_num+1, mo_num 
     do l = elec_beta_num+1, mo_num 
      !                               1 2 1 2                                 1
      !                               1 2 2 
      two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) = 0d0 
     enddo
    enddo
   enddo
  enddo
 
  do i = 1,elec_alpha_num
   do j = elec_beta_num+1, mo_num 
    do k = elec_alpha_num+1, mo_num
     do l = elec_beta_num+1, mo_num
      !                               1 2 1 2                                 1
      two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) = 0d0
     enddo
    enddo
   enddo
  enddo
!!!!! tout le monde virtuel
 
 
  do i = elec_alpha_num+1, mo_num 
   do j = elec_beta_num+1, mo_num 
    do k = elec_alpha_num+1, mo_num
     do l = elec_beta_num+1, mo_num
      !                               1 2 1 2                                 1
      two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) = 0d0
     enddo
    enddo
   enddo
  enddo
!!!!!!!comptage HF/Virtuelle

  do i = 1, mo_num
   do j = elec_beta_num+1, mo_num
    do k =1, mo_num
     do l = elec_beta_num+1, mo_num
      !                               1 2 1 2                                 1
      two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) = 0d0
     enddo
    enddo
   enddo
  enddo

  do i = elec_alpha_num+1, mo_num
   do j = 1, mo_num
    do k = elec_alpha_num+1, mo_num
     do l = 1, mo_num
      !                               1 2 1 2                                 1
      two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) = 0d0
     enddo
    enddo
   enddo
  enddo




 enddo 


 touch two_bod_alpha_beta_mo_physicist

end


 subroutine print_gamma 
 implicit none
 integer :: istate,i,j,k,l
 print*,'*****************************'
 do istate = 1, N_states
!!aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa 
  do i = 1, mo_num
   do j = 1, mo_num
    do k = 1, mo_num
     do l = 1, mo_num
      if (dabs(two_bod_alpha_beta_mo_physicist(i,j,k,l,istate)) .gt. 1.d-15) then
       print*, 'i,j,k,l,Gamma =   ',i,' ',j,'  ',k, '  ',l,'  ',two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) 
      endif
     enddo
    enddo
   enddo
  enddo
 enddo


end



 subroutine comparaison_decomp_tensor
 implicit none

!double precision :: accu_tucker
!accu_tucker= E_cor_tot_normal_prov-integral_on_top_of_r_tucker(1)

 double precision :: accu_Manu
 accu_Manu= E_cor_tot_normal_prov-int_on_top_of_r_approx_svd(1)

 double precision :: accu_Manu_corr
 accu_Manu_corr= E_cor_tot_normal_prov-int_on_top_of_r_approx_svd_correl(1)

 double precision :: accu_Manu_hf_s_d
 accu_Manu_hf_s_d=E_cor_tot_normal_prov-int_on_top_of_r_approx_svd_hf_s_d(1)

 print*, '**************'
 !print*, 'Absolute error tucker          =', accu
 print*, 'Absolute error manual               =', accu_Manu
 print*, 'Absolute error manual correlation   =', accu_Manu_corr 
 print*, 'Absolute error manual HF S D        =', accu_Manu_hf_s_d
 print*, '**************'

 double precision :: accu_Manu_rela
 accu_Manu_rela = (E_cor_tot_normal_prov-int_on_top_of_r_approx_svd(1))/E_cor_tot_normal_prov

 double precision :: accu_Manu_corr_rela
 accu_Manu_corr_rela = (E_cor_tot_normal_prov-int_on_top_of_r_approx_svd_correl(1))/E_cor_tot_normal_prov

 double precision :: accu_Manu_hf_s_d_rela
 accu_Manu_hf_s_d_rela = (E_cor_tot_normal_prov-int_on_top_of_r_approx_svd_hf_s_d(1))/E_cor_tot_normal_prov

 print*, '**************'
 !print*, 'Absolute error tucker          =', accu
 print*, 'Relative error manual               =', accu_Manu_rela
 print*, 'Relative error manual correlation   =', accu_Manu_corr_rela 
 print*, 'Relative error manual HF S D        =', accu_Manu_hf_s_d_rela
 print*, '**************'

 print*, '**************'
 print*, 'E_cor_tot_normal_provider         =', E_cor_tot_normal_prov
 print*, 'E_cor_tot_manual_provider         =', int_on_top_of_r_approx_svd(1)
 !print*, 'E_cor_tot_tucker_provider         =', integral_on_top_of_r_tucker(1) 
 print*, '**************'
 print*, 'E_cor_tot_approx svd correlation  =',int_on_top_of_r_approx_svd_correl(1)
 print*, 'E_cor_tot_approx svd HF S D       =',int_on_top_of_r_approx_svd_hf_s_d(1)
 end



subroutine test_on_top
 implicit none
 integer :: i,j
 double precision :: accu, on_top_of_r_from_provider,core_inact_act_on_top_of_r_from_provider 
 double precision :: accuex,accuc
 accuex = 0.d0
 accuc  = 0.d0
 accu = 0.d0
 do i = 1, n_points_final_grid
  accuc += core_inact_act_on_top_of_r(i,1) * final_weight_at_r_vector(i) 
  accuex += on_top_of_r_vector(i,1) * final_weight_at_r_vector(i)
  accu += dabs(accuex - accuc) * final_weight_at_r_vector(i)
 enddo
 print*,'accuex = ',accuex
 print*,'accuc  = ',accuc
 print*,'accu   = ',accu

end
