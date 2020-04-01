program pouet
 read_wf = .True.
 touch read_wf
 !call test_delta      
 !call test_delta_derivative
 !call test_delta_derivative_xi
 !call test_delta_total_derivative 
 !call test_delta_total_derivative_bis 
 !call test_epsilon_ueg
 !call test_f_pbeueg
 !call test_total_derivative_beta
 !call test_total_derivative_bis_PBE_UEG
  call test_pot_pbe_ueg
end

!subroutine test_provider
!implicit none
!integer :: i,j
!double precision :: accu_1,accu_2,accu_3,accu_4
!accu_1 = 0d0
!accu_1 = 0d0
!accu_3 = 0d0
!accu_4 = 0d0
!do i = 1, ao_num
! do j = 1, ao_num
!  accu_1 += dabs(potential_deltarho_ecmd_alpha_ao(j,i,1) - potential_deltarho_ecmd_alpha_ao_2(j,i,1))
!  accu_2 += dabs(potential_deltarho_ecmd_beta_ao(j,i,1) - potential_deltarho_ecmd_beta_ao_2(j,i,1))
!  accu_3 += dabs(potential_e_c_lda_ecmd_alpha_ao(j,i,1) - potential_e_c_lda_ecmd_alpha_ao_2(j,i,1))
!  accu_4 += dabs(potential_e_c_lda_ecmd_beta_ao(j,i,1) - potential_e_c_lda_ecmd_beta_ao_2(j,i,1))
! enddo
!enddo
!print*,'dn=n*',accu_1,accu_2,accu_3,accu_4
!end



subroutine test_delta_derivative
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rs,rsplus,rsmoins,drs,accu
 double precision :: d_delta_barth_bis,delta_barth,d_delta_barth,delta_plus,delta_moins
 double precision :: aos_array(ao_num)
 double precision :: dn,rhot,rhos,rhoa,rhob,xi,mu,pi,denominator_delta,wignerseitz_radius,d_wignerseitz_radius

!!!!!!!!!!!delta2!!!!!!!!
 double precision :: accu_2,d_delta_barth_bis_2,delta_2,delta_plus_delta_2,delta_moins_delta_2,d_1st_deltaterm
!!!!!!!!!!!delta3!!!!!!!!
 double precision :: accu_3,d_delta_barth_bis_3,delta_3,delta_plus_delta_3,delta_moins_delta_3,d_2nd_deltaterm 
!!!!!!!!!!!delta4!!!!!!!!
 double precision :: accu_4,d_delta_barth_bis_4,delta_4,delta_plus_delta_4,delta_moins_delta_4,d_3rd_deltaterm
!!!!!!!!!!!delta5!!!!!!!!
 double precision :: accu_5,d_delta_barth_bis_5,delta_5,delta_plus_delta_5,delta_moins_delta_5,d_4th_deltaterm
!!!!!!!!!!!delta6!!!!!!!!
 double precision :: accu_6,d_delta_barth_bis_6,delta_6,delta_plus_delta_6,delta_moins_delta_6,d_5th_deltaterm


 pi=dacos(-1.d0)
 mu=mu_erf

 do n = 1, 16
!dn = 10d0**(-n)
 accu = 0d0
 accu_2= 0d0
 accu_3 = 0d0
 accu_5 = 0d0
 accu_4 = 0d0
 accu_6 = 0d0
  do j = 1,nucl_num
   do k = 1, n_points_radial_grid -1
    do l = 1 , n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call dm_dft_alpha_beta_and_all_aos_at_r(r,rhoa,rhob,aos_array)
     rhot = rhoa + rhob
     rhos = rhoa - rhob
     
     dn = rhot*10d0**(-n)

     xi = (rhoa-rhob)/(rhoa+rhob)
     rs = wignerseitz_radius(rhot) 
     drs= d_wignerseitz_radius(rhot)

     rsplus = wignerseitz_radius(rhot+dn)
     rsmoins = wignerseitz_radius(rhot-dn) 

!!!!!!!!!!!!!!!!!testsss delta2!!!!!!!!!
     
    delta_plus_delta_2 = delta_2(rsplus)*mu**2/denominator_delta(rsplus,xi,mu)
    delta_moins_delta_2 = delta_2(rsmoins)*mu**2/denominator_delta(rsmoins,xi,mu)
     
    d_delta_barth_bis_2= (delta_plus_delta_2-delta_moins_delta_2)/(2.d0 * dn)
 
    accu_2 += dabs(d_delta_barth_bis_2 - d_1st_deltaterm(rs,mu,xi)*drs)*final_weight_at_r(l,k,j)

!!!!!!!!!!!!!!!!!testsss delta3!!!!!!!!!

    delta_plus_delta_3 = delta_3(rsplus,xi)*mu**3/denominator_delta(rsplus,xi,mu)
    delta_moins_delta_3= delta_3(rsmoins,xi)*mu**3/denominator_delta(rsmoins,xi,mu)


    d_delta_barth_bis_3= (delta_plus_delta_3-delta_moins_delta_3)/(2.d0 * dn)

    accu_3 += dabs(d_delta_barth_bis_3-d_2nd_deltaterm(rs,mu,xi)*drs)*final_weight_at_r(l,k,j)


!!!!!!!!!!!!!!!!!TEST delta_4!!!!!!!!!
 
    delta_plus_delta_4 = delta_4(rsplus,xi)*mu**4/denominator_delta(rsplus,xi,mu)
    delta_moins_delta_4= delta_4(rsmoins,xi)*mu**4/denominator_delta(rsmoins,xi,mu)
         
 
    d_delta_barth_bis_4= (delta_plus_delta_4-delta_moins_delta_4)/(2.d0 * dn)
 
    accu_4 += dabs(d_delta_barth_bis_4-d_3rd_deltaterm(rs,mu,xi)*drs)*final_weight_at_r(l,k,j)



!!!!!!!!!!!!!!!!!TEST delta_5!!!!!!!!!

    delta_plus_delta_5 = delta_5(rsplus,xi)*mu**5/denominator_delta(rsplus,xi,mu)
    delta_moins_delta_5= delta_5(rsmoins,xi)*mu**5/denominator_delta(rsmoins,xi,mu)



    d_delta_barth_bis_5= (delta_plus_delta_5-delta_moins_delta_5)/(2.d0 * dn)

    accu_5 += dabs(d_delta_barth_bis_5-d_4th_deltaterm(rs,mu,xi)*drs)*final_weight_at_r(l,k,j)

!!!!!!!!!!!!!!!!!TEST delta_6!!!!!!!!!

    delta_plus_delta_6 = delta_6(rsplus,xi)*mu**6/denominator_delta(rsplus,xi,mu)
    delta_moins_delta_6= delta_6(rsmoins,xi)*mu**6/denominator_delta(rsmoins,xi,mu)
 
    d_delta_barth_bis_6= (delta_plus_delta_6-delta_moins_delta_6)/(2.d0 * dn)
 
    accu_6 += dabs(d_delta_barth_bis_6-d_5th_deltaterm(rs,mu,xi)*drs)*final_weight_at_r(l,k,j)
 
!!!!!!!!!!!!!!!!!TOTAL TEST!!!!!!!!!

    delta_plus = delta_barth(rsplus,xi,mu)
    delta_moins = delta_barth(rsmoins,xi,mu)

    d_delta_barth_bis = (delta_plus - delta_moins) /(2.d0 * dn)

    accu += dabs(d_delta_barth_bis - d_delta_barth(rs,xi,mu)*drs)*rhot*final_weight_at_r(l,k,j)

    enddo
   enddo
  enddo
  !print*,'dn=n*',10d0**(-n),accu_2,accu_3,accu_4,accu_5,accu_6
  print*,'dn=n*',10d0**(-n),accu
 enddo

end



subroutine test_delta_derivative_xi
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rs,xiplus,ximoins,drs,accu,dxi_dna,drhoa
 double precision :: d_delta_barth_bis,delta_barth,d_xi_delta_barth,delta_plus,delta_moins
 double precision :: aos_array(ao_num)
 double precision :: dxi,rhot,rhos,rhoa,rhob,xi,mu,pi,denominator_delta,wignerseitz_radius,d_wignerseitz_radius

!!!!!!!!!!!delta2!!!!!!!!
 double precision :: accu_2,d_delta_barth_bis_2,delta_2,delta_plus_delta_2,delta_moins_delta_2,d_xi_1st_deltaterm
!!!!!!!!!!!delta3!!!!!!!!
 double precision :: accu_3,d_delta_barth_bis_3,delta_3,delta_plus_delta_3,delta_moins_delta_3,d_xi_2nd_deltaterm
!!!!!!!!!!!delta4!!!!!!!!
 double precision :: accu_4,d_delta_barth_bis_4,delta_4,delta_plus_delta_4,delta_moins_delta_4,d_xi_3rd_deltaterm
!!!!!!!!!!!delta5!!!!!!!!
 double precision :: accu_5,d_delta_barth_bis_5,delta_5,delta_plus_delta_5,delta_moins_delta_5,d_xi_4th_deltaterm
!!!!!!!!!!!delta6!!!!!!!!
 double precision :: accu_6,d_delta_barth_bis_6,delta_6,delta_plus_delta_6,delta_moins_delta_6,d_xi_5th_deltaterm


 pi=dacos(-1.d0)
 mu=mu_erf

 do n = 1, 16
!dn = 10d0**(-n)
 accu = 0d0
 accu_2= 0d0
 accu_3 = 0d0
 accu_5 = 0d0
 accu_4 = 0d0
 accu_6 = 0d0
  do j = 1,nucl_num
   do k = 1, n_points_radial_grid -1
    do l = 1 , n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call dm_dft_alpha_beta_and_all_aos_at_r(r,rhoa,rhob,aos_array)
     rhot = rhoa + rhob
     rhos = rhoa - rhob

     drhoa = rhoa*10d0**(-n)

     xi = (rhoa-rhob)/(rhoa+rhob)
   
     dxi_dna=2d0*rhob/(rhoa+rhob)**2      
  
     rs = wignerseitz_radius(rhot)

     xiplus = (rhoa+drhoa-rhob)/(rhoa+drhoa+rhob)
     ximoins = (rhoa-drhoa-rhob)/(rhoa-drhoa+rhob) 

!!!!!!!!!!!!!!!!!testsss delta2!!!!!!!!!
    delta_plus_delta_2 = delta_2(rs)*mu**2/denominator_delta(rs,xiplus,mu)
    delta_moins_delta_2 = delta_2(rs)*mu**2/denominator_delta(rs,ximoins,mu)

    d_delta_barth_bis_2= (delta_plus_delta_2-delta_moins_delta_2)/(2.d0 * drhoa)
 
    accu_2 += dabs(d_delta_barth_bis_2 - d_xi_1st_deltaterm(rs,mu,xi)*dxi_dna)*final_weight_at_r(l,k,j)
!!!!!!!!!!!!!!!!!testsss delta3!!!!!!!!!
    delta_plus_delta_3 = delta_3(rs,xiplus)*mu**3/denominator_delta(rs,xiplus,mu)
    delta_moins_delta_3= delta_3(rs,ximoins)*mu**3/denominator_delta(rs,ximoins,mu)

    d_delta_barth_bis_3= (delta_plus_delta_3-delta_moins_delta_3)/(2.d0 * drhoa)

    accu_3 += dabs(d_delta_barth_bis_3-d_xi_2nd_deltaterm(rs,mu,xi)*dxi_dna)*final_weight_at_r(l,k,j)
!!!!!!!!!!!!!!!!!TEST delta_4!!!!!!!!!
    delta_plus_delta_4 = delta_4(rs,xiplus)*mu**4/denominator_delta(rs,xiplus,mu)
    delta_moins_delta_4= delta_4(rs,ximoins)*mu**4/denominator_delta(rs,ximoins,mu)

    d_delta_barth_bis_4= (delta_plus_delta_4-delta_moins_delta_4)/(2.d0 * drhoa)

    accu_4 += dabs(d_delta_barth_bis_4-d_xi_3rd_deltaterm(rs,mu,xi)*dxi_dna)*final_weight_at_r(l,k,j)
!!!!!!!!!!!!!!!!!TEST delta_5!!!!!!!!!
    delta_plus_delta_5 = delta_5(rs,xiplus)*mu**5/denominator_delta(rs,xiplus,mu)
    delta_moins_delta_5= delta_5(rs,ximoins)*mu**5/denominator_delta(rs,ximoins,mu)

    d_delta_barth_bis_5= (delta_plus_delta_5-delta_moins_delta_5)/(2.d0 * drhoa)

    accu_5 += dabs(d_delta_barth_bis_5-d_xi_4th_deltaterm(rs,mu,xi)*dxi_dna)*final_weight_at_r(l,k,j)
!!!!!!!!!!!!!!!!!TEST delta_6!!!!!!!!!
    delta_plus_delta_6 = delta_6(rs,xiplus)*mu**6/denominator_delta(rs,xiplus,mu)
    delta_moins_delta_6= delta_6(rs,ximoins)*mu**6/denominator_delta(rs,ximoins,mu)

    d_delta_barth_bis_6= (delta_plus_delta_6-delta_moins_delta_6)/(2.d0 * drhoa)

    accu_6 += dabs(d_delta_barth_bis_6-d_xi_5th_deltaterm(rs,mu,xi)*dxi_dna)*final_weight_at_r(l,k,j)
!!!!!!!!!!!!!!!!!TOTAL TEST!!!!!!!!!
    delta_plus = delta_barth(rs,xiplus,mu)
    delta_moins = delta_barth(rs,ximoins,mu)

    d_delta_barth_bis = (delta_plus - delta_moins) /(2.d0 * drhoa)

    accu += dabs(d_delta_barth_bis - d_xi_delta_barth(rs,xi,mu)*dxi_dna)*final_weight_at_r(l,k,j)
    enddo
   enddo
  enddo
  !print*,'dxi=xi*',10d0**(-n),accu_2,accu_3,accu_4,accu_5,accu_6
  print*,'dn=n*',10d0**(-n),accu
 enddo

end



subroutine test_delta_total_derivative
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rs,xiplus,ximoins,accu,accu_2,dxi_dna,drhoa,drs,rhotplus_a,rhotmoins_a,rsmoins,rsplus
 double precision :: d_delta_barth_bis,delta_barth,d_xi_delta_barth,delta_plus,delta_moins,d_delta_barth
 double precision :: aos_array(ao_num)
 double precision :: dxi,rhot,rhos,rhoa,rhob,xi,mu,pi,denominator_delta,wignerseitz_radius,d_wignerseitz_radius
 double precision :: delta_moins_b,delta_plus_b,d_delta_barth_bis_b,rsplus_b,rsmoins_b,rhotplus_b,rhotmoins_b,drhob,dxi_dnb,xiplus_b,ximoins_b 
 pi=dacos(-1.d0)
 mu=mu_erf
 
 do n = 1, 16
!dn = 10d0**(-n)
 accu = 0d0
 accu_2= 0d0
  do j = 1,nucl_num
   do k = 1, n_points_radial_grid -1
    do l = 1 , n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call dm_dft_alpha_beta_and_all_aos_at_r(r,rhoa,rhob,aos_array)
     rhot = rhoa + rhob
     rhos = rhoa - rhob

     drhoa = rhoa*10d0**(-n)
     drhob = rhob*10d0**(-n)
    
     xi = (rhoa-rhob)/(rhoa+rhob)

     dxi_dna=2d0*rhob/(rhoa+rhob)**2
     dxi_dnb=2d0*rhoa/(rhoa+rhob)**2
     
     rs = wignerseitz_radius(rhot)

     drs= d_wignerseitz_radius(rhot)

     rhotplus_a = rhoa+drhoa + rhob
     rhotmoins_a= rhoa-drhoa + rhob

     rsplus = wignerseitz_radius(rhotplus_a)
     rsmoins = wignerseitz_radius(rhotmoins_a)

     rhotplus_b = rhoa+drhob + rhob
     rhotmoins_b= rhoa-drhob + rhob

     rsplus_b = wignerseitz_radius(rhotplus_b)
     rsmoins_b = wignerseitz_radius(rhotmoins_b)


     xiplus = (rhoa+drhoa-rhob)/(rhoa+drhoa+rhob)
     ximoins = (rhoa-drhoa-rhob)/(rhoa-drhoa+rhob)

     xiplus_b = (rhoa-drhob-rhob)/(rhoa+drhob+rhob)
     ximoins_b = (rhoa+drhob-rhob)/(rhoa-drhob+rhob)

!!!!!!!!!!!!!!!!!ALPHA!!!!!!!!!
    delta_plus = delta_barth(rsplus,xiplus,mu)
    delta_moins = delta_barth(rsmoins,ximoins,mu)

    d_delta_barth_bis = (delta_plus - delta_moins) /(2.d0 * drhoa)

    accu += dabs(d_delta_barth_bis - (d_delta_barth(rs,xi,mu)*drs+d_xi_delta_barth(rs,xi,mu)*dxi_dna))*final_weight_at_r(l,k,j)

!!!!!!!!!!!!!!!!!BETA!!!!!!!!!

    delta_plus_b = delta_barth(rsplus_b,xiplus_b,mu)
    delta_moins_b = delta_barth(rsmoins_b,ximoins_b,mu)

    d_delta_barth_bis_b = (delta_plus_b - delta_moins_b) /(2.d0 * drhob)

    accu_2 += dabs(d_delta_barth_bis - (d_delta_barth(rs,xi,mu)*drs+d_xi_delta_barth(rs,xi,mu)*dxi_dnb))*final_weight_at_r(l,k,j)

    enddo
   enddo
  enddo
  !print*,'dxi=xi*',10d0**(-n),accu_2,accu_3,accu_4,accu_5,accu_6
  print*,'dn=n*',10d0**(-n),accu,accu_2
 enddo

end


subroutine test_delta_total_derivative_bis
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rs,xiplus,ximoins,accu,accu_2,dxi_dna,drhoa,drs,rhotplus_a,rhotmoins_a,rsmoins,rsplus
 double precision :: d_delta_barth_bis,delta_barth,d_xi_delta_barth,delta_plus,delta_moins,d_delta_barth
 double precision :: aos_array(ao_num)
 double precision :: dxi,rhot,rhos,rhoa,rhob,xi,mu,pi,denominator_delta,wignerseitz_radius,d_wignerseitz_radius
 double precision :: d_total_deltarho_rhob,d_total_deltarho_rhoa,d_xi_rhoa,d_xi_rhob,delta_moins_b,delta_plus_b,d_delta_barth_bis_b,rsplus_b,rsmoins_b,rhotplus_b,rhotmoins_b,drhob,dxi_dnb,xiplus_b,ximoins_b
 pi=dacos(-1.d0)
 mu=mu_erf

 do n = 1, 16
!dn = 10d0**(-n)
 accu = 0d0
 accu_2= 0d0
  do j = 1,nucl_num
   do k = 1, n_points_radial_grid -1
    do l = 1 , n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call dm_dft_alpha_beta_and_all_aos_at_r(r,rhoa,rhob,aos_array)
     rhot = rhoa + rhob
     rhos = rhoa - rhob

     drhoa = rhoa*10d0**(-n)
     drhob = rhob*10d0**(-n)

     xi = (rhoa-rhob)/(rhoa+rhob)

     rhotplus_a = rhoa+drhoa + rhob
     rhotmoins_a= rhoa-drhoa + rhob

     rsplus = wignerseitz_radius(rhotplus_a)
     rsmoins = wignerseitz_radius(rhotmoins_a)

     rhotplus_b = rhoa+drhob + rhob
     rhotmoins_b= rhoa-drhob + rhob

     rsplus_b = wignerseitz_radius(rhotplus_b)
     rsmoins_b = wignerseitz_radius(rhotmoins_b)


     xiplus = (rhoa+drhoa-rhob)/(rhoa+drhoa+rhob)
     ximoins = (rhoa-drhoa-rhob)/(rhoa-drhoa+rhob)

     xiplus_b = (rhoa-drhob-rhob)/(rhoa+drhob+rhob)
     ximoins_b = (rhoa+drhob-rhob)/(rhoa-drhob+rhob)


!!!!!!!!!!!!!!!!!ALPHA!!!!!!!!!
    delta_plus = delta_barth(rsplus,xiplus,mu)*(rhot+drhoa)
    delta_moins = delta_barth(rsmoins,ximoins,mu)*(rhot-drhoa)

    d_delta_barth_bis = (delta_plus - delta_moins) /(2.d0 * drhoa)

    accu += dabs(d_delta_barth_bis - d_total_deltarho_rhoa(rhoa,rhob,mu))*final_weight_at_r(l,k,j)

!!!!!!!!!!!!!!!!!BETA!!!!!!!!!

    delta_plus_b = delta_barth(rsplus_b,xiplus_b,mu)*(rhot+drhob)
    delta_moins_b = delta_barth(rsmoins_b,ximoins_b,mu)*(rhot-drhob)

    d_delta_barth_bis_b = (delta_plus_b - delta_moins_b) /(2.d0 * drhob)

    accu_2 += dabs(d_delta_barth_bis - d_total_deltarho_rhob(rhoa,rhob,mu))*final_weight_at_r(l,k,j)

    enddo
   enddo
  enddo
  !print*,'dxi=xi*',10d0**(-n),accu_2,accu_3,accu_4,accu_5,accu_6
  print*,'dn=n*',10d0**(-n),accu,accu_2
 enddo

end

subroutine test_epsilon_ueg 
 implicit none
 integer :: i,istate
 double precision :: mu,Energy_julien(N_states),Energy_barth_g0_new(N_states),Energy_barth_g0_old(N_states),r(3),rhoa,rhob,weight,Energy_barth_test(N_states)
 double precision :: aos_array(ao_num),grad_rho_a(3,N_states),grad_rho_b(3,N_states),grad_aos_array(3,ao_num) 
 double precision :: eps_c_md_PBE(N_states),eps_c_md_PBE_2(N_states),epsilon_c_pbeueg_bart(N_states),beta_test(N_states),beta_test2(N_states),epsilon_c_pbeueg_bartg0f(N_states)
 mu = mu_erf
 Energy_julien = 0.d0
 Energy_barth_g0_new = 0.d0
 Energy_barth_g0_old = 0.d0
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rhoa,rhob, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
  weight=final_weight_at_r_vector(i)
  call eps_c_md_PBE_at_grid_pt(mu,i,eps_c_md_PBE)
 !call give_epsilon_pbe_effective_spin_dens_provider_barth(mu,i,eps_c_md_PBE_2)
  call give_Ec_pbeueg_test(mu,rhoa,rhob,grad_rho_a,grad_rho_b,epsilon_c_pbeueg_bart,epsilon_c_pbeueg_bartg0f,beta_test,beta_test2) 
  do istate = 1, N_states
   Energy_julien(istate) += eps_c_md_PBE(istate) * weight
   Energy_barth_g0_new(istate) += epsilon_c_pbeueg_bart(istate) * weight
   Energy_barth_g0_old(istate) += epsilon_c_pbeueg_bartg0f(istate) * weight
  enddo
 enddo

 print*,'*****************$******************'
 print*,'ECMD PBE UEG REF     =',Energy_julien
 print*,'ECMD PBE UEG g0f JUL =',Energy_barth_g0_new
 print*,'ECMD PBE UEG g0f OLD =',Energy_barth_g0_old
 print*,'*****************$******************'
end


subroutine test_f_pbeueg
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),accu,accu_2
 double precision :: aos_array(ao_num),grad_rho_a(3,N_states),grad_rho_b(3,N_states),grad_aos_array(3,ao_num) 
 double precision :: rhot,rhos,rhoa,rhob,rhoa_plus,rhoa_moins,drhoa,drhob,pi,mu
 double precision :: threshold, weight,accu_dens

 pi=dacos(-1.d0)
 mu=mu_erf
 threshold = 1d-10

 do n = 1, 16
 accu = 0d0
 accu_2= 0d0
 accu_dens=0.d0
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight=final_weight_at_r_vector(i)
   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rhoa,rhob, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
  !if(dabs(rhob).lt.threshold) cycle 
   rhot = rhoa + rhob
 ! if(dabs(rhot).lt.threshold) cycle 
   rhos = rhoa - rhob
 
   drhoa = rhoa*10d0**(-n)
   drhob = rhob*10d0**(-n)
   rhoa_plus = rhoa + drhoa
   rhoa_moins = rhoa - drhoa
   

!!!!!!!!!!!!!!!!!!ALPHA!!!!!!!!!
   double precision :: xi_plus,xi_moins,rs_plus,rs_moins,f_plus,f_moins,d_f_pbeueg,d_f_pbeueg_ana
   double precision :: wignerseitz_radius,f_pbeueg,d_f_pbeueg_rhoa

   xi_plus= (rhoa+ drhoa -rhob)/(rhoa+rhob)
   xi_moins= (rhoa- drhoa -rhob)/(rhoa+rhob)
   
   rs_plus = wignerseitz_radius(rhot+drhoa)
   rs_moins = wignerseitz_radius(rhot-drhoa)
   
   f_plus =  f_pbeueg(rs_plus,xi_plus)
   f_moins =  f_pbeueg(rs_moins,xi_moins)
   
   d_f_pbeueg = (f_plus - f_moins)/(2.d0 * drhoa) 

   d_f_pbeueg_ana = d_f_pbeueg_rhoa(rhoa,rhob)

  !print*,'d_ec1 num, dec_ ana         =',d_Ec_barth_num_1,d_ec_pbeueg_rhoa(1)
   accu += dabs(d_f_pbeueg_ana - d_f_pbeueg)*weight
  !accu_2 += dabs(d_beta_rhoa(1) - d_beta_num_2)*weight 
   accu_dens += rhot*weight
!!!!!!!!!!!!!!!!BETA!!!!!!!!!
 
 
  enddo
  print*,'dn=n*',10d0**(-n),accu,accu_dens
 enddo

end

subroutine test_total_derivative_beta
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),accu,accu_2
 double precision :: aos_array(ao_num),grad_rho_a(3,N_states),grad_rho_b(3,N_states),grad_aos_array(3,ao_num) 
 double precision :: grad_rho_a_plus(3,N_states),grad_rho_a_moins(3,N_states)
 double precision :: rhot,rhos,rhoa,rhob,rhoa_plus,rhoa_moins,drhoa,drhob,pi,mu
 double precision :: beta_plus,beta_moins,threshold
 double precision :: weight

 pi=dacos(-1.d0)
 mu=mu_erf
 threshold = 1d-10

 double precision :: d_beta_rhoa(N_states),d_beta_rhob(N_states)
 double precision :: accu_dens_test_bg,accu_dens_test_looser 
 double precision :: accu_looser,accu_BG

!do n = 1, 16
 do n = 3, 3 
 accu = 0d0
 accu_2= 0d0
 accu_looser = 0.d0
 accu_BG = 0.d0
 accu_dens_test_bg = 0.d0
 accu_dens_test_looser = 0.d0
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight=final_weight_at_r_vector(i)
   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rhoa,rhob, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
 
  !if(dabs(rhob).lt.threshold) cycle 
   rhot = rhoa + rhob
  !if(dabs(rhot).lt.threshold) cycle 
   rhos = rhoa - rhob
 
   drhoa = rhoa*10d0**(-n)
   drhob = rhob*10d0**(-n)
   rhoa_plus = rhoa + drhoa
   rhoa_moins = rhoa - drhoa
 
  !call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rhoa_plus,rhob, grad_rho_a_plus, grad_rho_b, aos_array, grad_aos_array)
  !call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rhoa_moins,rhob, grad_rho_a_moins, grad_rho_b, aos_array, grad_aos_array)

!!!!!!!!!!!!!!!!!!ALPHA!!!!!!!!!

   double precision :: epsilon_c_pbeueg_bart1(N_states),epsilon_c_pbeueg_bart2(N_states)
   double precision :: beta_test_plus(N_states),beta_test2_plus(N_states)
   double precision :: beta_test_moins(N_states),beta_test2_moins(N_states)
   double precision :: d_beta_num_1,d_beta_num_2

   call give_Ec_pbeueg_test(mu,rhoa_plus,rhob,grad_rho_a,grad_rho_b,epsilon_c_pbeueg_bart1,epsilon_c_pbeueg_bart2,beta_test_plus,beta_test2_plus)
   call give_Ec_pbeueg_test(mu,rhoa_moins,rhob,grad_rho_a,grad_rho_b,epsilon_c_pbeueg_bart1,epsilon_c_pbeueg_bart2,beta_test_moins,beta_test2_moins)

   d_beta_num_1 = (beta_test_plus(1) - beta_test_moins(1)) /(2.d0 * drhoa)
   d_beta_num_2 = (beta_test2_plus(1) - beta_test2_moins(1)) /(2.d0 * drhoa)

  !print*,'epilon1,2 plus =',epsilon_c_pbeueg_bart1_plus(1),  epsilon_c_pbeueg_bart2_plus(1)
  !print*,'d_ec1,2         =',d_Ec_barth_num_1,d_Ec_barth_num_2
  !accu += dabs( d_Ec_barth_num_2 - d_Ec_barth_num_2)*weight


   call d_beta_pbeueg_rho(rhoa,rhob,grad_rho_a,grad_rho_b,d_beta_rhoa,d_beta_rhob)
 
   if (dabs(d_beta_rhoa(1) - d_beta_num_1)/dabs(d_beta_rhoa(1)) > 10.d0*drhoa ) then
    double precision :: xi, rs, wignerseitz_radius,f_pbeueg
    xi = (rhoa- rhob)/(rhoa+rhob)
    rs = wignerseitz_radius(rhot)

    print*,'****************************************'
    print*,'d_beta 1,2 num/ ana                  =',d_beta_num_1,d_beta_num_2,'/',d_beta_rhoa(1)
    print*,'rhot / rhot**2 / fpbeueg             =',rhot,'/',rhot**2.d0,'/',f_pbeueg(rs,xi)
!   print*,'Numerateur ana                       =',
    print*,'(rhot**2 * fpbeueg)**2 (denom ana)   =',(f_pbeueg(rs,xi)*rhot**2.d0)**2.d0
    print*,'Weight                               =',weight

    accu_dens_test_looser += rhot * weight
    accu_looser += dabs(d_beta_rhoa(1) - d_beta_num_1)*weight
   else
    accu_BG += dabs(d_beta_rhoa(1) - d_beta_num_1)*weight
    accu_dens_test_bg += rhot * weight
   endif
  
   accu += dabs(d_beta_rhoa(1) - d_beta_num_1)*weight
   accu_2 += dabs(d_beta_rhoa(1) - d_beta_num_2)*weight 

!!!!!!!!!!!!!!!!BETA!!!!!!!!!
 
  !delta_plus_b = delta_barth(rsplus_b,xiplus_b,mu)*(rhot+drhob)
  !delta_moins_b = delta_barth(rsmoins_b,ximoins_b,mu)*(rhot-drhob)
 
  !d_delta_barth_bis_b = (delta_plus_b - delta_moins_b) /(2.d0 * drhob)
 
  !accu_2 += dabs(d_delta_barth_bis - d_total_deltarho_rhob(rhoa,rhob,mu))*final_weight_at_r(l,k,j)
 
  enddo
  print*,'dn=n*                   =',10d0**(-n)
  print*,'************** Densite*********************' 
  print*,'accu_dens_test_looser   =', accu_dens_test_looser
  print*,'accu_dens_test_bg       =', accu_dens_test_bg 
  print*,'************** Accu derive*********************' 
  print*,'accu_looser             =', accu_looser 
  print*,'accu_BG                 =', accu_BG 
  print*,'accu tot                =',accu,accu_2
 enddo

end


subroutine test_total_derivative_bis_PBE_UEG
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),accu,accu_2
 double precision :: aos_array(ao_num),grad_rho_a(3,N_states),grad_rho_b(3,N_states),grad_aos_array(3,ao_num) 
 double precision :: rhot,rhos,rhoa,rhob,rhoa_plus,rhoa_moins,drhoa,drhob,pi,mu
 double precision :: beta_plus,beta_moins,threshold
 double precision :: d_beta_barth_bis,d_beta_rhoa,d_beta_rhob 
 double precision :: weight,accu_dens_test_looser,accu_looser,accu_BG,accu_dens_test_bg

 pi=dacos(-1.d0)
 mu=mu_erf
 threshold = 1d-10

!do n = 1, 1
 do n = 1, 16 
 accu = 0d0
 accu_2= 0d0
 accu_dens_test_looser  = 0.d0
 accu_looser  =0.d0 
 accu_BG  = 0.d0
 accu_dens_test_bg  = 0.d0

  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight=final_weight_at_r_vector(i)
  !call dm_dft_alpha_beta_and_all_aos_at_r(r,rhoa,rhob,aos_array)
   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rhoa,rhob, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
 
  !if(dabs(rhob).lt.threshold) cycle 
   rhot = rhoa + rhob
  !if(dabs(rhot).lt.threshold) cycle 
   rhos = rhoa - rhob
 
   drhoa = rhoa*10d0**(-n)
   drhob = rhob*10d0**(-n)
   rhoa_plus = rhoa + drhoa
   rhoa_moins = rhoa - drhoa
 
 
!!!!!!!!!!!!!!!!!!ALPHA!!!!!!!!!
  
   double precision :: epsilon_c_pbeueg_bart1_plus(N_states),epsilon_c_pbeueg_bart2_plus(N_states),beta_test(N_states),beta_test2(N_states)
   double precision :: epsilon_c_pbeueg_bart1_moins(N_states),epsilon_c_pbeueg_bart2_moins(N_states)
   double precision :: d_Ec_barth_num_1,d_Ec_barth_num_2
 
   call give_Ec_pbeueg_test(mu,rhoa_plus,rhob,grad_rho_a,grad_rho_b,epsilon_c_pbeueg_bart1_plus,epsilon_c_pbeueg_bart2_plus,beta_test,beta_test2)
   call give_Ec_pbeueg_test(mu,rhoa_moins,rhob,grad_rho_a,grad_rho_b,epsilon_c_pbeueg_bart1_moins,epsilon_c_pbeueg_bart2_moins,beta_test,beta_test2)
 
   d_Ec_barth_num_1 = (epsilon_c_pbeueg_bart1_plus(1) - epsilon_c_pbeueg_bart1_moins(1)) /(2.d0 * drhoa)
   d_Ec_barth_num_2 = (epsilon_c_pbeueg_bart2_plus(1) - epsilon_c_pbeueg_bart2_moins(1)) /(2.d0 * drhoa)
 
  double precision :: d_ec_pbeueg_rhoa(N_states),d_ec_pbeueg_rhob(N_states)
  call give_d_Ec_pbeueg_rho(rhoa,rhob,grad_rho_a,grad_rho_b,mu,d_ec_pbeueg_rhoa,d_ec_pbeueg_rhob)
 
   if (dabs(d_ec_pbeueg_rhoa(1) - d_Ec_barth_num_1) > 10.d0*drhoa ) then
    print*,'****************************************'
    print*,'d_ec_pbe 1,2 num/ ana                =',d_Ec_barth_num_1,d_Ec_barth_num_2,'/',d_ec_pbeueg_rhoa(1)
    print*,'rhot / rhot**2 / fpbeueg             =',beta_test,beta_test2
    print*,'(1+beta*mu**3)**2 (denom ana)        =',(1.d0+beta_test*mu**3.d0)**2.d0,'/',(1.d0+beta_test2*mu**3.d0)**2.d0
    print*,'Weight                               =',weight

    accu_dens_test_looser += rhot * weight
    accu_looser += dabs(d_ec_pbeueg_rhoa(1) - d_Ec_barth_num_1)*weight
   else
    accu_BG += dabs(d_ec_pbeueg_rhoa(1) - d_Ec_barth_num_1)*weight 
    accu_dens_test_bg += rhot * weight
   endif


   accu += dabs(d_ec_pbeueg_rhoa(1) - d_Ec_barth_num_1)*weight
   accu_2 += dabs(d_ec_pbeueg_rhoa(1) - d_Ec_barth_num_2)*weight

!! accu += dabs(d_beta_barth_bis)*final_weight_at_r(l,k,j)
 
!! print*,'beta_plus,beta_moins,d_beta_barth_bis,d_beta_rhoa',beta_plus,beta_moins,d_beta_barth_bis,d_beta_rhoa
!!!!!!!!!!!!!!!!BETA!!!!!!!!!
 
  !delta_plus_b = delta_barth(rsplus_b,xiplus_b,mu)*(rhot+drhob)
  !delta_moins_b = delta_barth(rsmoins_b,ximoins_b,mu)*(rhot-drhob)
 
  !d_delta_barth_bis_b = (delta_plus_b - delta_moins_b) /(2.d0 * drhob)
 
  !accu_2 += dabs(d_delta_barth_bis - d_total_deltarho_rhob(rhoa,rhob,mu))*final_weight_at_r(l,k,j)
 
  enddo
  !print*,'dxi=xi*',10d0**(-n),accu_2,accu_3,accu_4,accu_5,accu_6
  print*,'************** YOLO*********************' 
  print*,'dn=n*                   =',10d0**(-n)
  print*,'************** Densite*********************' 
  print*,'accu_dens_test_looser   =', accu_dens_test_looser
  print*,'accu_dens_test_bg       =', accu_dens_test_bg 
  print*,'************** Accu derive*********************' 
  print*,'accu_looser             =', accu_looser 
  print*,'accu_BG                 =', accu_BG 
  print*,'accu tot                =',accu,accu_2
 enddo

end


subroutine test_pot_pbe_ueg
 implicit none
 integer :: i,j,istate
 double precision :: accu_alpha,accu_beta,accu_contri_grad_alpha,accu_contri_grad_beta 
 
 accu_alpha = 0d0
 accu_beta = 0d0

 do istate = 1, N_states

  do i = 1, ao_num
   do j = 1, ao_num
    accu_alpha += dabs(potential_c_alpha_ao_sr_pbe_ueg(j,i,istate)-potential_c_alpha_ao_sr_pbe(j,i,istate))
    accu_beta  += dabs(potential_c_beta_ao_sr_pbe_ueg(j,i,istate)-potential_c_beta_ao_sr_pbe(j,i,istate))
   enddo
  enddo
 
  print*,'*************************'
  print*,'Limite en 0'
  print*,'mu_erf_dft  =', mu_erf_dft 
  print*,'accu_alpha / accu_beta  =', accu_alpha,'/',accu_beta


  accu_alpha = 0.d0
  accu_beta = 0.d0
  accu_contri_grad_alpha =0.d0
  accu_contri_grad_beta  =0.d0

  do i = 1, ao_num
   do j = 1, ao_num
    accu_contri_grad_alpha += dabs(pot_sr_grad_c_alpha_ao_pbe_ueg(j,i,istate)) 
    accu_contri_grad_beta  += dabs(pot_sr_grad_c_beta_ao_pbe_ueg(j,i,istate)) 
    accu_alpha += dabs(potential_c_alpha_ao_sr_pbe_ueg(j,i,istate)-potential_c_alpha_ao_sr_pbe_ueg_grandmu(j,i,istate))
    accu_beta  += dabs(potential_c_beta_ao_sr_pbe_ueg(j,i,istate)-potential_c_beta_ao_sr_pbe_ueg_grandmu(j,i,istate))
   !print*,'***********************************************'
   !print*,'sr_pbe_ueg / sr_pbe    =', potential_c_alpha_ao_sr_pbe_ueg(j,i,istate),'/',potential_c_alpha_ao_sr_pbe(j,i,istate)
   !print*,'sr_pbe_ueg scal / grad =', pot_sr_scal_c_alpha_ao_pbe_ueg(j,i,istate),'/',pot_sr_grad_c_alpha_ao_pbe_ueg(j,i,istate)
   enddo
  enddo

  print*,'*************************'
  print*,'Limite en linfini'
  print*,'mu_erf_dft  =', mu_erf_dft 
  print*,'accu_grad_contriub_alpha/beta  =', accu_contri_grad_alpha,'/',accu_contri_grad_beta
  print*,'accu_alpha / accu_beta  =', accu_alpha,'/',accu_beta

 enddo

end
