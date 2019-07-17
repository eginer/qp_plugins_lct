program pouet
 read_wf = .True.
 touch read_wf
 !call test_delta      
 !call test_delta_derivative
 !call test_delta_derivative_xi
 !call test_delta_total_derivative 
  call test_delta_total_derivative_bis 
 !call test_provider
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



subroutine test_delta  
 implicit none
 integer :: i,j,k,l,m,n
 print*,'*****************$******************'
 print*,'ECMD LDA standart =',Energy_c_md_LDA
 print*,'ECMD LDA mu_of_r  =',Energy_c_md_LDA_mu_of_r 
 print*,'ECMD LDA barth    =',Energy_c_md_LDA_barth
 print*,'*****************$******************'
end


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
