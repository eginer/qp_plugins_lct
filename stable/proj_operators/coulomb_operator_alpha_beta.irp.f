subroutine f_PSI_ab_routine(r1,r2,coulomb,two_body_dm)
 BEGIN_DOC
 ! 2 * f_{\Psi}(r1,r2) = Equation A16 of J. Chem. Phys. 149, 194301 (2018) for alpha/beta electrons
 ! int dr1 dr2 f_{\Psi}(r1,r2) = <\Psi | Wee_{alpha/beta} | \Psi>
 END_DOC
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: coulomb, two_body_dm
 integer :: i,j,k,l,m,n  
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: a1,a2,a3,a4
 double precision :: c1,c2,c3,c4
 double precision :: threshold
 double precision :: integral,get_mo_bielec_integral
 threshold = 0.d0

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 two_body_dm = 0.d0
 coulomb = 0.d0
 double precision :: tmp1
 do i = 1, mo_tot_num ! 1
  do j = 1, mo_tot_num ! 2 
   do m = 1, mo_tot_num ! 1
    do n = 1, mo_tot_num ! 2
     !                                               2 1 2 1
     tmp1 = mos_array_r2(n) * mos_array_r1(m) * mos_array_r2(j) * mos_array_r1(i)
     coulomb     += two_bod_alpha_beta_mo_contracted(n,m,j,i,1) * tmp1 
     two_body_dm += two_bod_alpha_beta_mo_physicist (n,m,j,i,1) * tmp1 
    enddo
   enddo
  enddo
 enddo

end

double precision function integral_of_f_PSI_ab_over_2(r)
 implicit none
! integral over 2 of f_PSI_ab(1,2) 
 double precision, intent(in) :: r(3)
 double precision :: mos_array_r(mo_tot_num)
 integer :: i,j,m,n
 call give_all_mos_at_r(r,mos_array_r)
 integral_of_f_PSI_ab_over_2 = 0.d0
 do i = 1, mo_tot_num ! 1
  do j = 1, mo_tot_num ! 2 
   do m = 1, mo_tot_num ! 1
    do n = j,j           ! 2 :: kronecker(n,j) because of orthogonality of MOs
     !                                                                                 2 1 2 1
     integral_of_f_PSI_ab_over_2 += two_bod_alpha_beta_mo_contracted(n,m,j,i,1) * mos_array_r(m) * mos_array_r(i) 
    enddo
   enddo
  enddo
 enddo
end

subroutine f_HF_ab(r1,r2,integral_psi,two_bod)
 implicit none
 BEGIN_DOC
! f_HF_ab(X1,X2) = function f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for alpha beta spins and an HF wave function
! < HF | wee_{\alpha\alpha} | HF > = 0.5 * \int (X1,X2) f_HF_aa(X1,X2)
 END_DOC
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: integral_psi,two_bod
 integer :: i,j,m,n
 integer :: ii,jj,mm,nn
 double precision :: mo_bielec_integral
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: get_mo_bielec_integral
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 integral_psi = 0.d0
 two_bod = 0.d0
 do m = 1, elec_alpha_num
  do n = 1, elec_beta_num
   two_bod += mos_array_r1(n) * mos_array_r1(n) * mos_array_r2(m) * mos_array_r2(m) 
   do i = 1, mo_tot_num
    do j = 1, mo_tot_num
     integral_psi +=  integrals_for_hf_potential(j,i,n,m) * mos_array_r1(i) * mos_array_r2(j) * mos_array_r2(n) * mos_array_r1(m) 
    enddo
   enddo
  enddo
 enddo
end


BEGIN_PROVIDER [double precision, integrals_for_hf_potential, (mo_tot_num,mo_tot_num,elec_alpha_num,elec_alpha_num)]
 implicit none
 integer :: i,j,m,n
 double precision :: get_mo_bielec_integral
 do m = 1, elec_alpha_num ! electron 1 
  do n = 1, elec_alpha_num ! electron 2 
   do i = 1, mo_tot_num   ! electron 1 
    do j = 1, mo_tot_num  ! electron 2 
     integrals_for_hf_potential(j,i,n,m) = get_mo_bielec_integral(m,n,i,j,mo_integrals_map) 
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, integrals_for_hf_potential_integrated_on_beta, (mo_tot_num,elec_alpha_num)]
 implicit none
 integer :: i,j,k
 double precision :: get_mo_bielec_integral
 integrals_for_hf_potential_integrated_on_beta = 0.d0
 do k = 1, elec_alpha_num ! electron 1 alpha 
  do i = 1, mo_tot_num   ! electron 1 alpha
   do j = 1, elec_beta_num ! electron 2 beta 
    integrals_for_hf_potential_integrated_on_beta(i,k) += get_mo_bielec_integral(i,j,k,j,mo_integrals_map) 
   enddo
  enddo
 enddo
END_PROVIDER 


subroutine integral_of_f_12_hf_over_beta(r1,integral_f)
 implicit none
 BEGIN_DOC
! computes the following ANALYTICAL integral
! integral_(r2,beta) f_HF(r1,r2) with f_HF being Equation A16 of J. Chem. Phys. 149, 194301 (2018) with \Psi = HF
 END_DOC
 double precision, intent(in) :: r1(3)
 double precision, intent(out):: integral_f
 integer :: i,k
 double precision :: mos_array_r1(mo_tot_num)
 call give_all_mos_at_r(r1,mos_array_r1)

 integral_f = 0.d0
 do k = 1, elec_alpha_num
  do i = 1, mo_tot_num
   integral_f += integrals_for_hf_potential_integrated_on_beta(i,k) * mos_array_r1(i) * mos_array_r1(k)
  enddo
 enddo

end

 BEGIN_PROVIDER [double precision, integral_f_hf ]
 implicit none
 BEGIN_DOC 
! int dr1 dr2 f_HF_ab(r1,r2) with f_HF being Equation A16 of J. Chem. Phys. 149, 194301 (2018) with \Psi = HF
! should be the average value of the alpha/beta bielectronic repulsion over the HF wave function
 END_DOC
 double precision :: r(3),integral_f,weight
 integer :: i_point
 integral_f_hf= 0.D0
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call integral_of_f_12_hf_over_beta(r,integral_f)
  weight = final_weight_functions_at_final_grid_points(i_point)
  integral_f_hf += weight * integral_f
 enddo
 END_PROVIDER 

