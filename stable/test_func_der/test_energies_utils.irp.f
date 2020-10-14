subroutine delta_gamma_i_j_for_energy_test (delta_rho_11,delta_gamma_i_j)
 implicit none
 BEGIN_DOC
!
 END_DOC

 double precision, intent(in)  :: delta_rho_11
 double precision, intent(out) :: delta_gamma_i_j(mo_num,mo_num)
 
 integer :: istate,i,j
 double precision :: r(3), mo_i(mo_num)


!  do i = 1, n_points_final_grid
!   r(1) = final_grid_points(1,i)   
!   r(2) = final_grid_points(2,i) 
!   r(3) = final_grid_points(3,i)
!   call give_all_mos_at_r(r, mo_i)
   do i=1, mo_num
    do j=1, mo_num
     delta_gamma_i_j(i,j) = 0.d0
    enddo
   enddo
   delta_gamma_i_j(1,1) = delta_rho_11
end

subroutine delta_gamma_i_j_k_l_for_energy_test (delta_rho_11,delta_gamma_i_j_k_l)
 implicit none
 BEGIN_DOC
!
 END_DOC

 double precision, intent(in)  :: delta_rho_11
 double precision, intent(out) :: delta_gamma_i_j_k_l(mo_num,mo_num,mo_num,mo_num)
 
 integer :: istate,j,k,l,m
 double precision :: r(3), mo_i(mo_num)


!  do i = 1, n_points_final_grid
!   r(1) = final_grid_points(1,i)   
!   r(2) = final_grid_points(2,i) 
!   r(3) = final_grid_points(3,i)
!   call give_all_mos_at_r(r, mo_i)
   do j=1, mo_num
    do k=1, mo_num
     do l=1, mo_num
      do m=1, mo_num
       delta_gamma_i_j_k_l(j,k,l,m) = 0.d0
      enddo
     enddo
    enddo
   enddo
   delta_gamma_i_j_k_l(1,1,1,1) = delta_rho_11
end

subroutine delta_density_for_energy_test (delta_rho_a, delta_rho_b, delta_rho_11)
 implicit none
 BEGIN_DOC
!
 END_DOC

 double precision, intent(in)  :: delta_rho_11
 double precision, intent(out) :: delta_rho_a(n_points_final_grid), delta_rho_b(n_points_final_grid)
 
 integer :: istate,i,j
 double precision :: r(3), mo_i(mo_num)


  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   call give_all_mos_at_r(r, mo_i)

   delta_rho_a(i) = delta_rho_11  * mo_i(1) * mo_i(1)   !delta_rho_ij = delta_rho_11 = kro(i1)kro(j1)*epsilon
   delta_rho_b(i) = delta_rho_11  * mo_i(1) * mo_i(1)  

  enddo
end

subroutine delta_grad_density_for_energy_test (delta_grad_rho_a, delta_grad_rho_b, delta_grad_rho_11)
 implicit none

 double precision, intent(in)  :: delta_grad_rho_11
 double precision, intent(out) :: delta_grad_rho_a(3,n_points_final_grid), delta_grad_rho_b(3,n_points_final_grid)
 integer :: i,m
 double precision :: r(3), mo_i(mo_num), mo_grad_i(3,mo_num)


  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i)   
   r(3) = final_grid_points(3,i)   
 
   call give_all_mos_and_grad_at_r(r,mo_i,mo_grad_i)

   do m=1,3
     delta_grad_rho_a(m,i) = 2.d0 * delta_grad_rho_11 * mo_i(1) * mo_grad_i(m,1) 
     delta_grad_rho_b(m,i) = 2.d0 * delta_grad_rho_11 * mo_i(1) * mo_grad_i(m,1)
   enddo

  enddo

end

subroutine delta_n2_for_energy_test (delta_n2, delta_n2_11)
 implicit none
 BEGIN_DOC
!
 END_DOC

 double precision, intent(in)  :: delta_n2_11
 double precision, intent(out) :: delta_n2(n_points_final_grid)
 
 integer :: istate,i,j
 double precision :: r(3), mo_i(mo_num)


  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   call give_all_mos_at_r(r, mo_i)

   delta_n2(i) = delta_n2_11  * mo_i(1) * mo_i(1) * mo_i(1) * mo_i(1)

  enddo
end

subroutine write_energies_test_on_file(delta_rho_11, ex_lda_test, ex_lda_test_plus_delta,int_vx_lda_test, ex_pbe_test, ex_pbe_test_plus_delta,int_vx_pbe_test, ec_pbe_test, ec_pbe_test_plus_delta,int_vc_pbe_test,ex_pbeUEG_test, ex_pbeUEG_test_plus_delta,int_vx_pbeUEG_test, ec_pbeUEG_test, ec_pbeUEG_test_plus_delta,int_vc_pbeUEG_test,ec_pben2_test, ec_pben2_test_plus_delta,int_vc_pben2_one_e_test,int_vc_pben2_two_e_test)

  implicit none
  BEGIN_DOC
  ! Write small variations of energies for the test of the energies potentials
  END_DOC

  double precision, intent(in) :: delta_rho_11
  double precision, intent(in) :: ex_lda_test, ex_lda_test_plus_delta,int_vx_lda_test
  double precision, intent(in) :: ex_pbe_test, ex_pbe_test_plus_delta,int_vx_pbe_test
  double precision, intent(in) :: ec_pbe_test, ec_pbe_test_plus_delta,int_vc_pbe_test
  double precision, intent(in) :: ex_pbeUEG_test, ex_pbeUEG_test_plus_delta,int_vx_pbeUEG_test
  double precision, intent(in) :: ec_pbeUEG_test, ec_pbeUEG_test_plus_delta,int_vc_pbeUEG_test
  double precision, intent(in) :: ec_pben2_test, ec_pben2_test_plus_delta,int_vc_pben2_one_e_test,int_vc_pben2_two_e_test
  write(34,*) 'gamma_i_j = epsilon*delta_kronecker(i,1)*delta_kronecker(1,j) with epsilon =', delta_rho_11
   write(34,*) '----------------------------------------------LDA--------------------------------------------' 
   write(34,*) 'ex_lda_test[n]                                                          =', ex_lda_test
  ! write(34,*) 'ex_lda(provider)[n]                                                     =', energy_x_lda(1)
   write(34,*) 'ex_lda_test[n + delta_n]                                                =', ex_lda_test_plus_delta
   write(34,*) '------------------------------------------------'
  ! write(34,*) 'Theoretical value : int (dr delta_n * [delta (E) / delta(n(r))] )       =', integral_potential_lda
  ! write(34,*) 'We test the equality between the two following results :'
   write(34,*) '1st_method = ex_lda [n + delta_n] - ex_lda[n]                           =', ex_lda_test_plus_delta - ex_lda_test
   write(34,*) '2nd method = sum(i,j) delta_gamma_i_j * int (dr phi_i(r) v(r) phi_j(r)) =', int_vx_lda_test
   write(34,*) 'Relative error : (1st_method - 2nd_method)/1st_method                   =', (int_vx_lda_test - ex_lda_test_plus_delta + ex_lda_test)/(ex_lda_test_plus_delta - ex_lda_test)
   write(34,*) '------------------------------------------------'
   write(34,*) '------------------------------------------------'
   write(34,*) '----------------------------------------------PBE--------------------------------------------'
   write(34,*) '1st_method = ex_pbe [n + delta_n] - ex_pbe[n]                           =', ex_pbe_test_plus_delta - ex_pbe_test
   write(34,*) '2nd method = sum(i,j) delta_gamma_i_j * int (dr phi_i(r) v(r) phi_j(r)) =', int_vx_pbe_test
   write(34,*) 'Relative error : (1st_method - 2nd_method)/1st_method                   =', (int_vx_pbe_test - ex_pbe_test_plus_delta + ex_pbe_test)/(ex_pbe_test_plus_delta - ex_pbe_test)

   write(34,*) '------------------------------------------------'
   write(34,*) '1st_method = ec_pbe [n + delta_n] - ec_pbe[n]                           =', ec_pbe_test_plus_delta - ec_pbe_test
   write(34,*) '2nd method = sum(i,j) delta_gamma_i_j * int (dr phi_i(r) v(r) phi_j(r)) =', int_vc_pbe_test
   write(34,*) 'Relative error : (1st_method - 2nd_method)/1st_method                   =', (int_vc_pbe_test - ec_pbe_test_plus_delta + ec_pbe_test)/(ec_pbe_test_plus_delta - ec_pbe_test)


   write(34,*) '----------------------------------------------PBEUEG--------------------------------------------'
   write(34,*) '1st_method = ex_pbeUEG [n + delta_n] - ex_pbeUEG[n]                           =', ex_pbeUEG_test_plus_delta - ex_pbeUEG_test
   write(34,*) '2nd method = sum(i,j) delta_gamma_i_j * int (dr phi_i(r) v(r) phi_j(r)) =', int_vx_pbeUEG_test
   write(34,*) 'Relative error : (1st_method - 2nd_method)/1st_method                   =', (int_vx_pbeUEG_test - ex_pbeUEG_test_plus_delta + ex_pbeUEG_test)/(ex_pbeUEG_test_plus_delta - ex_pbeUEG_test)
   
   write(34,*) '------------------------------------------------'                                        
   write(34,*) '1st_method = ec_pbeUEG [n + delta_n] - ec_pbeUEG[n]                           =', ec_pbeUEG_test_plus_delta - ec_pbeUEG_test
   write(34,*) '2nd method = sum(i,j) delta_gamma_i_j * int (dr phi_i(r) v(r) phi_j(r)) =', int_vc_pbeUEG_test
   write(34,*) 'Relative error : (1st_method - 2nd_method)/1st_method                   =', (int_vc_pbeUEG_test - ec_pbeUEG_test_plus_delta + ec_pbeUEG_test)/(ec_pbeUEG_test_plus_delta - ec_pbeUEG_test)

   write(34,*) '----------------------------------------------PBEn2--------------------------------------------'
   write(34,*) '1st_method = ec_pben2 [n + delta_n] - ec_pben2[n]                           =', ec_pben2_test_plus_delta - ec_pben2_test
   write(34,*) '2nd method =                                                                =', int_vc_pben2_one_e_test + int_vc_pben2_two_e_test
   write(34,*) 'Relative error : (1st_method - 2nd_method)/1st_method                   =', (int_vc_pben2_one_e_test + int_vc_pben2_two_e_test - ec_pben2_test_plus_delta + ec_pben2_test)/(ec_pben2_test_plus_delta - ec_pben2_test)

end
