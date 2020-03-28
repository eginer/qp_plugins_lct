 BEGIN_PROVIDER [double precision, Energy_c_md_LDA, (N_states)]
 implicit none
 BEGIN_DOC
 ! Corelation energy for the multi determinent short range LDA. PRB 73 155111 2006
 !
 ! Depends on the total density and spin polarization
 END_DOC
 integer :: ipoint,istate 
 double precision, allocatable :: aos_array(:), r(:), rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: r2(3),dr2(3), local_potential,r12,dx2,mu,coulomb,two_body_dm,weight
 double precision :: threshold
 double precision :: cpu0,cpu1
 dospin = .True. 
 threshold = 1.d-07
 Energy_c_md_LDA = 0.d0
 allocate(aos_array(ao_num),r(3), rho_a(N_states), rho_b(N_states), ec(N_states))
 call cpu_time(cpu0)
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
   weight = final_weight_at_r_vector(ipoint)
   do istate = 1, N_states
    call ESRC_MD_LDAERF (mu_erf_dft,rho_a(istate),rho_b(istate),dospin,ec(istate))
    Energy_c_md_LDA(istate) += weight * ec(istate)
   enddo
  enddo
 deallocate(aos_array,r,rho_a,rho_b, ec)
 call cpu_time(cpu1)
 print*,'Time for the Energy_c_md_LDA integration :',cpu1-cpu0
END_PROVIDER

