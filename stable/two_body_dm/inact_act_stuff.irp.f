 BEGIN_PROVIDER[double precision, inact_act_mos_in_r_array, (mo_num,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, inact_act_mos_in_r_array_transp,(n_points_final_grid,mo_num)]
 implicit none
 integer :: i,j,k
 do i = 1, n_inact_act
  j = list_inact_act(i) 
  do k = 1, n_points_final_grid
   inact_act_mos_in_r_array_transp(k,i) = mos_in_r_array_transp(k,j)
  enddo
 enddo

 do k = 1, n_points_final_grid
  do i = 1, n_inact_act
   inact_act_mos_in_r_array(i,k) = inact_act_mos_in_r_array_transp(k,i)
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER[double precision, core_inact_act_mos_in_r_array, (n_core_inact_act_orb,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, core_inact_act_mos_in_r_array_transp,(n_points_final_grid,n_core_inact_act_orb)]
 implicit none
 integer :: i,j,k
 do i = 1, n_core_inact_act_orb
  j = list_core_inact_act(i) 
  do k = 1, n_points_final_grid
   core_inact_act_mos_in_r_array_transp(k,i) = mos_in_r_array_transp(k,j)
  enddo
 enddo

 do k = 1, n_points_final_grid
  do i = 1, n_core_inact_act_orb
   core_inact_act_mos_in_r_array(i,k) = core_inact_act_mos_in_r_array_transp(k,i)
  enddo
 enddo

END_PROVIDER 


 BEGIN_PROVIDER [double precision, core_inact_act_two_bod_alpha_beta_mo, (n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,N_states)]
 implicit none
 BEGIN_DOC
 !  core_inact_act_two_bod_alpha_beta(i,j,k,l) = <Psi| a^{dagger}_{j,alpha} a^{dagger}_{l,beta} a_{k,beta} a_{i,alpha} | Psi>
 !                     1 1 2 2  = chemist notations 
 !  note that no 1/2 factor is introduced in order to take into acccount for the spin symmetry
 !  
 END_DOC
 integer :: dim1,dim2,dim3,dim4
 double precision :: cpu_0,cpu_1
 dim1 = n_core_inact_act_orb
 dim2 = n_core_inact_act_orb
 dim3 = n_core_inact_act_orb
 dim4 = n_core_inact_act_orb
 core_inact_act_two_bod_alpha_beta_mo = 0.d0
 print*,'providing core_inact_act_two_bod_alpha_beta_mo ...'
 call wall_time(cpu_0)
 call two_body_dm_nstates_openmp(core_inact_act_two_bod_alpha_beta_mo,dim1,dim2,dim3,dim4,psi_coef,size(psi_coef,2),size(psi_coef,1))
 call wall_time(cpu_1)
 print*,'core_inact_act_two_bod_alpha_beta_mo provided in',dabs(cpu_1-cpu_0)

 integer :: ii,jj,i,j,k,l
 if(no_core_density .EQ. "no_core_dm")then
  print*,'USING THE VALENCE ONLY TWO BODY DENSITY'

  do ii = 1, n_core_orb ! 1 
   i = list_core(ii) 
   do j = 1, n_core_inact_act_orb ! 2 
    do k = 1, n_core_inact_act_orb  ! 1 
     do l = 1, n_core_inact_act_orb ! 2 
      !                     2 2 1 1
      core_inact_act_two_bod_alpha_beta_mo(l,j,k,i,:) = 0.d0
      core_inact_act_two_bod_alpha_beta_mo(j,l,k,i,:) = 0.d0
      core_inact_act_two_bod_alpha_beta_mo(l,j,i,k,:) = 0.d0
      core_inact_act_two_bod_alpha_beta_mo(j,l,i,k,:) = 0.d0

      core_inact_act_two_bod_alpha_beta_mo(k,i,l,j,:) = 0.d0
      core_inact_act_two_bod_alpha_beta_mo(k,i,j,l,:) = 0.d0
      core_inact_act_two_bod_alpha_beta_mo(i,k,l,j,:) = 0.d0
      core_inact_act_two_bod_alpha_beta_mo(i,k,j,l,:) = 0.d0
     enddo
    enddo
   enddo
  enddo


 endif

 END_PROVIDER 



 BEGIN_PROVIDER [double precision, core_inact_act_two_bod_alpha_beta_mo_physicist, (n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,N_states)]
 implicit none
 BEGIN_DOC
 !  core_inact_act_two_bod_alpha_beta_mo_physicist,(i,j,k,l) = <Psi| a^{dagger}_{k,alpha} a^{dagger}_{l,beta} a_{j,beta} a_{i,alpha} | Psi>
 !                                   1 2 1 2  = physicist notations 
 !  note that no 1/2 factor is introduced in order to take into acccount for the spin symmetry
 !  
 END_DOC
 integer :: i,j,k,l,istate
 double precision :: cpu_0,cpu_1
 core_inact_act_two_bod_alpha_beta_mo_physicist = 0.d0
 print*,'providing core_inact_act_two_bod_alpha_beta_mo_physicist ...'
 call wall_time(cpu_0)
 do istate = 1, N_states 
  do i = 1, n_core_inact_act_orb
   do j = 1, n_core_inact_act_orb
    do k = 1, n_core_inact_act_orb
     do l = 1, n_core_inact_act_orb
      !                               1 2 1 2                                 1 1 2 2 
      core_inact_act_two_bod_alpha_beta_mo_physicist(l,k,i,j,istate) = core_inact_act_two_bod_alpha_beta_mo(i,l,j,k,istate)
     enddo
    enddo
   enddo
  enddo
 enddo
 call wall_time(cpu_1)
 print*,'core_inact_act_two_bod_alpha_beta_mo_physicist provided in',dabs(cpu_1-cpu_0)

 END_PROVIDER 


 double precision function core_inact_act_on_top_of_r_from_provider(ipoint,istate)
 implicit none
 BEGIN_DOC
 ! on top pair density evaluated at a given point of the grid 
 END_DOC
 integer, intent(in) :: ipoint,istate
 integer :: i,j,k,l
 core_inact_act_on_top_of_r_from_provider = 0.d0
 do l = 1, n_core_inact_act_orb
  do k = 1, n_core_inact_act_orb
    do j = 1, n_core_inact_act_orb
     do i = 1, n_core_inact_act_orb
     !                                                                                          1 2 1 2 
     core_inact_act_on_top_of_r_from_provider += core_inact_act_two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) * core_inact_act_mos_in_r_array(j,ipoint) * core_inact_act_mos_in_r_array(i,ipoint) * core_inact_act_mos_in_r_array(l,ipoint) * core_inact_act_mos_in_r_array(k,ipoint)
    enddo
   enddo
  enddo
 enddo
 end

 BEGIN_PROVIDER [double precision, core_inact_act_on_top_of_r,(n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, wall_time_core_inact_act_on_top_of_r ]
 implicit none
 BEGIN_DOC
 ! on top pair density at each grid point computed using the full two-body density matrix 
 END_DOC
 integer :: i_point,i_state
 double precision :: wall_0,wall_1
 double precision :: core_inact_act_on_top_of_r_from_provider

 print*,'providing the core_inact_act_on_top_of_r'
 i_point = 1
 provide core_inact_act_two_bod_alpha_beta_mo_physicist 
 i_state = 1
 core_inact_act_on_top_of_r(i_point,i_state) = core_inact_act_on_top_of_r_from_provider(i_point,i_state)
 call wall_time(wall_0)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,i_state) & 
 !$OMP SHARED(core_inact_act_on_top_of_r,n_points_final_grid,N_states)
 do i_point = 1, n_points_final_grid
  do i_state = 1, N_states
   core_inact_act_on_top_of_r(i_point,i_state) = core_inact_act_on_top_of_r_from_provider(i_point,i_state)
  enddo
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall_1)
 print*,'provided the core_inact_act_on_top_of_r'
 print*,'Time to provide :',wall_1 - wall_0
 wall_time_core_inact_act_on_top_of_r = wall_1 - wall_0

 END_PROVIDER 

 
