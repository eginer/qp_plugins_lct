
 BEGIN_PROVIDER [double precision, one_e_act_dm_alpha_mo_for_dft, (n_act_orb,n_act_orb,N_states)]
 implicit none
 integer :: i,j,ii,jj,istate
 do istate = 1, N_states
  do ii = 1, n_act_orb
   i = list_act(ii)
   do jj = 1, n_act_orb
    j = list_act(jj)
    one_e_act_dm_alpha_mo_for_dft(jj,ii,istate) = one_e_dm_mo_alpha_for_dft(j,i,istate)
   enddo
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, one_e_act_density_alpha,(n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
 ! one_e_act_density = pure act part of the STATE AVERAGED alpha density
 END_DOC
 one_e_act_density_alpha = 0.d0
 integer :: ipoint,i,j,istate
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   do i = 1, n_act_orb
    do j = 1, n_act_orb
     one_e_act_density_alpha(ipoint,istate) += one_e_act_dm_alpha_mo_for_dft(j,i,istate) * act_mos_in_r_array(j,ipoint) * act_mos_in_r_array(i,ipoint)
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER 




 BEGIN_PROVIDER [double precision, one_e_act_dm_beta_mo_for_dft, (n_act_orb,n_act_orb,N_states)]
 implicit none
 integer :: i,j,ii,jj,istate
 do istate = 1, N_states
  do ii = 1, n_act_orb
   i = list_act(ii)
   do jj = 1, n_act_orb
    j = list_act(jj)
    one_e_act_dm_beta_mo_for_dft(jj,ii,istate) = one_e_dm_mo_beta_for_dft(j,i,istate)
   enddo
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, one_e_act_density_beta,(n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
 ! one_e_act_density = pure act part of the STATE AVERAGED beta density
 END_DOC
 one_e_act_density_beta = 0.d0
 integer :: ipoint,i,j,istate
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   do i = 1, n_act_orb
    do j = 1, n_act_orb
     one_e_act_density_beta(ipoint,istate) += one_e_act_dm_beta_mo_for_dft(j,i,istate) * act_mos_in_r_array(j,ipoint) * act_mos_in_r_array(i,ipoint)
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER 



