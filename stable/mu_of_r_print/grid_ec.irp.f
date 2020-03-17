
 BEGIN_PROVIDER [double precision, E_c_md_mur_grid_n_and_on_top_PBE, (N_states,n_points_print_mur)]
 BEGIN_DOC
  ! E_c_md_mur_grid_n_and_on_top_PBE = PBE-on_top multi determinant functional with exact on top extrapolated for large mu using a mu(r) interaction and spin polarization computed only with on-top and total density
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: mu,pi
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_on_top_PBE(:),two_dm(:)
 pi = 4.d0 * datan(1.d0)

 allocate(eps_c_md_on_top_PBE(N_states),two_dm(N_states))
 E_c_md_mur_grid_n_and_on_top_PBE = 0.d0
  
 print*,'Providing E_c_md_mur_grid_n_and_on_top_PBE ...'
 call wall_time(wall0)
 do i = 1, n_points_print_mur
  r(1) = grid_points_mur(1,i)
  r(2) = grid_points_mur(2,i)
  r(3) = grid_points_mur(3,i)
  two_dm(:) = total_cas_on_top_density_grid_mur(i,:)
  mu = cas_full_mu_of_r_grid_mur_psi_coal_vector(i)

  call give_epsilon_c_md_n_and_on_top_PBE_mu_corrected_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  do istate = 1, N_states
   E_c_md_mur_grid_n_and_on_top_PBE(istate,i) = eps_c_md_on_top_PBE(istate) 
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the E_c_md_mur_grid_n_and_on_top_PBE :',wall1-wall0

 END_PROVIDER

BEGIN_PROVIDER [double precision, E_c_md_mur_grid_PBE, (N_states,n_points_print_mur)]
 BEGIN_DOC
  ! E_c_md_mur_grid_PBE_vector           = PBE multi determinant functional with UEG on top for large mu using a mu(r) interaction (JT)
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_PBE(:)

 allocate(eps_c_md_PBE(N_states))
 E_c_md_mur_grid_PBE = 0.d0
  
 print*,'Providing E_c_md_mur_grid_PBE ...'
 call wall_time(wall0)
 do i = 1, n_points_print_mur
  r(1) = grid_points_mur(1,i)
  r(2) = grid_points_mur(2,i)
  r(3) = grid_points_mur(3,i)
  mu = cas_full_mu_of_r_grid_mur_psi_coal_vector(i)

  call give_epsilon_c_md_PBE_mu(mu,r,eps_c_md_PBE)
  do istate = 1, N_states
   E_c_md_mur_grid_PBE(istate,i) += eps_c_md_PBE(istate) 
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the E_c_md_mur_grid_PBE:',wall1-wall0

 END_PROVIDER

