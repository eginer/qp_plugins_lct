
 BEGIN_PROVIDER [double precision, ecmd_pbe_ueg_prov, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a large mu behaviour in function of the UEG on top pair density coupled to the PBE correlation energy at mu=0
  ! Ec_md_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_PBE  with beta chosen to recover the UEG large mu behaviour (JT)
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_PBE(:)
 allocate(eps_c_md_PBE(N_states))
 mu = mu_erf_dft
 ecmd_pbe_ueg_prov = 0.d0
  
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight = final_weight_at_r_vector(i)
  call ecmd_pbe_ueg_at_r(mu,r,eps_c_md_PBE)
  do istate = 1, N_states
   ecmd_pbe_ueg_prov(istate) += eps_c_md_PBE(istate) * weight
  enddo
 enddo
 END_PROVIDER


