
 BEGIN_PROVIDER [double precision, average_mu_rs_c     ]
&BEGIN_PROVIDER [double precision, mu_of_r_rs_c, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_mu_of_r_rs_c, (3,n_points_final_grid,N_states) ]
 implicit none
 integer :: ipoint,i,m
 include 'constants.include.F'
 double precision :: weight, rho_a_hf, rho_b_hf, rho_hf, mu_rs_c,r(3)
 average_mu_rs_c    = 0.d0
 double precision :: elec_a,elec_b,grad_mu(3)
 double precision :: dx,mu_lda, mu_tmp,mu_min,mu,damped_mu
 dx = 1.d-5
 elec_a = 0.d0
 elec_b = 0.d0
 mu_min = mu_erf
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  r(:) = final_grid_points_extra(:,ipoint)
  rho_a_hf = one_e_dm_and_grad_alpha_in_r(4,ipoint,1)
  rho_b_hf = one_e_dm_and_grad_beta_in_r(4,ipoint,1)
  rho_hf = rho_a_hf + rho_b_hf
  mu_tmp = mu_rs_c(rho_hf)
  if(damped_mu_of_r)then
   mu = damped_mu(mu_tmp,mu_min)
  else
   mu = max(mu_tmp,mu_of_r_min)
  endif
  mu_of_r_rs_c(ipoint,1) = mu
  average_mu_rs_c   +=  mu_of_r_rs_c(ipoint,1) * rho_hf * weight
  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight

  if(damped_mu_of_r)then
   call get_grad_damped_mu_rsc(r,mu_min,dx,grad_mu) 
  else
   call get_grad_mu_rsc(r,dx,grad_mu)
  endif
  grad_mu_of_r_rs_c(:,ipoint,1) = grad_mu(:)

 enddo 
 average_mu_rs_c   = average_mu_rs_c  / dble(elec_a+ elec_b)

END_PROVIDER 



double precision function mu_rs_c(rho)
 implicit none
 double precision, intent(in) :: rho
 include 'constants.include.F'
 double precision :: cst_rs,alpha_rs,rs
 cst_rs   = (4.d0 * dacos(-1.d0)/3.d0)**(-1.d0/3.d0)
 alpha_rs = 2.d0 * dsqrt((9.d0 * dacos(-1.d0)/4.d0)**(-1.d0/3.d0)) / sqpi
  
 rs = cst_rs * rho**(-1.d0/3.d0)
 mu_rs_c =  alpha_rs/dsqrt(rs) 

end

double precision function damped_mu_rs_c(rho,mu_min)
 implicit none
 double precision, intent(in) :: rho,mu_min
 double precision :: mu_rs_c,mu_tmp,damped_mu
 mu_tmp = mu_rs_c(rho)
 damped_mu_rs_c = damped_mu(mu_tmp, mu_min)
end

subroutine get_grad_mu_rsc(r,dx,grad_mu)
 implicit none
 double precision, intent(in) :: r(3),dx
 double precision, intent(out):: grad_mu(3)
 integer :: m
 double precision :: r_tmp(3),mu_plus, mu_minus,mu_rs_c,rho_a_hf,rho_b_hf,rho_hf
  do m = 1, 3
   r_tmp = r 
   r_tmp(m) += dx 
   call dm_dft_alpha_beta_at_r(r_tmp,rho_a_hf,rho_b_hf)
   rho_hf = rho_a_hf + rho_b_hf
   mu_plus = mu_rs_c(rho_hf)
   r_tmp = r 
   r_tmp(m) -= dx 
   call dm_dft_alpha_beta_at_r(r_tmp,rho_a_hf,rho_b_hf)
   rho_hf = rho_a_hf + rho_b_hf
   mu_minus = mu_rs_c(rho_hf)

   grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
  enddo
end

subroutine get_grad_damped_mu_rsc(r,mu_min,dx,grad_mu)
 implicit none
 double precision, intent(in) :: r(3),mu_min,dx
 double precision, intent(out):: grad_mu(3)
 integer :: m
 double precision :: r_tmp(3),mu_plus, mu_minus,mu_rs_c
 double precision :: rho_a_hf,rho_b_hf,rho_hf,damped_mu,rho,damped_mu_rs_c
  do m = 1, 3
   r_tmp = r 
   r_tmp(m) += dx 
   call dm_dft_alpha_beta_at_r(r_tmp,rho_a_hf,rho_b_hf)
   rho_hf = rho_a_hf + rho_b_hf
   mu_plus = damped_mu_rs_c(rho,mu_min)

   r_tmp = r 
   r_tmp(m) -= dx 
   call dm_dft_alpha_beta_at_r(r_tmp,rho_a_hf,rho_b_hf)
   rho_hf = rho_a_hf + rho_b_hf
   mu_minus = damped_mu_rs_c(rho,mu_min)

   grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
  enddo
end


 BEGIN_PROVIDER [double precision, mu_of_r_extra_grid_rs_c, (n_points_extra_final_grid,N_states) ]
 implicit none
 integer :: ipoint,i,m
 include 'constants.include.F'
 double precision :: weight, rho_a_hf, rho_b_hf, rho_hf, mu_rs_c
 average_mu_rs_c    = 0.d0
 double precision :: elec_a,elec_b,r(3)
 double precision :: mu_tmp ,damped_mu, mu ,dx,mu_min
 dx = 1.d-5
 elec_a = 0.d0
 elec_b = 0.d0
 mu_min = mu_erf 
 do ipoint = 1, n_points_extra_final_grid
  weight = final_weight_at_r_vector_extra(ipoint)
  r(:) = final_grid_points_extra(:,ipoint)
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  rho_hf = rho_a_hf + rho_b_hf
  mu_tmp = mu_rs_c(rho_hf)
  if(damped_mu_of_r)then
   mu = damped_mu(mu_tmp,mu_min)
  else
   mu = max(mu_tmp,mu_of_r_min)
  endif
  mu_of_r_extra_grid_rs_c(ipoint,1) = mu
  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight
 enddo 

END_PROVIDER 


 BEGIN_PROVIDER [double precision, mu_of_r_extra_grid_hf, (n_points_extra_final_grid) ]
 implicit none 
 BEGIN_DOC
 ! mu(r) computed with a HF wave function (assumes that HF MOs are stored in the EZFIO)
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) but for \Psi^B = HF^B
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals 
 !
 ! in the two-body density matrix are excluded
 END_DOC
 integer :: ipoint
 double precision :: wall0,wall1,f_hf,on_top,w_hf,sqpi,r(3)
 PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals 
 print*,'providing mu_of_r_extra_grid_hf ...'
 call wall_time(wall0)
 sqpi = dsqrt(dacos(-1.d0))
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,f_hf,on_top,w_hf,r) & 
 !$OMP ShARED (n_points_extra_final_grid,mu_of_r_extra_grid_hf,sqpi,final_grid_points_extra) 
 do ipoint = 1, n_points_extra_final_grid
  r(:) = final_grid_points_extra(:,ipoint)
  call f_HF_valence_ab(r,r,f_hf,on_top)
  if(on_top.le.1.d-12.or.f_hf.le.0.d0.or.f_hf * on_top.lt.0.d0)then
    w_hf   = 1.d+10
  else 
    w_hf  = f_hf /  on_top
  endif
  mu_of_r_extra_grid_hf(ipoint) =  w_hf * sqpi * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall1)
 print*,'Time to provide mu_of_r_extra_grid_hf = ',wall1-wall0
 END_PROVIDER 





