 BEGIN_PROVIDER [double precision, ecmd_lda_mu_of_r, (N_states)]                                                                      
 BEGIN_DOC
! ecmd_lda_mu_of_r = multi-determinantal Ecmd within the LDA approximation with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none
 integer :: ipoint,istate
 double precision :: rho_a, rho_b, ec
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing ecmd_lda_mu_of_r ...'

 ecmd_lda_mu_of_r = 0.d0
 call wall_time(wall0)
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   mu     = mu_of_r_prov(ipoint,istate)
   weight = final_weight_at_r_vector(ipoint)
   rho_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   call ESRC_MD_LDAERF (mu,rho_a,rho_b,dospin,ec)
   if(isnan(ec))then
    print*,'ec is nan'
    stop
   endif
   ecmd_lda_mu_of_r(istate) += weight * ec
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for ecmd_lda_mu_of_r :',wall1-wall0
 END_PROVIDER 

