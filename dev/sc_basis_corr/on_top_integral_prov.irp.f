!mu :
!mu_of_r_prov            (ipoint, istate) -> sqrt(pi)/pi* W_PsiB
!!!!!!!
!rho_2 :
!mu_correction_of_on_top (ipoint, istate) -> extrapolation 
!on_top_cas_mu_r         (ipoint, istate) -> 

!---------------------------------------------------------------------
BEGIN_PROVIDER[double precision, n2_on_top_int, (N_states)]
&BEGIN_PROVIDER[double precision, n2_on_top_extrapolated_int, (N_states)]
 implicit none
 BEGIN_DOC
 ! Give the integral of on top pair density
 END_DOC

  integer          :: ipoint, istate
  double precision :: weight
  double precision :: rho2, rho2_ex, mu
  double precision :: mu_correction_of_on_top

  n2_on_top_int = 0.d0
  n2_on_top_extrapolated_int = 0.d0

  do istate=1, N_states
   do ipoint=1, n_points_final_grid
    weight = final_weight_at_r_vector(ipoint)
    mu = mu_of_r_prov(ipoint,istate)

    rho2 = on_top_cas_mu_r(ipoint,istate)
    rho2_ex = mu_correction_of_on_top(mu,rho2)

    n2_on_top_int(istate)              += weight*rho2
    n2_on_top_extrapolated_int(istate) += weight*rho2_ex  
   enddo 
  enddo
END_PROVIDER
