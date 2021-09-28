BEGIN_PROVIDER [double precision, denom_non_ad_coupling , (n_states,n_states)]
 implicit none
 integer :: i,j,m
 double precision :: delta_mu(3),d_mu,mu_12(3), mu_12_sq
 denom_non_ad_coupling = 0.d0
 do i = 1, n_states
  do j = 1, n_states
   ! mu_11 - mu_22
   delta_mu(1) = trans_dipole_x(i,i) - trans_dipole_x(j,j)  
   delta_mu(2) = trans_dipole_y(i,i) - trans_dipole_y(j,j)  
   delta_mu(3) = trans_dipole_z(i,i) - trans_dipole_z(j,j)  
   ! mu_12 
   mu_12(1) = trans_dipole_x(j,i) 
   mu_12(2) = trans_dipole_y(j,i) 
   mu_12(3) = trans_dipole_z(j,i) 
   d_mu = 0.d0
   mu_12_sq = 0.d0
   do m = 1, 3
    d_mu     += delta_mu(m) * delta_mu(m)
    mu_12_sq += mu_12(m)    * mu_12(m)
   enddo
   denom_non_ad_coupling(j,i) = dsqrt(d_mu + 4.d0 * mu_12_sq)
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, non_ad_coupling, (n_states,n_states)]
 implicit none
 integer :: i,j,m
 double precision :: delta_mu(3),d_mu,mu_12(3), mu_12_norm, delta_E
 non_ad_coupling = 0.d0
 do i = 1, n_states
  do j = 1, n_states
   if(i==j)cycle
   ! mu_12 
   mu_12(1) = trans_dipole_x(j,i) 
   mu_12(2) = trans_dipole_y(j,i) 
   mu_12(3) = trans_dipole_z(j,i) 
   mu_12_norm = 0.d0
   do m = 1, 3
    mu_12_norm += mu_12(m) * mu_12(m)
   enddo
   mu_12_norm = dsqrt(mu_12_norm)
   delta_E = ci_electronic_energy(j) - ci_electronic_energy(i)
   non_ad_coupling(j,i) = delta_E * mu_12_norm / denom_non_ad_coupling(j,i)
  enddo
 enddo
END_PROVIDER 
