
BEGIN_PROVIDER [ double precision, erfc_extra_grid_rk, ( ao_num, ao_num,n_points_extra_final_grid)]
 implicit none
 BEGIN_DOC
! int dr phi_i(r) phi_j(r) (erf(mu(R) |r - R| - 1)/|r - R|
 END_DOC
 integer :: i,j,ipoint
 double precision :: mu,r(3),NAI_pol_mult_erf_ao
 double precision :: int_mu, int_coulomb
 provide mu_erf final_grid_points_extra 
 double precision :: wall0, wall1
 call wall_time(wall0)
 provide mu_of_r_for_ints
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,int_coulomb) & 
 !$OMP SHARED (ao_num,n_points_extra_final_grid,mu_of_r_for_ints,erfc_extra_grid_rk,final_grid_points_extra)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_extra_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_for_ints(ipoint,1)
     r(1) = final_grid_points_extra(1,ipoint)
     r(2) = final_grid_points_extra(2,ipoint)
     r(3) = final_grid_points_extra(3,ipoint)
     int_mu = NAI_pol_mult_erf_ao(i,j,mu,r)
     int_coulomb = NAI_pol_mult_erf_ao(i,j,1.d+9,r)
     erfc_extra_grid_rk(j,i,ipoint)= (int_mu - int_coulomb )
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 do ipoint = 1, n_points_extra_final_grid
  do i = 1, ao_num
   do j = 1, i-1
     erfc_extra_grid_rk(j,i,ipoint)= erfc_extra_grid_rk(i,j,ipoint)
    enddo
   enddo
  enddo

 call wall_time(wall1)
 print*,'wall time for erfc_extra_grid_rk  ',wall1 - wall0
END_PROVIDER 


BEGIN_PROVIDER [double precision, erf_mu_r12_inv_r12_rk_extra,( ao_num, ao_num,n_points_extra_final_grid)]
 implicit none
 BEGIN_DOC
! erf_mu_r12_inv_r12_rk_extra(i,j,r) = \int dr' phi_i(r') phi_j(r') erf(mu(r) |r - r'|)/|r-r'|
 END_DOC
 double precision :: mu,r(3),int_mu,delta,wall0,wall1,phi_j_erf_mu_r_phi
 integer :: i,j,ipoint
  provide mu_of_r_for_ints
 call wall_time(wall0)
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta) & 
  !$OMP SHARED (ao_num,n_points_extra_final_grid,mu_of_r_for_ints,erf_mu_r12_inv_r12_rk_extra,final_grid_points_extra)
  !$OMP DO SCHEDULE (dynamic)
  do ipoint = 1, n_points_extra_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_for_ints(ipoint,1)
     r(1) = final_grid_points_extra(1,ipoint)
     r(2) = final_grid_points_extra(2,ipoint)
     r(3) = final_grid_points_extra(3,ipoint)
     int_mu = phi_j_erf_mu_r_phi(i,j,mu, r)
     erf_mu_r12_inv_r12_rk_extra(j,i,ipoint) = int_mu 
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

 do ipoint = 1, n_points_extra_final_grid
  do i = 1, ao_num
   do j = 1, i-1
    erf_mu_r12_inv_r12_rk_extra(j,i,ipoint)= erf_mu_r12_inv_r12_rk_extra(i,j,ipoint)
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for erf_mu_r12_inv_r12_rk_extra  ',wall1 - wall0

END_PROVIDER 

