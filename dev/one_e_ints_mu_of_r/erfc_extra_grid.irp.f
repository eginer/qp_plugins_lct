
BEGIN_PROVIDER [ double precision, erfc_extra_grid_rk, ( n_points_extra_final_grid , ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! erfc_extra_grid_rk(j,l,R) = int dr phi_j(r) phi_l(r) (erf(mu(R) |r - R| - 1)/|r - R|
 END_DOC
 integer :: j,l,ipoint
 double precision :: mu,r(3),NAI_pol_mult_erf_ao
 double precision :: int_mu, int_coulomb
 double precision :: thr
 thr = 1.d-12
 provide mu_erf final_grid_points_extra 
 double precision :: wall0, wall1
 call wall_time(wall0)
 provide mu_of_r_for_ints
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (j,l,ipoint,mu,r,int_mu,int_coulomb) & 
 !$OMP SHARED (aos_in_r_array_transp,ao_num,n_points_extra_final_grid,mu_of_r_for_ints,erfc_extra_grid_rk,final_grid_points_extra,thr)
 !$OMP DO SCHEDULE (dynamic)
 do l = 1, ao_num
  do j = l, ao_num
   do ipoint = 1, n_points_extra_final_grid
     erfc_extra_grid_rk(ipoint,j,l) = 0.d0
     if(dabs( aos_in_r_array_transp(ipoint,j) * aos_in_r_array_transp(ipoint,l)).lt.thr)cycle
     mu = mu_of_r_for_ints(ipoint,1)
     r(1) = final_grid_points_extra(1,ipoint)
     r(2) = final_grid_points_extra(2,ipoint)
     r(3) = final_grid_points_extra(3,ipoint)
     int_mu = NAI_pol_mult_erf_ao(j,l,mu,r)
     int_coulomb = NAI_pol_mult_erf_ao(j,l,1.d+9,r)
     erfc_extra_grid_rk(ipoint,j,l)= (int_mu - int_coulomb )
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 do l = 1, ao_num
  do j = 1, l-1
   do ipoint = 1, n_points_extra_final_grid
     erfc_extra_grid_rk(ipoint,j,l)= erfc_extra_grid_rk(ipoint,l,j)
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for erfc_extra_grid_rk  ',wall1 - wall0
END_PROVIDER 


BEGIN_PROVIDER [double precision, erf_mu_r12_inv_r12_rk_extra,(n_points_extra_final_grid , ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! erf_mu_r12_inv_r12_rk_extra(j,l,r) = \int dr' phi_j(r') phi_l(r') erf(mu(r) |r - r'|)/|r-r'|
 END_DOC
 double precision :: mu,r(3),int_mu,delta,wall0,wall1,phi_j_erf_mu_r_phi,thr
 thr = 1.d-12
 integer :: l,j,ipoint
  provide mu_of_r_for_ints
 call wall_time(wall0)
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (l,j,ipoint,mu,r,int_mu,delta) & 
  !$OMP SHARED (ao_num,aos_in_r_array_transp,n_points_extra_final_grid,mu_of_r_for_ints,erf_mu_r12_inv_r12_rk_extra,final_grid_points_extra,thr)
  !$OMP DO SCHEDULE (dynamic)
 do l = 1, ao_num
  do j = l, ao_num
   do ipoint = 1, n_points_extra_final_grid
     mu = mu_of_r_for_ints(ipoint,1)
     erf_mu_r12_inv_r12_rk_extra(ipoint,j,l) = 0.d0
     if(dabs( aos_in_r_array_transp(ipoint,j) * aos_in_r_array_transp(ipoint,l)).lt.thr)cycle
     r(1) = final_grid_points_extra(1,ipoint)
     r(2) = final_grid_points_extra(2,ipoint)
     r(3) = final_grid_points_extra(3,ipoint)
     int_mu = phi_j_erf_mu_r_phi(j,l,mu, r)
     erf_mu_r12_inv_r12_rk_extra(ipoint,j,l) = int_mu 
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

 do l = 1, ao_num
  do j = 1, l-1
   do ipoint = 1, n_points_extra_final_grid
    erf_mu_r12_inv_r12_rk_extra(ipoint,j,l)= erf_mu_r12_inv_r12_rk_extra(ipoint,l,j)
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for erf_mu_r12_inv_r12_rk_extra  ',wall1 - wall0

END_PROVIDER 

