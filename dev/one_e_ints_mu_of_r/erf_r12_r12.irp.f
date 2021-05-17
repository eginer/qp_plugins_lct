
BEGIN_PROVIDER [double precision, erf_mu_r12_inv_r12_rk,( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! erf_mu_r12_inv_r12_rk(i,j,r) = \int dr' phi_i(r') phi_j(r') erf(mu(r) |r - r'|)/|r-r'|
 END_DOC
 double precision :: mu,r(3),int_mu,delta,wall0,wall1,phi_j_erf_mu_r_phi
 integer :: i,j,ipoint
  provide mu_of_r_for_ints
 call wall_time(wall0)
 if(constant_mu)then
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta) & 
  !$OMP SHARED (ao_num,n_points_final_grid,mu_erf,erf_mu_r12_inv_r12_rk,final_grid_points)
  !$OMP DO SCHEDULE (dynamic)
  do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_erf
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = phi_j_erf_mu_r_phi(i,j,mu, r)
     erf_mu_r12_inv_r12_rk(j,i,ipoint) = int_mu 
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
 else
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta) & 
  !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,erf_mu_r12_inv_r12_rk,final_grid_points)
  !$OMP DO SCHEDULE (dynamic)
  do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_for_ints(ipoint,1)
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = phi_j_erf_mu_r_phi(i,j,mu, r)
     erf_mu_r12_inv_r12_rk(j,i,ipoint) = int_mu 
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif

 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
    erf_mu_r12_inv_r12_rk(j,i,ipoint)= erf_mu_squared_ij_rk(i,j,ipoint)
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for erf_mu_r12_inv_r12_rk  ',wall1 - wall0

END_PROVIDER 


BEGIN_PROVIDER [double precision, coulomb_rk,( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! coulomb_rk(i,j,r) = \int dr' phi_i(r') phi_j(r') 1/|r-r'|
 END_DOC
 double precision :: r(3),int_mu,delta,wall0,wall1,mu,phi_j_erf_mu_r_phi
 integer :: i,j,ipoint
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,r,int_mu,delta,mu) & 
 !$OMP SHARED (ao_num,n_points_final_grid,coulomb_rk,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = i, ao_num
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)
    mu = 1.d+10
    int_mu = phi_j_erf_mu_r_phi(i,j,mu, r)
    coulomb_rk(j,i,ipoint) = int_mu
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL


 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
    coulomb_rk(j,i,ipoint)= erf_mu_squared_ij_rk(i,j,ipoint)
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for coulomb_rk  ',wall1 - wall0

END_PROVIDER 


