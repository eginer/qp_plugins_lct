BEGIN_PROVIDER [ double precision, v_ij_erf_rk, ( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr phi_i(r) phi_j(r) erf(mu(R) |r - R|)/|r - R|
 END_DOC
 integer :: i,j,ipoint
 double precision :: mu,r(3),NAI_pol_mult_erf_ao
 double precision :: int_mu, int_coulomb
 provide mu_erf final_grid_points constant_mu 
 if(constant_mu)then
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,int_coulomb) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_erf,v_ij_erf_rk,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_erf 
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = NAI_pol_mult_erf_ao(i,j,mu,r)
!     int_coulomb = NAI_pol_mult_erf_ao(i,j,1.d+9,r)
     int_coulomb = 0.d0
     v_ij_erf_rk(j,i,ipoint)= (int_mu - int_coulomb )
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 else 
  provide mu_of_r_prov
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,int_coulomb) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_prov,v_ij_erf_rk,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_prov(ipoint,1)
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = NAI_pol_mult_erf_ao(i,j,mu,r)
!     int_coulomb = NAI_pol_mult_erf_ao(i,j,1.d+9,r)
     int_coulomb = 0.d0
     v_ij_erf_rk(j,i,ipoint)= (int_mu - int_coulomb )
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 endif

 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
     v_ij_erf_rk(j,i,ipoint)= v_ij_erf_rk(i,j,ipoint)
    enddo
   enddo
  enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, x_v_ij_erf_rk, (3,ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr x * phi_i(r) phi_j(r) erf(mu(R) |r - R|)/|r - R|
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: mu,r(3),ints(3),ints_coulomb(3)
 if(constant_mu)then
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,ints,m,ints_coulomb) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_erf,x_v_ij_erf_rk,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_erf 
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     call NAI_pol_x_mult_erf_ao(i,j,mu,r,ints)
!     call NAI_pol_x_mult_erf_ao(i,j,1.d+9,r,ints_coulomb)
     ints_coulomb = 0.d0
     do m = 1, 3
      x_v_ij_erf_rk(m,j,i,ipoint)= ( ints(m) - ints_coulomb(m))
     enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 else
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,ints,m,ints_coulomb) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_prov,x_v_ij_erf_rk,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_prov(ipoint,1)
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     call NAI_pol_x_mult_erf_ao(i,j,mu,r,ints)
!     call NAI_pol_x_mult_erf_ao(i,j,1.d+9,r,ints_coulomb)
     ints_coulomb = 0.d0
     do m = 1, 3
      x_v_ij_erf_rk(m,j,i,ipoint)= ( ints(m) - ints_coulomb(m))
     enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 endif
 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
     x_v_ij_erf_rk(:,j,i,ipoint)= x_v_ij_erf_rk(:,i,j,ipoint)
    enddo
   enddo
  enddo


END_PROVIDER 



subroutine NAI_pol_x_mult_erf_ao(i_ao,j_ao,mu_in,C_center,ints)
 implicit none
  BEGIN_DOC
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr y * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr z * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  END_DOC
 include 'utils/constants.include.F'                                                                                                                                  
 integer, intent(in) :: i_ao,j_ao
 double precision, intent(in) :: mu_in, C_center(3)
 double precision, intent(out):: ints(3)
 double precision               :: A_center(3), B_center(3),integral, alpha,beta
 double precision               :: NAI_pol_mult_erf
 integer                        :: i,j,num_A,num_B, power_A(3), power_B(3), n_pt_in, power_xA(3),m
 num_A = ao_nucl(i_ao)
 power_A(1:3)= ao_power(i_ao,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j_ao)
 power_B(1:3)= ao_power(j_ao,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 n_pt_in = n_pt_max_integrals

 ints = 0.d0

 do i = 1, ao_prim_num(i_ao)
  alpha = ao_expo_ordered_transp(i,i_ao)
   do m = 1, 3
    power_xA = power_A
    ! x * phi_i(r) = x * (x-Ax)**ax = (x-Ax)**(ax+1) + Ax * (x-Ax)**ax
    power_xA(m) += 1
    do j = 1, ao_prim_num(j_ao)
      beta = ao_expo_ordered_transp(j,j_ao)
      ! First term = (x-Ax)**(ax+1)
      integral =  NAI_pol_mult_erf(A_center,B_center,power_xA,power_B,alpha,beta,C_center,n_pt_in,mu_in)
      ints(m) += integral * ao_coef_normalized_ordered_transp(j,j_ao)*ao_coef_normalized_ordered_transp(i,i_ao)
      ! Second term = Ax * (x-Ax)**(ax)
      integral =  NAI_pol_mult_erf(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in,mu_in)
      ints(m) += A_center(m) * integral * ao_coef_normalized_ordered_transp(j,j_ao)*ao_coef_normalized_ordered_transp(i,i_ao)
    enddo
  enddo
 enddo
end

BEGIN_PROVIDER [ double precision, x_v_ij_erf_rk_trans, (ao_num, ao_num,n_points_final_grid,3)]
 implicit none
 BEGIN_DOC
! int dr x * phi_i(r) phi_j(r) erf(mu(R) |r - R|)/|r - R|
 END_DOC
 integer :: i,j,ipoint
 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, ao_num
     x_v_ij_erf_rk_trans(j,i,ipoint,:)= x_v_ij_erf_rk(:,j,i,ipoint)
    enddo
   enddo
  enddo

END_PROVIDER 



