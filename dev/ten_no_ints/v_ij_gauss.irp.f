
BEGIN_PROVIDER [ double precision, v_ij_gauss_rk, ( ao_num, ao_num,n_points_final_grid,n_max_fit_ten_no_slat)]
 implicit none
 BEGIN_DOC
! int dr phi_i(r) phi_j(r) exp(-alpha_ten_no * |r - R|)
 END_DOC
 integer :: i,j,ipoint,pp
 double precision :: r(3),delta
 double precision :: int_gauss,overlap_gauss_r12_ao
 double precision :: wall0, wall1
 provide final_grid_points 
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,r,int_gauss,pp,delta) & 
 !$OMP SHARED (ao_num,n_points_final_grid,v_ij_gauss_rk,final_grid_points,n_max_fit_ten_no_slat,expo_fit_ten_no_slat_gauss)
 !$OMP DO SCHEDULE (dynamic)
 do pp = 1, n_max_fit_ten_no_slat
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
     do j = i, ao_num
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      delta = expo_fit_ten_no_slat_gauss(pp)
      int_gauss = overlap_gauss_r12_ao(r,delta,i,j)
      v_ij_gauss_rk(j,i,ipoint,pp)= int_gauss
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 do pp = 1, n_max_fit_ten_no_slat
  do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = 1, i-1
     v_ij_gauss_rk(j,i,ipoint,pp)= v_ij_gauss_rk(i,j,ipoint,pp)
    enddo
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for v_ij_gauss_rk ',wall1 - wall0
END_PROVIDER 

BEGIN_PROVIDER [ double precision, x_v_ij_gauss_rk, (ao_num, ao_num,n_points_final_grid,3,n_max_fit_ten_no_slat)]
 implicit none
 BEGIN_DOC
! int dr x * phi_i(r) phi_j(r) gauss(mu(R) |r - R|)/|r - R|
 END_DOC
 integer :: i,j,ipoint,m,pp
 double precision :: r(3),ints,delta
 double precision :: wall0, wall1
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,r,m,ints,pp,delta) & 
 !$OMP SHARED (ao_num,n_points_final_grid,x_v_ij_gauss_rk,final_grid_points,n_max_fit_ten_no_slat,expo_fit_ten_no_slat_gauss)
 !$OMP DO SCHEDULE (dynamic)
 do pp = 1, n_max_fit_ten_no_slat
  do m = 1, 3
   do ipoint = 1, n_points_final_grid
     do i = 1, ao_num
      do j = i, ao_num
       r(1) = final_grid_points(1,ipoint)
       r(2) = final_grid_points(2,ipoint)
       r(3) = final_grid_points(3,ipoint)
       delta = expo_fit_ten_no_slat_gauss(pp)
       call gauss_int_x_ao(i,j,delta,r,m,ints)
       x_v_ij_gauss_rk(j,i,ipoint,m,pp) =  ints 
     enddo
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 do pp = 1, n_max_fit_ten_no_slat
  do m = 1, 3
   do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
     do j = 1, i-1
       x_v_ij_gauss_rk(j,i,ipoint,m,pp)= x_v_ij_gauss_rk(i,j,ipoint,m,pp)
     enddo
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall1)
 print*,'wall time for x_v_ij_gauss_rk',wall1 - wall0
END_PROVIDER 


BEGIN_PROVIDER [ double precision, v_ij_gauss_rk_dble_alpha, ( ao_num, ao_num,n_points_final_grid,n_max_fit_ten_no_slat,n_max_fit_ten_no_slat)]
 implicit none
 BEGIN_DOC
! int dr phi_i(r) phi_j(r) exp(-(alpha_ten_no(i)+alpha_ten_no(j)) * |r - R|)
 END_DOC
 integer :: i,j,ipoint,pp,qq
 double precision :: r(3),delta,delta2
 double precision :: int_gauss,overlap_gauss_r12_ao
 double precision :: wall0, wall1
 provide final_grid_points 
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,r,int_gauss,pp,qq,delta,delta2) & 
 !$OMP SHARED (ao_num,n_points_final_grid,v_ij_gauss_rk_dble_alpha,final_grid_points,n_max_fit_ten_no_slat,expo_fit_ten_no_slat_gauss)
 !$OMP DO SCHEDULE (dynamic)
 do qq = 1, n_max_fit_ten_no_slat
  do pp = 1, n_max_fit_ten_no_slat
   do ipoint = 1, n_points_final_grid
     do i = 1, ao_num
      do j = i, ao_num
       r(1) = final_grid_points(1,ipoint)
       r(2) = final_grid_points(2,ipoint)
       r(3) = final_grid_points(3,ipoint)
       delta = expo_fit_ten_no_slat_gauss(pp)
       delta2 = expo_fit_ten_no_slat_gauss(qq)
       delta += delta2
       int_gauss = overlap_gauss_r12_ao(r,delta,i,j)
       v_ij_gauss_rk_dble_alpha(j,i,ipoint,pp,qq)= int_gauss
     enddo
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 do qq = 1, n_max_fit_ten_no_slat
  do pp = 1, n_max_fit_ten_no_slat
   do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
     do j = 1, i-1
      v_ij_gauss_rk_dble_alpha(j,i,ipoint,pp,qq)= v_ij_gauss_rk_dble_alpha(i,j,ipoint,pp,qq)
     enddo
    enddo
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for v_ij_gauss_rk_dble_alpha ',wall1 - wall0
END_PROVIDER 

BEGIN_PROVIDER [ double precision, x_v_ij_gauss_rk_dble_alpha, (ao_num, ao_num,n_points_final_grid,3,n_max_fit_ten_no_slat,n_max_fit_ten_no_slat)]
 implicit none
 BEGIN_DOC
! int dr x * phi_i(r) phi_j(r) gauss(mu(R) |r - R|)/|r - R|
 END_DOC
 integer :: i,j,ipoint,m,pp,qq
 double precision :: r(3),ints,delta,delta2
 double precision :: wall0, wall1
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,r,m,ints,pp,qq,delta,delta2) & 
 !$OMP SHARED (ao_num,n_points_final_grid,x_v_ij_gauss_rk_dble_alpha,final_grid_points,n_max_fit_ten_no_slat,expo_fit_ten_no_slat_gauss)
 !$OMP DO SCHEDULE (dynamic)
 do qq = 1, n_max_fit_ten_no_slat
  do pp = 1, n_max_fit_ten_no_slat
   do m = 1, 3
    do ipoint = 1, n_points_final_grid
      do i = 1, ao_num
       do j = i, ao_num
        r(1) = final_grid_points(1,ipoint)
        r(2) = final_grid_points(2,ipoint)
        r(3) = final_grid_points(3,ipoint)
        delta = expo_fit_ten_no_slat_gauss(pp)
        delta2 = expo_fit_ten_no_slat_gauss(qq)
        delta += delta2
        call gauss_int_x_ao(i,j,delta,r,m,ints)
        x_v_ij_gauss_rk_dble_alpha(j,i,ipoint,m,pp,qq) =  ints 
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 do qq = 1, n_max_fit_ten_no_slat
  do pp = 1, n_max_fit_ten_no_slat
   do m = 1, 3
    do ipoint = 1, n_points_final_grid
     do i = 1, ao_num
      do j = 1, i-1
        x_v_ij_gauss_rk_dble_alpha(j,i,ipoint,m,pp,qq)= x_v_ij_gauss_rk_dble_alpha(i,j,ipoint,m,pp,qq)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall1)
 print*,'wall time for x_v_ij_gauss_rk_dble_alpha',wall1 - wall0
END_PROVIDER 


subroutine gauss_int_x_ao(i_ao,j_ao,delta,C_center,m,ints)
 implicit none
  BEGIN_DOC
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) exp(-delta * |r-C_center|^2)$.
  END_DOC
 include 'utils/constants.include.F'                                                                                                                                  
 integer, intent(in) :: i_ao,j_ao,m
 double precision, intent(in) :: delta, C_center(3)
 double precision, intent(out):: ints
 double precision               :: overlap_gauss_r12
 double precision               :: A_center(3), B_center(3),integral, alpha,beta
 integer                        :: i,j,num_A,num_B, power_A(3), power_B(3), n_pt_in, power_xA(3)
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
    power_xA = power_A
    ! x * phi_i(r) = x * (x-Ax)**ax = (x-Ax)**(ax+1) + Ax * (x-Ax)**ax
    power_xA(m) += 1
    do j = 1, ao_prim_num(j_ao)
      beta = ao_expo_ordered_transp(j,j_ao)
      ! First term = (x-Ax)**(ax+1)
      integral = overlap_gauss_r12(C_center,delta,A_center,B_center,power_xA,power_B,alpha,beta)
      ints += integral * ao_coef_normalized_ordered_transp(j,j_ao)*ao_coef_normalized_ordered_transp(i,i_ao)
      ! Second term = Ax * (x-Ax)**(ax)
      integral = overlap_gauss_r12(C_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
      ints += A_center(m) * integral * ao_coef_normalized_ordered_transp(j,j_ao)*ao_coef_normalized_ordered_transp(i,i_ao)
    enddo
 enddo
end


