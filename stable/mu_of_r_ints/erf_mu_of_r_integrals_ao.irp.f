double precision function erf_mu_of_r_ao_old(i,j,k,l)
 implicit none
 BEGIN_DOC
! computes the following integral :
! int dr1 dr2 AO_i(r1) AO_j(r1) erf(mu(r1) |r1-r2|) / |r1-r2| AO_k(r2) AO_l(r2)
 END_DOC
 integer, intent(in) :: i,j,k,l
 integer :: i_atom,i_radial,i_angular
 double precision, allocatable :: aos_array(:)
 double precision :: r(3),NAI_pol_mult_erf_ao
 allocate(aos_array(ao_num))
 erf_mu_of_r_ao_old = 0.d0
  do i_atom = 1,nucl_num
   do i_radial = 1, n_points_radial_grid -1
    do i_angular = 1 , n_points_integration_angular
     r(1) = grid_points_per_atom(1,i_angular,i_radial,i_atom)
     r(2) = grid_points_per_atom(2,i_angular,i_radial,i_atom)
     r(3) = grid_points_per_atom(3,i_angular,i_radial,i_atom)
     call give_all_aos_at_r(r,aos_array)
     erf_mu_of_r_ao_old += aos_array(i) * aos_array(j) * NAI_pol_mult_erf_ao(k,l,mu_of_r_prov(i_angular,i_radial,i_atom),r) * final_weight_functions_at_grid_points(i_angular,i_radial,i_atom)
    enddo
   enddo
  enddo
 deallocate(aos_array)

end

double precision function erf_mu_of_r_ao(i,j,k,l)
 implicit none
 BEGIN_DOC
! computes the following integral :
! int dr1 dr2 AO_i(r1) AO_j(r1) erf(mu(r1) |r1-r2|) / |r1-r2| AO_k(r2) AO_l(r2)
 END_DOC
 integer, intent(in) :: i,j,k,l
 integer :: i_point
 double precision, allocatable :: aos_array(:)
 double precision :: r(3),NAI_pol_mult_erf_ao
 allocate(aos_array(ao_num))
 erf_mu_of_r_ao = 0.d0
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point) 
  r(2) = final_grid_points(2,i_point) 
  r(3) = final_grid_points(3,i_point) 
  call give_all_aos_at_r(r,aos_array)
  erf_mu_of_r_ao += 0.5d0 * aos_array(i) * aos_array(j) * NAI_pol_mult_erf_ao(k,l,mu_of_r_prov_selected(i_point),r) * final_weight_functions_at_final_grid_points(i_point)
  erf_mu_of_r_ao += 0.5d0 * aos_array(k) * aos_array(l) * NAI_pol_mult_erf_ao(j,i,mu_of_r_prov_selected(i_point),r) * final_weight_functions_at_final_grid_points(i_point)
 enddo
 deallocate(aos_array)

end


 subroutine give_all_erf_mu_of_r_kl(k,l,integrals)
 implicit none 
 include 'Utils/constants.include.F'
 integer, intent(in) :: k,l
 double precision, intent(out) :: integrals(ao_num,ao_num)
 integer :: i_point,i,j
 double precision :: integrals_kl_of_r(n_points_final_grid),r(3),NAI_pol_mult_erf_ao,tmp
 double precision,allocatable :: v_array(:,:),v_vector(:)
 allocate(v_array(ao_num,n_points_final_grid),v_vector(n_points_final_grid))
 integrals = 0.d0
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point) 
  r(2) = final_grid_points(2,i_point) 
  r(3) = final_grid_points(3,i_point) 
  tmp = NAI_pol_mult_erf_ao(k,l,mu_of_r_prov_selected(i_point),r) 
  v_vector(i_point) = tmp
  do i = 1, ao_num
   v_array(i,i_point) = tmp * final_weight_functions_at_final_grid_points(i_point) * aos_in_r_array(i,i_point)
  enddo
 enddo
!call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,v_array,size(v_array,1),aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,integrals,size(integrals,1))
 do i_point = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, ao_num
    integrals(j,i) += v_vector(i_point) * aos_in_r_array(i,i_point) * aos_in_r_array(j,i_point) * final_weight_functions_at_final_grid_points(i_point)
   enddo
  enddo
 enddo
 end

