
 BEGIN_PROVIDER[double precision, aos_in_r_extra_grid_array, (ao_num,n_points_extra_final_grid)]
 implicit none
 BEGIN_DOC
 ! aos_in_r_extra_grid_array(i,j)        = value of the ith ao on the jth grid point
 END_DOC
 integer :: i,j
 double precision :: aos_array(ao_num), r(3)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i,r,aos_array,j) & 
 !$OMP SHARED(aos_in_r_extra_grid_array,n_points_extra_final_grid,ao_num,final_grid_points_extra)
 do i = 1, n_points_extra_final_grid
  r(1) = final_grid_points_extra(1,i)
  r(2) = final_grid_points_extra(2,i)
  r(3) = final_grid_points_extra(3,i)
  call give_all_aos_at_r(r,aos_array)
  do j = 1, ao_num
   aos_in_r_extra_grid_array(j,i) = aos_array(j)
  enddo
 enddo
 !$OMP END PARALLEL DO

 END_PROVIDER


 BEGIN_PROVIDER[double precision, aos_in_r_extra_grid_array_transp, (n_points_extra_final_grid,ao_num)]
 implicit none
 BEGIN_DOC
 ! aos_in_r_extra_grid_array_transp(i,j) = value of the jth ao on the ith grid point
 END_DOC
 integer :: i,j
 double precision :: aos_array(ao_num), r(3)
 do i = 1, n_points_extra_final_grid
  do j = 1, ao_num
   aos_in_r_extra_grid_array_transp(i,j) = aos_in_r_extra_grid_array(j,i) 
  enddo
 enddo

 END_PROVIDER



subroutine all_erf_mu_r1_lr_int_big_mat(big_mat)
 implicit none
! BEGIN_DOC
! enters with the array big_mat and add the following integrals 
!
! big_mat(i,k,j,l) += \int dr1 \phi_i(r1) \phi_k(r1) \int dr2 erf(\mu(r1) r12)/r12 \phi_j(r2) \phi_l(r2)
! END_DOC
 include 'constants.include.F'
 integer :: i,k,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
 double precision, intent(inout) :: big_mat(ao_num,ao_num,ao_num,ao_num)
 

 print*,'computing all_erf_mu_r1_lr_int_big_mat ...'
 call wall_time(wall0)
 double precision :: wall0,wall1

 allocate(a_mat(ao_num,ao_num,n_points_extra_final_grid))

 do ipoint = 1, n_points_extra_final_grid
  weight = final_weight_at_r_vector_extra(ipoint)
  do k = 1, ao_num
   do i = 1, ao_num
    a_mat(i,k, ipoint) = aos_in_r_extra_grid_array(k,ipoint) * aos_in_r_extra_grid_array(i,ipoint) * weight
   enddo
  enddo
 enddo

 ! erf(mu(r) * r12)/r12
 call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_extra_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                   ,erf_mu_r12_inv_r12_rk_extra(1,1,1),n_points_extra_final_grid,1.d0,big_mat,ao_num*ao_num)

 call wall_time(wall1)
 print*,'wall time for all_erf_mu_r1_lr_int_big_mat ',wall1 - wall0   

end

subroutine all_erf_mu_r1_lr_int_big_mat_bis(big_mat)
 implicit none
 BEGIN_DOC
! enters with the array big_mat and add the following integrals 
!
! big_mat(i,k,j,l) += \int dr1 \phi_i(r1) \phi_k(r1) \int dr2 erf(\mu(r1) r12)/r12 \phi_j(r2) \phi_l(r2)
 END_DOC
 include 'constants.include.F'
 integer :: i,k,j,l,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
 double precision, intent(inout) :: big_mat(ao_num,ao_num,ao_num,ao_num)

 print*,'computing all_erf_mu_r1_lr_int_big_mat ...'
 call wall_time(wall0)
 double precision :: wall0,wall1
 double precision:: ints_coulomb(ao_num)
 ! First you put the full interaction 1 /r12
 do l = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    call get_ao_two_e_integrals(j,k,l,ao_num,ints_coulomb)
    do i = 1, ao_num
     big_mat(i,k,j,l) +=  ints_coulomb(i)
    enddo
   enddo
  enddo
 enddo

 allocate(a_mat(ao_num,ao_num,n_points_extra_final_grid))

 do ipoint = 1, n_points_extra_final_grid
  weight = final_weight_at_r_vector_extra(ipoint)
  do k = 1, ao_num
   do i = 1, ao_num
    a_mat(i,k, ipoint) = aos_in_r_extra_grid_array(k,ipoint) * aos_in_r_extra_grid_array(i,ipoint) * weight
   enddo
  enddo
 enddo

 ! 1/r12 + (erf(mu(r1)) - 1)/r12
 call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_extra_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                   ,erfc_extra_grid_rk(1,1,1),n_points_extra_final_grid,1.d0,big_mat,ao_num*ao_num)

 call wall_time(wall1)
 print*,'wall time for all_erf_mu_r1_lr_int_big_mat ',wall1 - wall0   

end
