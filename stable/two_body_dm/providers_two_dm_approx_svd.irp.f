
 BEGIN_PROVIDER [integer, approx_two_dm_map, (mo_tot_num**2,2)]
&BEGIN_PROVIDER [integer, approx_two_dm_map_rev, (mo_tot_num,mo_tot_num)]
 implicit none
 BEGIN_DOC 
! map a couple of orbitals to an index  
 END_DOC
 integer :: i,j,ij
 ij = 0
 do i = 1,mo_tot_num
  do j=1,mo_tot_num
   ij += 1
   approx_two_dm_map(ij,1) = i
   approx_two_dm_map(ij,2) = j
   approx_two_dm_map_rev(i,j) = ij
  enddo
 enddo

END_PROVIDER

 BEGIN_PROVIDER [integer, n_singular_approx_svd_two_dm, (N_states)]
&BEGIN_PROVIDER [integer, n_max_singular_approx_svd_two_dm]
 BEGIN_DOC 
! Number of singular vectors (left and right) kept to represent the two-body tensor of each states
! This number depends on thr_ontop_approx 
 END_DOC
 implicit none
 integer :: i,j,k,l,ij,kl
 double precision :: thresh
 thresh= thr_ontop_approx

 print*,'************************'
 print*,'thresh     approx_two_dm    =',thresh
 print*,'************************'

 double precision, allocatable :: mat_i(:,:)
 allocate(mat_i(mo_tot_num**2,mo_tot_num**2))

 double precision, allocatable :: u_i(:,:),vt_i(:,:),D_i(:)
 allocate(u_i(mo_tot_num**2,mo_tot_num**2),vt_i(mo_tot_num**2,mo_tot_num**2),D_i(mo_tot_num**2))

 !!!!!!unfoldage!!!!!!!
 print*,'*** Approximated two-body tensor '
 integer :: istate
 do istate = 1, N_states
  ij = 0
  kl = 0
  print*,'*** For state',istate
  do i=1,mo_tot_num
   do j=1,mo_tot_num
    ij += 1
    kl = 0
    do k=1,mo_tot_num
     do l=1,mo_tot_num
      kl += 1
      mat_i(kl,ij) = two_bod_alpha_beta_mo_physician(l,k,j,i,istate)
     enddo
    enddo
   enddo
  enddo
 !!!!!!!test SVD!!!!!!!

  print*,'***** SVD of the two-body tensor *****'
  double precision :: wall_1,wall_2
  call wall_time(wall_1)
  call svd(mat_i,size(mat_i,1),u_i,size(u_i,1),D_i,vt_i,size(vt_i,1),size(mat_i,1),size(mat_i,2))
  call wall_time(wall_2)
 print*,'*****SVD is done ******' 
 print*,'Wall time for SVD = ',wall_2 - wall_1
  n_singular_approx_svd_two_dm(istate) =1 
  print*,n_singular_approx_svd_two_dm(istate),D_i(n_singular_approx_svd_two_dm(istate))
  do while ( (dabs(D_i(n_singular_approx_svd_two_dm(istate))) .gt. thresh) .AND. (n_singular_approx_svd_two_dm(istate) .lt. mo_tot_num**2)  )
   n_singular_approx_svd_two_dm(istate) += 1
   print*,n_singular_approx_svd_two_dm(istate),D_i(n_singular_approx_svd_two_dm(istate))
  enddo
  if (dabs(D_i(n_singular_approx_svd_two_dm(istate))) .lt. thresh) then
   n_singular_approx_svd_two_dm(istate) -= 1
  endif
  print*,'************************'
  print*,'n_eigen_approx_two_dm  =',n_singular_approx_svd_two_dm(istate)
  print*,'************************'
 enddo
 n_max_singular_approx_svd_two_dm = maxval(n_singular_approx_svd_two_dm)

END_PROVIDER

 BEGIN_PROVIDER [double precision, singular_approx_svd_two_dm, (n_max_singular_approx_svd_two_dm, N_states)]
&BEGIN_PROVIDER [double precision, singular_left_vec_approx_svd_two_dm, (mo_tot_num**2,n_max_singular_approx_svd_two_dm,N_states)]
&BEGIN_PROVIDER [double precision, singular_right_vec_approx_svd_two_dm, (mo_tot_num**2,n_max_singular_approx_svd_two_dm,N_states)]
 BEGIN_DOC 
! singular_approx_svd_two_dm(i,istate) = ith Singular values to represent the two-body tensor of the state istate 
! singular_left_vec_approx_svd_two_dm(I,i,istate) = component on the couple of MOs approx_two_dm_map(I,1) and approx_two_dm_map(I,2) of the ith Singular left vector of the state istate 
! singular_right_vec_approx_svd_two_dm(I,i,istate) = component on the couple of MOs approx_two_dm_map(I,1) and approx_two_dm_map(I,2) of the ith Singular right vector of the state istate 
 END_DOC
 implicit none
 integer :: i,j,k,l,ij,kl
 double precision, allocatable :: mat_i(:,:)
 allocate(mat_i(mo_tot_num**2,mo_tot_num**2))

 double precision, allocatable :: u_i(:,:),vt_i(:,:),D_i(:)
 allocate(u_i(mo_tot_num**2,mo_tot_num**2),vt_i(mo_tot_num**2,mo_tot_num**2),D_i(mo_tot_num**2))

 integer :: istate
 do istate = 1, N_states
  ij = 0
  kl = 0
  !!!!!!unfoldage!!!!!!!
  do i = 1,mo_tot_num
   do j=1,mo_tot_num
    ij += 1
    kl = 0
    do k=1,mo_tot_num
     do l=1,mo_tot_num
      kl += 1
      mat_i(kl,ij) = two_bod_alpha_beta_mo_physician(l,k,j,i,istate)
     enddo
    enddo
   enddo
  enddo
  !!!!!!!test SVD!!!!!!!
 
   call svd(mat_i,size(mat_i,1),u_i,size(u_i,1),D_i,vt_i,size(vt_i,1),size(mat_i,1),size(mat_i,2))
 
  do i = 1,n_singular_approx_svd_two_dm(istate)
   singular_approx_svd_two_dm(i,istate) = D_i(i)
   do j = 1,mo_tot_num**2
    singular_left_vec_approx_svd_two_dm(j,i,istate) = u_i(j,i)
    singular_right_vec_approx_svd_two_dm(j,i,istate) = vt_i(i,j)
   enddo 
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, singular_left_vec_approx_svd_two_dm_at_r , (n_max_singular_approx_svd_two_dm,n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, singular_right_vec_approx_svd_two_dm_at_r, (n_max_singular_approx_svd_two_dm,n_points_final_grid,N_states)]
 BEGIN_DOC 
! singular_left_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith left singular vector at the jth point of the grid for the state istate
! singular_right_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith right singular vector at the jth point of the grid for the state istate
 END_DOC
 implicit none
 integer :: i,j,k

 singular_left_vec_approx_svd_two_dm_at_r = 0d0
 singular_right_vec_approx_svd_two_dm_at_r = 0d0
 !!!!!!unfoldage!!!!!!!
 integer :: istate
 do istate = 1, N_states
  do i = 1,n_points_final_grid
   do j = 1,n_singular_approx_svd_two_dm(istate) 
    do k = 1,mo_tot_num**2
     singular_left_vec_approx_svd_two_dm_at_r(j,i,istate) +=  singular_left_vec_approx_svd_two_dm(k,j,istate) * mos_in_r_array(approx_two_dm_map(k,1),i) * mos_in_r_array(approx_two_dm_map(k,2),i)
     singular_right_vec_approx_svd_two_dm_at_r(j,i,istate)+= singular_right_vec_approx_svd_two_dm(k,j,istate) * mos_in_r_array(approx_two_dm_map(k,1),i) * mos_in_r_array(approx_two_dm_map(k,2),i) 
    enddo
   enddo
  enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, integral_on_top_of_r_approx_svd, (N_states)]
 implicit none
 BEGIN_DOC 
! Numerical integration of the on top pair density approximated by a svd 
 END_DOC 
 integer :: i,k
 double precision :: weight,wall_1,wall_2,wall_3,wall_4
 
 integral_on_top_of_r_approx_svd = 0.d0
!call cpu_time(wall_3)
 provide singular_left_vec_approx_svd_two_dm_at_r
!call cpu_time(wall_4)

!print*,'cpu time SVD provinding = ',wall_4 - wall_3

 call cpu_time(wall_1)
 integer :: istate
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight=final_weight_functions_at_final_grid_points(i)
   integral_on_top_of_r_approx_svd(istate) += on_top_of_r_approx_svd(i,istate) * weight
  enddo
 enddo
 call cpu_time(wall_2)
 print*,'cpu time SVD approx = ',wall_2 - wall_1
END_PROVIDER

BEGIN_PROVIDER [double precision, on_top_of_r_approx_svd, (n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC 
! on top pair density computed with the SVD approximation 
 END_DOC 
 integer :: i,istate
 double precision :: on_top_of_r_approx_svd_function
 ! initialization for paralellization 
 on_top_of_r_approx_svd(1,1) = on_top_of_r_approx_svd_function(1,1)
 do istate = 1, N_states
  do i= 1, n_points_final_grid
   on_top_of_r_approx_svd(i,istate) = on_top_of_r_approx_svd_function(i,istate)
  enddo
 enddo
END_PROVIDER 

double precision function on_top_of_r_approx_svd_function(i,istate)
 implicit none
 BEGIN_DOC
! evaluation at a given point of the grid of the on top pair density approximated by a svd
 END_DOC
 integer, intent(in) :: i,istate
 integer :: k 
 on_top_of_r_approx_svd_function = 0.d0
 do k =1,n_singular_approx_svd_two_dm(istate)
  on_top_of_r_approx_svd_function += singular_approx_svd_two_dm(k,istate) * singular_left_vec_approx_svd_two_dm_at_r(k,i,istate) * singular_right_vec_approx_svd_two_dm_at_r(k,i,istate) 
 enddo
end


 BEGIN_PROVIDER [integer, n_singular_approx_svd_two_dm_correlation, (N_states)]
&BEGIN_PROVIDER [integer, n_max_singular_approx_svd_two_dm_correlation]
 BEGIN_DOC 
! Number of singular vectors (left and right) kept to represent the two-body tensor of each states
! This number depends on thr_ontop_approx 
 END_DOC
 implicit none
 integer :: i,j,k,l,ij,kl
 double precision :: thresh
 thresh= thr_ontop_approx

 print*,'************************'
 print*,'thresh     approx_two_dm    =',thresh
 print*,'************************'

 double precision, allocatable :: mat_i(:,:)
 allocate(mat_i(mo_tot_num**2,mo_tot_num**2))

 double precision, allocatable :: u_i(:,:),vt_i(:,:),D_i(:)
 allocate(u_i(mo_tot_num**2,mo_tot_num**2),vt_i(mo_tot_num**2,mo_tot_num**2),D_i(mo_tot_num**2))

 !!!!!!unfoldage!!!!!!!
 print*,'*** Approximated two-body tensor '
 integer :: istate
 do istate = 1, N_states
  ij = 0
  kl = 0
  print*,'*** For state',istate
  do i=1,mo_tot_num
   do j=1,mo_tot_num
    ij += 1
    kl = 0
    do k=1,mo_tot_num
     do l=1,mo_tot_num
      kl += 1
      mat_i(kl,ij) = two_bod_alpha_beta_mo_physician(l,k,j,i,istate) - one_body_dm_mo_alpha_for_dft(k,i,istate)* one_body_dm_mo_beta_for_dft(l,j,istate)
     enddo
    enddo
   enddo
  enddo
 !!!!!!!test SVD!!!!!!!

  print*,'***** SVD of the two-body tensor *****'
  double precision :: wall_1,wall_2
  call wall_time(wall_1)
  call svd(mat_i,size(mat_i,1),u_i,size(u_i,1),D_i,vt_i,size(vt_i,1),size(mat_i,1),size(mat_i,2))
  call wall_time(wall_2)
 print*,'*****SVD is done ******' 
 print*,'Wall time for SVD = ',wall_2 - wall_1
  n_singular_approx_svd_two_dm_correlation(istate) =1 
  print*,n_singular_approx_svd_two_dm_correlation(istate),D_i(n_singular_approx_svd_two_dm_correlation(istate))
  do while ( (dabs(D_i(n_singular_approx_svd_two_dm_correlation(istate))) .gt. thresh) .AND. (n_singular_approx_svd_two_dm_correlation(istate) .lt. mo_tot_num**2)  )
   n_singular_approx_svd_two_dm_correlation(istate) += 1
   print*,n_singular_approx_svd_two_dm_correlation(istate),D_i(n_singular_approx_svd_two_dm_correlation(istate))
  enddo
  if (dabs(D_i(n_singular_approx_svd_two_dm_correlation(istate))) .lt. thresh) then
   n_singular_approx_svd_two_dm_correlation(istate) -= 1
  endif
  print*,'************************'
  print*,'n_eigen_approx_two_dm  =',n_singular_approx_svd_two_dm_correlation(istate)
  print*,'************************'
 enddo
 n_max_singular_approx_svd_two_dm_correlation = maxval(n_singular_approx_svd_two_dm_correlation)

END_PROVIDER

 BEGIN_PROVIDER [double precision, singular_approx_svd_two_dm_correlation, (n_max_singular_approx_svd_two_dm_correlation, N_states)]
&BEGIN_PROVIDER [double precision, singular_left_vec_approx_svd_two_dm_correlation, (mo_tot_num**2,n_max_singular_approx_svd_two_dm_correlation,N_states)]
&BEGIN_PROVIDER [double precision, singular_right_vec_approx_svd_two_dm_correlation, (mo_tot_num**2,n_max_singular_approx_svd_two_dm_correlation,N_states)]
 BEGIN_DOC 
 !singular_approx_svd_two_dm(i,istate) = ith Singular values to represent the two-body tensor of the state istate 
 !singular_left_vec_approx_svd_two_dm(I,i,istate) = component on the couple of MOs approx_two_dm_map(I,1) and approx_two_dm_map(I,2) of the ith Singular left vector of the state istate 
 !singular_right_vec_approx_svd_two_dm(I,i,istate) = component on the couple of MOs approx_two_dm_map(I,1) and approx_two_dm_map(I,2) of the ith Singular right vector of the state istate 
 END_DOC
 implicit none
 integer :: i,j,k,l,ij,kl
 double precision, allocatable :: mat_i(:,:)
 allocate(mat_i(mo_tot_num**2,mo_tot_num**2))

 double precision, allocatable :: u_i(:,:),vt_i(:,:),D_i(:)
 allocate(u_i(mo_tot_num**2,mo_tot_num**2),vt_i(mo_tot_num**2,mo_tot_num**2),D_i(mo_tot_num**2))

 integer :: istate
 do istate = 1, N_states
  ij = 0
  kl = 0
  !!!!!!unfoldage!!!!!!!
  do i = 1,mo_tot_num
   do j=1,mo_tot_num
    ij += 1
    kl = 0
    do k=1,mo_tot_num
     do l=1,mo_tot_num
      kl += 1
      mat_i(kl,ij) = two_bod_alpha_beta_mo_physician(l,k,j,i,istate) - one_body_dm_mo_alpha_for_dft(k,i,istate)*one_body_dm_mo_beta_for_dft(l,j,istate)
     enddo
    enddo
   enddo
  enddo
  !!!!!!!test SVD!!!!!!!
 
   call svd(mat_i,size(mat_i,1),u_i,size(u_i,1),D_i,vt_i,size(vt_i,1),size(mat_i,1),size(mat_i,2))
 
  do i = 1,n_singular_approx_svd_two_dm_correlation(istate)
   singular_approx_svd_two_dm_correlation(i,istate) = D_i(i)
   do j = 1,mo_tot_num**2
    singular_left_vec_approx_svd_two_dm_correlation(j,i,istate) = u_i(j,i)
    singular_right_vec_approx_svd_two_dm_correlation(j,i,istate) = vt_i(i,j)
   enddo 
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, singular_left_vec_approx_svd_two_dm_at_r_correlation , (n_max_singular_approx_svd_two_dm_correlation,n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, singular_right_vec_approx_svd_two_dm_at_r_correlation, (n_max_singular_approx_svd_two_dm_correlation,n_points_final_grid,N_states)]
 BEGIN_DOC 
 !singular_left_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith left singular vector at the jth point of the grid for the state istate
 !singular_right_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith right singular vector at the jth point of the grid for the state istate
 END_DOC
 implicit none
 integer :: i,j,k

 singular_left_vec_approx_svd_two_dm_at_r_correlation = 0d0
 singular_right_vec_approx_svd_two_dm_at_r_correlation = 0d0
 !!!!!!unfoldage!!!!!!!
 integer :: istate
 do istate = 1, N_states
  do i = 1,n_points_final_grid
   do j = 1,n_singular_approx_svd_two_dm_correlation(istate) 
    do k = 1,mo_tot_num**2
     singular_left_vec_approx_svd_two_dm_at_r_correlation(j,i,istate) +=  singular_left_vec_approx_svd_two_dm_correlation(k,j,istate) * mos_in_r_array(approx_two_dm_map(k,1),i) * mos_in_r_array(approx_two_dm_map(k,2),i)
     singular_right_vec_approx_svd_two_dm_at_r_correlation(j,i,istate)+= singular_right_vec_approx_svd_two_dm_correlation(k,j,istate) * mos_in_r_array(approx_two_dm_map(k,1),i) * mos_in_r_array(approx_two_dm_map(k,2),i) 
    enddo
   enddo
  enddo
 enddo

END_PROVIDER

 BEGIN_PROVIDER [double precision, integral_on_top_of_r_approx_svd_correlation, (N_states)]
 implicit none
 BEGIN_DOC 
 !Numerical integration of the on top pair density approximated by a svd 
 END_DOC 
 integer :: i,k
 double precision :: weight,wall_1,wall_2,wall_3,wall_4

 provide singular_right_vec_approx_svd_two_dm_at_r_correlation  
 
 call cpu_time(wall_1)
 integer :: istate
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight=final_weight_functions_at_final_grid_points(i)
   integral_on_top_of_r_approx_svd_correlation(istate) += on_top_of_r_approx_svd_correlation(i,istate) * weight
  enddo
 enddo
 call cpu_time(wall_2)
 print*,'cpu time SVD correlation approx = ',wall_2 - wall_1
END_PROVIDER

 BEGIN_PROVIDER [double precision, on_top_of_r_approx_svd_correlation, (n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC 
 !on top pair density computed with the SVD approximation 
 END_DOC 
 integer :: i,istate
 double precision :: on_top_of_r_approx_svd_correlation_function
 double precision :: r(3)
 double precision :: dm_a,dm_b 
 ! initialization for paralellization 
 on_top_of_r_approx_svd_correlation(1,1) = on_top_of_r_approx_svd_correlation_function(1,1)
 do istate = 1, N_states
  do i= 1, n_points_final_grid
   dm_a = one_body_dm_alpha_at_r(i,istate)
   dm_b = one_body_dm_beta_at_r(i,istate)
   on_top_of_r_approx_svd_correlation(i,istate) = dm_a*dm_b + on_top_of_r_approx_svd_correlation_function(i,istate)
  enddo
 enddo
END_PROVIDER 

double precision function on_top_of_r_approx_svd_correlation_function(i,istate)
 implicit none
 BEGIN_DOC
 !evaluation at a given point of the grid of the on top pair density approximated by a svd
 END_DOC
 integer, intent(in) :: i,istate
 integer :: k
 on_top_of_r_approx_svd_correlation_function = 0.d0
 do k =1,n_singular_approx_svd_two_dm_correlation(istate)
  on_top_of_r_approx_svd_correlation_function += singular_approx_svd_two_dm_correlation(k,istate) * singular_left_vec_approx_svd_two_dm_at_r_correlation(k,i,istate) * singular_right_vec_approx_svd_two_dm_at_r_correlation(k,i,istate) 
 enddo
end
