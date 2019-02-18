
 BEGIN_PROVIDER [integer, approx_two_dm_map, (mo_num**2,2)]
&BEGIN_PROVIDER [integer, approx_two_dm_map_rev, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC 
! map a couple of orbitals to an index  
 END_DOC
 integer :: i,j,ij
 ij = 0
 do i = 1,mo_num
  do j=1,mo_num
   ij += 1
   approx_two_dm_map(ij,1) = i
   approx_two_dm_map(ij,2) = j
   approx_two_dm_map_rev(i,j) = ij
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [integer, n_approx_svd_two_dm, (N_states)]
&BEGIN_PROVIDER [integer, n_max_approx_svd_two_dm]
 BEGIN_DOC 
! Number of svd vectors (l and r) kept to represent the two-body tensor of each states
! This number depends on thr_ontop_approx 
 END_DOC
 implicit none
 integer :: i,j,k,l,ij,kl
 double precision :: thresh
 thresh= thr_ontop_approx
 provide two_bod_alpha_beta_mo_physicist 

 print*,'************************'
 print*,'thresh     approx_two_dm    =',thresh
 print*,'************************'

 double precision, allocatable :: mat_i(:,:)
 allocate(mat_i(mo_num**2,mo_num**2))

 double precision, allocatable :: u_i(:,:),vt_i(:,:),D_i(:)
 allocate(u_i(mo_num**2,mo_num**2),vt_i(mo_num**2,mo_num**2),D_i(mo_num**2))

 !!!!!!unfoldage!!!!!!!
 print*,'*** Approximated two-body tensor '
 integer :: istate
 do istate = 1, N_states
  ij = 0
  kl = 0
  print*,'*** For state',istate
  do i=1,mo_num
   do j=1,mo_num
    ij += 1
    kl = 0
    do k=1,mo_num
     do l=1,mo_num
      kl += 1
      mat_i(kl,ij) = two_bod_alpha_beta_mo_physicist(l,k,j,i,istate)
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
  n_approx_svd_two_dm(istate) =1 
  print*,n_approx_svd_two_dm(istate),D_i(n_approx_svd_two_dm(istate))
  do while ( (dabs(D_i(n_approx_svd_two_dm(istate))) .gt. thresh) .AND. (n_approx_svd_two_dm(istate) .lt. mo_num**2)  )
   n_approx_svd_two_dm(istate) += 1
   print*,n_approx_svd_two_dm(istate),D_i(n_approx_svd_two_dm(istate))
  enddo
  if (dabs(D_i(n_approx_svd_two_dm(istate))) .lt. thresh) then
   n_approx_svd_two_dm(istate) -= 1
  endif
  print*,'************************'
  print*,'n_eigen_approx_two_dm  =',n_approx_svd_two_dm(istate)
  print*,'************************'
 enddo
 n_max_approx_svd_two_dm = maxval(n_approx_svd_two_dm)

END_PROVIDER

 BEGIN_PROVIDER [double precision, approx_svd_two_dm, (n_max_approx_svd_two_dm, N_states)]
&BEGIN_PROVIDER [double precision, l_vec_approx_svd_two_dm, (mo_num**2,n_max_approx_svd_two_dm,N_states)]
&BEGIN_PROVIDER [double precision, r_vec_approx_svd_two_dm, (mo_num**2,n_max_approx_svd_two_dm,N_states)]
&BEGIN_PROVIDER [double precision, wall_time_approx_svd_two_dm ]
 BEGIN_DOC 
! approx_svd_two_dm(i,istate) = ith Singular values to represent the two-body tensor of the state istate 
! l_vec_approx_svd_two_dm(I,i,istate) = component on the couple of MOs approx_two_dm_map(I,1) and approx_two_dm_map(I,2) of the ith Singular l vector of the state istate 
! r_vec_approx_svd_two_dm(I,i,istate) = component on the couple of MOs approx_two_dm_map(I,1) and approx_two_dm_map(I,2) of the ith Singular r vector of the state istate 
 END_DOC
 implicit none
 integer :: i,j,k,l,ij,kl
 double precision, allocatable :: mat_i(:,:)
 allocate(mat_i(mo_num**2,mo_num**2))

 double precision, allocatable :: u_i(:,:),vt_i(:,:),D_i(:)
 allocate(u_i(mo_num**2,mo_num**2),vt_i(mo_num**2,mo_num**2),D_i(mo_num**2))

 provide two_bod_alpha_beta_mo_physicist 
 integer :: istate
 double precision :: wall_1,wall_2
 call wall_time(wall_1)
 do istate = 1, N_states
  ij = 0
  kl = 0
  !!!!!!unfoldage!!!!!!!
  do i = 1,mo_num
   do j=1,mo_num
    ij += 1
    kl = 0
    do k=1,mo_num
     do l=1,mo_num
      kl += 1
      mat_i(kl,ij) = two_bod_alpha_beta_mo_physicist(l,k,j,i,istate)
     enddo
    enddo
   enddo
  enddo
  !!!!!!!test SVD!!!!!!!
 
   call svd(mat_i,size(mat_i,1),u_i,size(u_i,1),D_i,vt_i,size(vt_i,1),size(mat_i,1),size(mat_i,2))
 
  do i = 1,n_approx_svd_two_dm(istate)
   approx_svd_two_dm(i,istate) = D_i(i)
   do j = 1,mo_num**2
    l_vec_approx_svd_two_dm(j,i,istate) = u_i(j,i)
    r_vec_approx_svd_two_dm(j,i,istate) = vt_i(i,j)
   enddo 
  enddo
 enddo
  call wall_time(wall_2)
 print*,'*****important SV left and right vectors provided ******' 
 print*,'Wall time for SVD = ',wall_2 - wall_1
 wall_time_approx_svd_two_dm = wall_2 - wall_1

END_PROVIDER


 BEGIN_PROVIDER [double precision, l_vec_approx_svd_two_dm_at_r , (n_max_approx_svd_two_dm,n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, r_vec_approx_svd_two_dm_at_r, (n_max_approx_svd_two_dm,n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, wall_time_r_l_vec]
 BEGIN_DOC 
! l_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith l svd vector at the jth point of the grid for the state istate
! r_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith r svd vector at the jth point of the grid for the state istate
 END_DOC
 implicit none
 integer :: i,j,k

 l_vec_approx_svd_two_dm_at_r = 0d0
 r_vec_approx_svd_two_dm_at_r = 0d0
 !!!!!!unfoldage!!!!!!!
 integer :: istate
 double precision :: wall_1,wall_2
 provide l_vec_approx_svd_two_dm 
 if(mat_mul_svd_vectors)then
  l_vec_approx_svd_two_dm_at_r = l_vec_approx_svd_two_dm_at_r_mat_mul
  r_vec_approx_svd_two_dm_at_r = r_vec_approx_svd_two_dm_at_r_mat_mul
  wall_time_r_l_vec = wall_time_r_l_mat_mul 
 else
  l_vec_approx_svd_two_dm_at_r = l_vec_approx_svd_two_dm_at_r_vec_mul
  r_vec_approx_svd_two_dm_at_r = r_vec_approx_svd_two_dm_at_r_vec_mul
  wall_time_r_l_vec = wall_time_r_l_vec_mul 
 endif
 print*,'time to provide l_vec_approx_svd_two_dm_at_r = ',wall_2 - wall_1

END_PROVIDER

 BEGIN_PROVIDER [double precision, l_vec_approx_svd_two_dm_at_r_mat_mul , (n_max_approx_svd_two_dm,n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, r_vec_approx_svd_two_dm_at_r_mat_mul, (n_max_approx_svd_two_dm,n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, wall_time_r_l_mat_mul ]
 BEGIN_DOC 
! l_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith l svd vector at the jth point of the grid for the state istate
! r_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith r svd vector at the jth point of the grid for the state istate
 END_DOC
 implicit none

 l_vec_approx_svd_two_dm_at_r_mat_mul = 0d0
 r_vec_approx_svd_two_dm_at_r_mat_mul = 0d0

 provide l_vec_approx_svd_two_dm 
 double precision :: wall_1,wall_2
 call wall_time(wall_1)
 
 integer :: istate,n_vec_local
 double precision, allocatable :: array_local(:,:)
 do istate = 1, N_states
  allocate(array_local(mo_num*mo_num,n_approx_svd_two_dm(istate)))
  n_vec_local = n_approx_svd_two_dm(istate)

  array_local(:,:) = l_vec_approx_svd_two_dm(:,:,istate)
  call dgemm('T','N',n_vec_local,n_points_final_grid,mo_num*mo_num,1.d0,         & 
             array_local,size(array_local,1),                                    &
             mo_ij_at_r,size(mo_ij_at_r,1),0.d0,                                 &
             l_vec_approx_svd_two_dm_at_r_mat_mul,size(l_vec_approx_svd_two_dm_at_r_mat_mul,1)) 

  array_local(:,:) = r_vec_approx_svd_two_dm(:,:,istate)
  call dgemm('T','N',n_vec_local,n_points_final_grid,mo_num*mo_num,1.d0,         & 
             array_local,size(array_local,1),                                    &
             mo_ij_at_r,size(mo_ij_at_r,1),0.d0,                                 &
             r_vec_approx_svd_two_dm_at_r_mat_mul,size(r_vec_approx_svd_two_dm_at_r_mat_mul,1)) 

  deallocate(array_local)
 enddo
 call wall_time(wall_2)
 print*,'time to provide l_vec_approx_svd_two_dm_MAT  = ',wall_2 - wall_1
 wall_time_r_l_mat_mul = wall_2 - wall_1 
 free mo_ij_at_r 

END_PROVIDER


 BEGIN_PROVIDER [double precision, l_vec_approx_svd_two_dm_at_r_vec_mul , (n_max_approx_svd_two_dm,n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, r_vec_approx_svd_two_dm_at_r_vec_mul, (n_max_approx_svd_two_dm,n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, wall_time_r_l_vec_mul ]
 BEGIN_DOC 
! l_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith l svd vector at the jth point of the grid for the state istate
! r_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith r svd vector at the jth point of the grid for the state istate
 END_DOC
 implicit none
 integer :: i,j,ij,ig

 l_vec_approx_svd_two_dm_at_r_vec_mul = 0d0
 r_vec_approx_svd_two_dm_at_r_vec_mul = 0d0

 provide l_vec_approx_svd_two_dm 
 double precision :: wall_1,wall_2
 call wall_time(wall_1)
 
 integer :: istate,n_vec_local
 double precision, allocatable :: array_local_l(:,:),mo_ij_local(:), array_local_r(:,:)
 allocate(mo_ij_local(mo_num*mo_num))
 do istate = 1, N_states
  allocate(array_local_l(mo_num*mo_num,n_approx_svd_two_dm(istate)))
  allocate(array_local_r(mo_num*mo_num,n_approx_svd_two_dm(istate)))
  n_vec_local = n_approx_svd_two_dm(istate)

  array_local_l(:,:) = l_vec_approx_svd_two_dm(:,:,istate)
  array_local_r(:,:) = r_vec_approx_svd_two_dm(:,:,istate)

  do ig = 1, n_points_final_grid
   ij = 0
   do i = 1,mo_num
    do j=1,mo_num
     ij += 1
     mo_ij_local(ij) = mos_in_r_array(i,ig) * mos_in_r_array(j,ig)
    enddo
   enddo
   
  call dgemv('T',mo_num*mo_num,n_vec_local,1.d0,                      & 
             array_local_l,size(array_local_l,1),mo_ij_local,1,0.d0,  & 
             l_vec_approx_svd_two_dm_at_r_vec_mul(1,ig,istate),1)

  call dgemv('T',mo_num*mo_num,n_vec_local,1.d0,                      &
             array_local_r,size(array_local_r,1),mo_ij_local,1,0.d0,  &
             r_vec_approx_svd_two_dm_at_r_vec_mul(1,ig,istate),1)

  enddo
  deallocate(array_local_l,array_local_r)
 enddo
 call wall_time(wall_2)
 print*,'time to provide l_vec_approx_svd_two_dm_VEC  = ',wall_2 - wall_1
 wall_time_r_l_vec_mul = wall_2 - wall_1 

END_PROVIDER




BEGIN_PROVIDER [double precision, int_on_top_of_r_approx_svd, (N_states)]
 implicit none
 BEGIN_DOC 
! Numerical integration of the on top pair density approximated by a svd 
 END_DOC 
 integer :: i,k
 double precision :: weight,wall_1,wall_2,wall_3,wall_4
 
 int_on_top_of_r_approx_svd = 0.d0
!call cpu_time(wall_3)
 provide l_vec_approx_svd_two_dm_at_r
!call cpu_time(wall_4)

!print*,'cpu time SVD provinding = ',wall_4 - wall_3

 call wall_time(wall_1)
 integer :: istate
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   int_on_top_of_r_approx_svd(istate) += on_top_of_r_approx_svd(i,istate) * weight
  enddo
 enddo
 call wall_time(wall_2)
 print*,'cpu time SVD approx = ',wall_2 - wall_1
END_PROVIDER

 BEGIN_PROVIDER [double precision, on_top_of_r_approx_svd, (n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, wall_time_on_top_of_r_approx_svd ]
 implicit none
 BEGIN_DOC 
! on top pair density computed with the SVD approximation 
 END_DOC 
 integer :: i,istate
 double precision :: on_top_of_r_approx_svd_function
 double precision :: wall_1,wall_2
 ! initialization for paralellization 
 on_top_of_r_approx_svd(1,1) = on_top_of_r_approx_svd_function(1,1)
 call wall_time(wall_1) 
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i,istate) & 
 !$OMP SHARED(on_top_of_r_approx_svd,n_points_final_grid,N_states)
 do istate = 1, N_states
  do i= 1, n_points_final_grid
   on_top_of_r_approx_svd(i,istate) = on_top_of_r_approx_svd_function(i,istate)
  enddo
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall_2) 
 wall_time_on_top_of_r_approx_svd = wall_2 - wall_1
END_PROVIDER 

double precision function on_top_of_r_approx_svd_function(i,istate)
 implicit none
 BEGIN_DOC
! evaluation at a given point of the grid of the on top pair density approximated by a svd
 END_DOC
 integer, intent(in) :: i,istate
 integer :: k 
 on_top_of_r_approx_svd_function = 0.d0
 do k =1,n_approx_svd_two_dm(istate)
  on_top_of_r_approx_svd_function += approx_svd_two_dm(k,istate) * l_vec_approx_svd_two_dm_at_r(k,i,istate) * r_vec_approx_svd_two_dm_at_r(k,i,istate) 
 enddo
end


