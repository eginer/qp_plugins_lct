
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

 integer :: istate
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

END_PROVIDER


 BEGIN_PROVIDER [double precision, l_vec_approx_svd_two_dm_at_r , (n_max_approx_svd_two_dm,n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, r_vec_approx_svd_two_dm_at_r, (n_max_approx_svd_two_dm,n_points_final_grid,N_states)]
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
 do istate = 1, N_states
  do i = 1,n_points_final_grid
   do j = 1,n_approx_svd_two_dm(istate) 
    do k = 1,mo_num**2
     l_vec_approx_svd_two_dm_at_r(j,i,istate) +=  l_vec_approx_svd_two_dm(k,j,istate) * mos_in_r_array(approx_two_dm_map(k,1),i) * mos_in_r_array(approx_two_dm_map(k,2),i)
     r_vec_approx_svd_two_dm_at_r(j,i,istate)+= r_vec_approx_svd_two_dm(k,j,istate) * mos_in_r_array(approx_two_dm_map(k,1),i) * mos_in_r_array(approx_two_dm_map(k,2),i) 
    enddo
   enddo
  enddo
 enddo

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

 call cpu_time(wall_1)
 integer :: istate
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   int_on_top_of_r_approx_svd(istate) += on_top_of_r_approx_svd(i,istate) * weight
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
 do k =1,n_approx_svd_two_dm(istate)
  on_top_of_r_approx_svd_function += approx_svd_two_dm(k,istate) * l_vec_approx_svd_two_dm_at_r(k,i,istate) * r_vec_approx_svd_two_dm_at_r(k,i,istate) 
 enddo
end


 BEGIN_PROVIDER [integer, n_approx_svd_two_dm_correl, (N_states)]
&BEGIN_PROVIDER [integer, n_max_approx_svd_two_dm_correl]
 BEGIN_DOC 
! Number of svd vectors (l and r) kept to represent the two-body tensor of each states
! This number depends on thr_ontop_approx 
 END_DOC
 implicit none
 integer :: i,j,k,l,ij,kl
 double precision :: thresh
 thresh= thr_ontop_approx

 print*,'************************'
 print*,'thresh     approx_two_dm correl   =',thresh
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
      mat_i(kl,ij) = two_bod_alpha_beta_mo_physicist(l,k,j,i,istate) - one_e_dm_mo_alpha_for_dft(k,i,istate)* one_e_dm_mo_beta_for_dft(l,j,istate)
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
  n_approx_svd_two_dm_correl(istate) =1 
  print*,n_approx_svd_two_dm_correl(istate),D_i(n_approx_svd_two_dm_correl(istate))
  do while ( (dabs(D_i(n_approx_svd_two_dm_correl(istate))) .gt. thresh) .AND. (n_approx_svd_two_dm_correl(istate) .lt. mo_num**2)  )
   n_approx_svd_two_dm_correl(istate) += 1
   print*,n_approx_svd_two_dm_correl(istate),D_i(n_approx_svd_two_dm_correl(istate))
  enddo
  if (dabs(D_i(n_approx_svd_two_dm_correl(istate))) .lt. thresh) then
   n_approx_svd_two_dm_correl(istate) -= 1
  endif
  print*,'************************'
  print*,'n_eigen_approx_two_dm  =',n_approx_svd_two_dm_correl(istate)
  print*,'************************'
 enddo
 n_max_approx_svd_two_dm_correl = maxval(n_approx_svd_two_dm_correl)

END_PROVIDER

 BEGIN_PROVIDER [double precision, approx_svd_two_dm_correl, (n_max_approx_svd_two_dm_correl, N_states)]
&BEGIN_PROVIDER [double precision, l_vec_approx_svd_two_dm_correl, (mo_num**2,n_max_approx_svd_two_dm_correl,N_states)]
&BEGIN_PROVIDER [double precision, r_vec_approx_svd_two_dm_correl, (mo_num**2,n_max_approx_svd_two_dm_correl,N_states)]
 BEGIN_DOC 
 !approx_svd_two_dm(i,istate) = ith Singular values to represent the two-body tensor of the state istate 
 !l_vec_approx_svd_two_dm(I,i,istate) = component on the couple of MOs approx_two_dm_map(I,1) and approx_two_dm_map(I,2) of the ith Singular l vector of the state istate 
 !r_vec_approx_svd_two_dm(I,i,istate) = component on the couple of MOs approx_two_dm_map(I,1) and approx_two_dm_map(I,2) of the ith Singular r vector of the state istate 
 END_DOC
 implicit none
 integer :: i,j,k,l,ij,kl
 double precision, allocatable :: mat_i(:,:)
 allocate(mat_i(mo_num**2,mo_num**2))

 double precision, allocatable :: u_i(:,:),vt_i(:,:),D_i(:)
 allocate(u_i(mo_num**2,mo_num**2),vt_i(mo_num**2,mo_num**2),D_i(mo_num**2))

 integer :: istate
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
      mat_i(kl,ij) = two_bod_alpha_beta_mo_physicist(l,k,j,i,istate) - one_e_dm_mo_alpha_for_dft(k,i,istate)*one_e_dm_mo_beta_for_dft(l,j,istate)
     enddo
    enddo
   enddo
  enddo
  !!!!!!!test SVD!!!!!!!
 
   call svd(mat_i,size(mat_i,1),u_i,size(u_i,1),D_i,vt_i,size(vt_i,1),size(mat_i,1),size(mat_i,2))
 
  do i = 1,n_approx_svd_two_dm_correl(istate)
   approx_svd_two_dm_correl(i,istate) = D_i(i)
   do j = 1,mo_num**2
    l_vec_approx_svd_two_dm_correl(j,i,istate) = u_i(j,i)
    r_vec_approx_svd_two_dm_correl(j,i,istate) = vt_i(i,j)
   enddo 
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, l_vec_approx_svd_two_dm_at_r_correl , (n_max_approx_svd_two_dm_correl,n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, r_vec_approx_svd_two_dm_at_r_correl, (n_max_approx_svd_two_dm_correl,n_points_final_grid,N_states)]
 BEGIN_DOC 
 !l_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith l svd vector at the jth point of the grid for the state istate
 !r_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith r svd vector at the jth point of the grid for the state istate
 END_DOC
 implicit none
 integer :: i,j,k

 l_vec_approx_svd_two_dm_at_r_correl = 0d0
 r_vec_approx_svd_two_dm_at_r_correl = 0d0
 !!!!!!unfoldage!!!!!!!
 integer :: istate
 do istate = 1, N_states
  do i = 1,n_points_final_grid
   do j = 1,n_approx_svd_two_dm_correl(istate) 
    do k = 1,mo_num**2
     l_vec_approx_svd_two_dm_at_r_correl(j,i,istate) +=  l_vec_approx_svd_two_dm_correl(k,j,istate) * mos_in_r_array(approx_two_dm_map(k,1),i) * mos_in_r_array(approx_two_dm_map(k,2),i)
     r_vec_approx_svd_two_dm_at_r_correl(j,i,istate)+= r_vec_approx_svd_two_dm_correl(k,j,istate) * mos_in_r_array(approx_two_dm_map(k,1),i) * mos_in_r_array(approx_two_dm_map(k,2),i) 
    enddo
   enddo
  enddo
 enddo

END_PROVIDER

 BEGIN_PROVIDER [double precision, int_on_top_of_r_approx_svd_correl, (N_states)]
 implicit none
 BEGIN_DOC 
 !Numerical integration of the on top pair density approximated by a svd 
 END_DOC 
 integer :: i,k
 double precision :: weight,wall_1,wall_2,wall_3,wall_4

 provide r_vec_approx_svd_two_dm_at_r_correl  
 
 call cpu_time(wall_1)
 integer :: istate
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   int_on_top_of_r_approx_svd_correl(istate) += on_top_of_r_approx_svd_correl(i,istate) * weight
  enddo
 enddo
 call cpu_time(wall_2)
 print*,'cpu time SVD correl approx = ',wall_2 - wall_1
END_PROVIDER

 BEGIN_PROVIDER [double precision, on_top_of_r_approx_svd_correl, (n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC 
 !on top pair density computed with the SVD approximation 
 END_DOC 
 integer :: i,istate
 double precision :: on_top_of_r_approx_svd_correl_function
 double precision :: r(3)
 double precision :: dm_a,dm_b 
 ! initialization for paralellization 
 on_top_of_r_approx_svd_correl(1,1) = on_top_of_r_approx_svd_correl_function(1,1)
 do istate = 1, N_states
  do i= 1, n_points_final_grid
   dm_a = one_e_dm_alpha_at_r(i,istate)
   dm_b = one_e_dm_beta_at_r(i,istate)
   on_top_of_r_approx_svd_correl(i,istate) = dm_a*dm_b + on_top_of_r_approx_svd_correl_function(i,istate)
  enddo
 enddo
END_PROVIDER 

double precision function on_top_of_r_approx_svd_correl_function(i,istate)
 implicit none
 BEGIN_DOC
 !evaluation at a given point of the grid of the on top pair density approximated by a svd
 END_DOC
 integer, intent(in) :: i,istate
 integer :: k
 on_top_of_r_approx_svd_correl_function = 0.d0
 do k =1,n_approx_svd_two_dm_correl(istate)
  on_top_of_r_approx_svd_correl_function += approx_svd_two_dm_correl(k,istate) * l_vec_approx_svd_two_dm_at_r_correl(k,i,istate) * r_vec_approx_svd_two_dm_at_r_correl(k,i,istate) 
 enddo
end

double precision function two_dm_HF_S_D_in_r(m,istate)
 implicit none
 BEGIN_DOC
 ! Hartree fock + simple and double two body density evaluated ar two points in real space 
 END_DOC
 integer, intent(in) :: istate
 integer, intent(in) :: m
 integer :: i,j,k,l
 two_dm_HF_S_D_in_r = 0.d0
 !Hartree Fock part
 do j = 1, elec_beta_num 
  do i = 1, elec_alpha_num
  !                                                     1 2 1 2 
   two_dm_HF_S_D_in_r += two_bod_alpha_beta_mo_physicist(i,j,i,j,istate) * mos_in_r_array(i,m) * mos_in_r_array(i,m)*mos_in_r_array(j,m)*mos_in_r_array(j,m)
  enddo
 enddo

 !Simple contributions
 do l = elec_beta_num+1, mo_num
  do i = 1, elec_alpha_num 
   do j = 1, elec_beta_num 
    !                                                     1 2 1 2 
    two_dm_HF_S_D_in_r += two_bod_alpha_beta_mo_physicist(i,j,i,l,istate) * mos_in_r_array(i,m) * mos_in_r_array(i,m)*mos_in_r_array(j,m)*mos_in_r_array(l,m) 
   enddo
  enddo
 enddo

 do l = 1, elec_beta_num 
  do i = 1, elec_alpha_num 
   do j = elec_beta_num+1, mo_num
    !                                                     1 2 1 2 
    two_dm_HF_S_D_in_r += two_bod_alpha_beta_mo_physicist(i,j,i,l,istate) * mos_in_r_array(i,m) * mos_in_r_array(i,m)*mos_in_r_array(j,m)*mos_in_r_array(l,m)

   enddo
  enddo
 enddo

 do j = 1, elec_beta_num 
  do k = elec_alpha_num+1, mo_num
   do i = 1, elec_alpha_num 
    !                                                     1 2 1 2 
    two_dm_HF_S_D_in_r += two_bod_alpha_beta_mo_physicist(i,j,k,j,istate) * mos_in_r_array(i,m) * mos_in_r_array(k,m)*mos_in_r_array(j,m)*mos_in_r_array(j,m) 
   enddo
  enddo
 enddo

 do j = 1, elec_beta_num 
  do k = 1, elec_alpha_num 
   do i = elec_alpha_num+1, mo_num
    !                                                     1 2 1 2 
    two_dm_HF_S_D_in_r += two_bod_alpha_beta_mo_physicist(i,j,k,j,istate) * mos_in_r_array(i,m) * mos_in_r_array(k,m)*mos_in_r_array(j,m)*mos_in_r_array(j,m)
   enddo
  enddo
 enddo
 !double contributions

 do l = elec_beta_num+1, mo_num
  do k = elec_alpha_num+1, mo_num
   do j = 1, elec_beta_num 
    do i = 1, elec_alpha_num
     !                                                     1 2 1 2 
     two_dm_HF_S_D_in_r += two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) * mos_in_r_array(i,m) * mos_in_r_array(k,m)*mos_in_r_array(j,m)*mos_in_r_array(l,m)

    enddo
   enddo
  enddo
 enddo

 do l = 1, elec_beta_num 
  do k = 1, elec_alpha_num
    do j = elec_beta_num+1, mo_num
     do i = elec_alpha_num+1, mo_num
     !                                                     1 2 1 2 
     two_dm_HF_S_D_in_r += two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) * mos_in_r_array(i,m) *mos_in_r_array(k,m)*mos_in_r_array(j,m)*mos_in_r_array(l,m)
    enddo
   enddo
  enddo
 enddo

 double precision :: test
 if(dabs(two_dm_HF_S_D_in_r).lt.1.d-15)then
  two_dm_HF_S_D_in_r = sign(1.d-15,two_dm_HF_S_D_in_r)
 endif

end


BEGIN_PROVIDER [double precision, two_bod_alpha_beta_without_HF_S_D, (mo_num,mo_num,mo_num,mo_num,N_states)]
 implicit none
 BEGIN_DOC
 !  Blablablablablabl
 END_DOC
 integer :: i,j,k,l,istate
 double precision :: cpu_0,cpu_1
 !two_bod_alpha_beta_mo_physicist = 0.d0
 print*,'providing two_bod_alpha_beta_without_HF_S_D ...'
 call cpu_time(cpu_0)
 two_bod_alpha_beta_without_HF_S_D = two_bod_alpha_beta_mo_physicist

 do istate = 1, N_states

  !Hartree Fock part
  do j = 1, elec_beta_num 
   do i = 1, elec_alpha_num
    two_bod_alpha_beta_without_HF_S_D(i,j,i,j,istate) = 0.d0 
   enddo
  enddo
 
  !Simple contributions
  do l = elec_beta_num+1, mo_num
   do i = 1, elec_alpha_num 
    do j = 1, elec_beta_num 
     two_bod_alpha_beta_without_HF_S_D(i,j,i,l,istate) = 0.d0 
    enddo
   enddo
  enddo
 
  do l = 1, elec_beta_num 
   do i = 1, elec_alpha_num 
    do j = elec_beta_num+1, mo_num
     two_bod_alpha_beta_without_HF_S_D(i,j,i,l,istate) = 0.d0 
    enddo
   enddo
  enddo
 
  do j = 1, elec_beta_num 
   do k = elec_alpha_num+1, mo_num
    do i = 1, elec_alpha_num 
     two_bod_alpha_beta_without_HF_S_D(i,j,k,j,istate) = 0.d0 
    enddo
   enddo
  enddo
 
  do j = 1, elec_beta_num 
   do k = 1, elec_alpha_num 
    do i = elec_alpha_num+1, mo_num
     two_bod_alpha_beta_without_HF_S_D(i,j,k,j,istate) = 0.d0 
    enddo
   enddo
  enddo
  !double contributions
 
  do l = elec_beta_num+1, mo_num
   do k = elec_alpha_num+1, mo_num
    do j = 1, elec_beta_num 
     do i = 1, elec_alpha_num
      two_bod_alpha_beta_without_HF_S_D(i,j,k,l,istate) = 0.d0 
     enddo
    enddo
   enddo
  enddo
 
  do l = 1, elec_beta_num 
   do k = 1, elec_alpha_num
    do j = elec_beta_num+1, mo_num
     do i = elec_alpha_num+1, mo_num
      two_bod_alpha_beta_without_HF_S_D(i,j,k,l,istate) = 0.d0 
     enddo
    enddo
   enddo
  enddo

 enddo
 call cpu_time(cpu_1)
 print*,'two_bod_alpha_beta_without_HF_S_D provided in',dabs(cpu_1-cpu_0)

 END_PROVIDER


 BEGIN_PROVIDER [integer, n_approx_svd_two_dm_hf_s_d, (N_states)]
&BEGIN_PROVIDER [integer, n_max_approx_svd_two_dm_hf_s_d]
 BEGIN_DOC 
! Number of svd vectors (l and r) kept to represent the two-body tensor of each states
! This number depends on thr_ontop_approx 
 END_DOC
 implicit none
 integer :: i,j,k,l,ij,kl
 double precision :: thresh
 thresh= thr_ontop_approx

 print*,'************************'
 print*,'thresh     approx_two_dm   HF S and D  =',thresh
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
      mat_i(kl,ij) = two_bod_alpha_beta_without_HF_S_D(i,j,k,l,istate) 
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
  n_approx_svd_two_dm_hf_s_d(istate) =1 
  print*,n_approx_svd_two_dm_hf_s_d(istate),D_i(n_approx_svd_two_dm_hf_s_d(istate))
  do while ( (dabs(D_i(n_approx_svd_two_dm_hf_s_d(istate))) .gt. thresh) .AND. (n_approx_svd_two_dm_hf_s_d(istate) .lt. mo_num**2)  )
   n_approx_svd_two_dm_hf_s_d(istate) += 1
   print*,n_approx_svd_two_dm_hf_s_d(istate),D_i(n_approx_svd_two_dm_hf_s_d(istate))
  enddo
  if (dabs(D_i(n_approx_svd_two_dm_hf_s_d(istate))) .lt. thresh) then
   n_approx_svd_two_dm_hf_s_d(istate) -= 1
  endif
  print*,'************************'
  print*,'n_eigen_approx_two_dm_hf_s_d  =',n_approx_svd_two_dm_hf_s_d(istate)
  print*,'************************'
 enddo
 n_max_approx_svd_two_dm_hf_s_d = maxval(n_approx_svd_two_dm_hf_s_d)

END_PROVIDER

 BEGIN_PROVIDER [double precision, approx_svd_two_dm_hf_s_d, (n_max_approx_svd_two_dm_hf_s_d, N_states)]
&BEGIN_PROVIDER [double precision, l_vec_approx_svd_two_dm_hf_s_d, (mo_num**2,n_max_approx_svd_two_dm_hf_s_d,N_states)]
&BEGIN_PROVIDER [double precision, r_vec_approx_svd_two_dm_hf_s_d, (mo_num**2,n_max_approx_svd_two_dm_hf_s_d,N_states)]
 BEGIN_DOC 
 !approx_svd_two_dm(i,istate) = ith Singular values to represent the two-body tensor of the state istate 
 !l_vec_approx_svd_two_dm(I,i,istate) = component on the couple of MOs approx_two_dm_map(I,1) and approx_two_dm_map(I,2) of the ith Singular l vector of the state istate 
 !r_vec_approx_svd_two_dm(I,i,istate) = component on the couple of MOs approx_two_dm_map(I,1) and approx_two_dm_map(I,2) of the ith Singular r vector of the state istate 
 END_DOC
 implicit none
 integer :: i,j,k,l,ij,kl
 double precision, allocatable :: mat_i(:,:)
 allocate(mat_i(mo_num**2,mo_num**2))

 double precision, allocatable :: u_i(:,:),vt_i(:,:),D_i(:)
 allocate(u_i(mo_num**2,mo_num**2),vt_i(mo_num**2,mo_num**2),D_i(mo_num**2))

 integer :: istate
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
      mat_i(kl,ij) = two_bod_alpha_beta_without_HF_S_D(i,j,k,l,istate) 
     enddo
    enddo
   enddo
  enddo
  !!!!!!!test SVD!!!!!!!
 
   call svd(mat_i,size(mat_i,1),u_i,size(u_i,1),D_i,vt_i,size(vt_i,1),size(mat_i,1),size(mat_i,2))
 
  do i = 1,n_approx_svd_two_dm_hf_s_d(istate)
   approx_svd_two_dm_hf_s_d(i,istate) = D_i(i)
   do j = 1,mo_num**2
    l_vec_approx_svd_two_dm_hf_s_d(j,i,istate) = u_i(j,i)
    r_vec_approx_svd_two_dm_hf_s_d(j,i,istate) = vt_i(i,j)
   enddo 
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, l_vec_approx_svd_two_dm_at_r_hf_s_d , (n_max_approx_svd_two_dm_hf_s_d,n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, r_vec_approx_svd_two_dm_at_r_hf_s_d, (n_max_approx_svd_two_dm_hf_s_d,n_points_final_grid,N_states)]
 BEGIN_DOC 
 !l_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith l svd vector at the jth point of the grid for the state istate
 !r_vec_approx_svd_two_dm_at_r(i,j,istate) = value in real space of the ith r svd vector at the jth point of the grid for the state istate
 END_DOC
 implicit none
 integer :: i,j,k

 l_vec_approx_svd_two_dm_at_r_hf_s_d = 0d0
 r_vec_approx_svd_two_dm_at_r_hf_s_d = 0d0
 !!!!!!unfoldage!!!!!!!
 integer :: istate
 do istate = 1, N_states
  do i = 1,n_points_final_grid
   do j = 1,n_approx_svd_two_dm_hf_s_d(istate) 
    do k = 1,mo_num**2
     l_vec_approx_svd_two_dm_at_r_hf_s_d(j,i,istate) +=  l_vec_approx_svd_two_dm_hf_s_d(k,j,istate) * mos_in_r_array(approx_two_dm_map(k,1),i) * mos_in_r_array(approx_two_dm_map(k,2),i)
     r_vec_approx_svd_two_dm_at_r_hf_s_d(j,i,istate)+= r_vec_approx_svd_two_dm_hf_s_d(k,j,istate) * mos_in_r_array(approx_two_dm_map(k,1),i) * mos_in_r_array(approx_two_dm_map(k,2),i) 
    enddo
   enddo
  enddo
 enddo

END_PROVIDER

 BEGIN_PROVIDER [double precision, int_on_top_of_r_approx_svd_hf_s_d, (N_states)]
 implicit none
 BEGIN_DOC 
 !Numerical integration of the on top pair density approximated by a svd 
 END_DOC 
 integer :: i,k
 double precision :: weight,wall_1,wall_2,wall_3,wall_4

 provide r_vec_approx_svd_two_dm_at_r_hf_s_d  
 
 call cpu_time(wall_1)
 integer :: istate
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   int_on_top_of_r_approx_svd_hf_s_d(istate) += on_top_of_r_approx_svd_hf_s_d(i,istate) * weight
  enddo
 enddo
 call cpu_time(wall_2)
 print*,'cpu time SVD hf_s_d approx = ',wall_2 - wall_1
END_PROVIDER

 BEGIN_PROVIDER [double precision, on_top_of_r_approx_svd_hf_s_d, (n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC 
 !on top pair density computed with the SVD approximation 
 END_DOC 
 integer :: i,istate
 double precision :: on_top_of_r_approx_svd_hf_s_d_function
 double precision :: two_dm_HF_S_D_in_r
 ! initialization for paralellization 
 do istate = 1, N_states
  do i= 1, n_points_final_grid
   on_top_of_r_approx_svd_hf_s_d(i,istate) = two_dm_HF_S_D_in_r(i,istate)  + on_top_of_r_approx_svd_hf_s_d_function(i,istate)
  enddo
 enddo
END_PROVIDER 

double precision function on_top_of_r_approx_svd_hf_s_d_function(i,istate)
 implicit none
 BEGIN_DOC
 !evaluation at a given point of the grid of the on top pair density approximated by a svd
 END_DOC
 integer, intent(in) :: i,istate
 integer :: k
 on_top_of_r_approx_svd_hf_s_d_function = 0.d0
 do k =1,n_approx_svd_two_dm_hf_s_d(istate)
  on_top_of_r_approx_svd_hf_s_d_function += approx_svd_two_dm_hf_s_d(k,istate) * l_vec_approx_svd_two_dm_at_r_hf_s_d(k,i,istate) * r_vec_approx_svd_two_dm_at_r_hf_s_d(k,i,istate) 
 enddo
end


  BEGIN_PROVIDER [double precision, E_cor_tot_normal_prov]
 implicit none
 integer :: i,m,n,p,q
 double precision :: weight,wall_4,wall_3
 E_cor_tot_normal_prov = 0d0
 provide two_bod_alpha_beta_mo_physicist 

 call cpu_time(wall_3)
 do i = 1, n_points_final_grid
  weight=final_weight_at_r_vector(i)

  do m = 1,mo_num
   do n = 1,mo_num
    do p = 1,mo_num
     do q =1,mo_num
      E_cor_tot_normal_prov += two_bod_alpha_beta_mo_physicist(q,p,n,m,1) * mos_in_r_array(m,i) * mos_in_r_array(n,i) * mos_in_r_array(p,i) * mos_in_r_array(q,i) * weight
     enddo
    enddo
   enddo
  enddo

 enddo
 call cpu_time(wall_4)
print*,'cpu time exact  = ',wall_4 - wall_3
END_PROVIDER

