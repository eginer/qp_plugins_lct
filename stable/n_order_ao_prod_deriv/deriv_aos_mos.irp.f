 BEGIN_PROVIDER[double precision, aos_nabla_2_in_r_array, (ao_num,ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
!    aos_nabla_2_in_r_array(i,j,ipoint) = (d/d_x)^2 (ao_i(ipoint),ao_j(ipoint)) 
!        
!                                       + (d/d_y)^2 (ao_i(ipoint),ao_j(ipoint)) 
!                                       
!                                       + (d/d_z)^2 (ao_i(ipoint),ao_j(ipoint)) 
!
! where ipoint corresponds to a point grid
 END_DOC
 integer :: i,j,k
 double precision :: r(3)
 double precision :: nabla_2_at_r(3,ao_num,ao_num)
 double precision :: nabla_2_tot_at_r(ao_num,ao_num)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_nabla_2_at_r(r,nabla_2_at_r,nabla_2_tot_at_r)
  do j = 1, ao_num
   do k=1,ao_num
    aos_nabla_2_in_r_array(k,j,i) = nabla_2_tot_at_r(k,j)
   enddo
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, aos_nabla_4_in_r_array, (ao_num,ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
!    aos_nabla_4_in_r_array(i,j,ipoint) = (d/d_x)^4 (ao_i(ipoint),ao_j(ipoint)) 
!        
!                                       + (d/d_y)^4 (ao_i(ipoint),ao_j(ipoint)) 
!                                       
!                                       + (d/d_z)^4 (ao_i(ipoint),ao_j(ipoint)) 
!
! where ipoint corresponds to a point grid
 END_DOC
 integer :: i,j,k
 double precision :: r(3)
 double precision :: nabla_4_at_r(ao_num,ao_num)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_nabla_4_at_r(r,nabla_4_at_r) 
  do j = 1, ao_num
   do k=1,ao_num
    aos_nabla_4_in_r_array(k,j,i) = nabla_4_at_r(k,j)
   enddo
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_nabla_4_in_r_array_slow, (mo_num,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k,m,n
 double precision :: r(3)
 double precision :: nabla_4_at_r(ao_num,ao_num)
 mos_nabla_4_in_r_array_slow =0.d0
 do i = 1, n_points_final_grid
  do j = 1, mo_num 
   do k = 1, mo_num
    do m = 1, ao_num
     do n = 1, ao_num 
      mos_nabla_4_in_r_array_slow(k,j,i) += mo_coef(m,j)*mo_coef(n,k)* aos_nabla_4_in_r_array(n,m,i)  
     enddo
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER[double precision, mos_nabla_4_in_r_array, (mo_num,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k,m,n
 double precision :: tempo_matrix(mo_num,ao_num)
 do i = 1, n_points_final_grid
  call ao_to_mo(aos_nabla_4_in_r_array(1,1,i),ao_num,mos_nabla_4_in_r_array(1,1,i),mo_num)
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_nabla_2_in_r_array_2, (mo_num,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k,m,n
 double precision :: tempo_matrix(mo_num,ao_num)
 do i = 1, n_points_final_grid
  call ao_to_mo(aos_nabla_2_in_r_array(1,1,i),ao_num,mos_nabla_2_in_r_array_2(1,1,i),mo_num)
 enddo
 END_PROVIDER

