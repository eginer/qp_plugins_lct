BEGIN_PROVIDER [ double precision, ao_prod_on_grid, (n_points_final_grid,ao_num,ao_num)]
 implicit none
 integer :: ipoint, i,j,k,l,m
 double precision :: weight
 do i = 1, ao_num
  do j = 1, ao_num
   do ipoint = 1, n_points_final_grid
    weight = final_weight_at_r_vector(ipoint)
    ao_prod_on_grid(ipoint,j,i) = aos_in_r_array_transp(ipoint,j) * aos_in_r_array_transp(ipoint,i) * weight
   enddo
  enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [ double precision, ao_prod_on_grid_transp, (ao_num,ao_num,n_points_final_grid)]
 implicit none
 integer :: ipoint, i,j,k,l,m
 double precision :: weight
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do i = 1, ao_num
   do j = 1, ao_num
    ao_prod_on_grid_transp(j,i,ipoint) = aos_in_r_array(j,ipoint) * aos_in_r_array(i,ipoint) * weight
   enddo
  enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [ double precision, big_array_naive, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 integer :: ipoint, i,j,k,l,m
 big_array_naive = 0.d0
 do ipoint = 1, n_points_final_grid
  do l = 1, ao_num
   do j = 1, ao_num
    do k = 1, ao_num
     do i = 1, ao_num
      big_array_naive(i,k,j,l) += ao_prod_on_grid_transp(i,k,ipoint) * erf_mu_r12_inv_r12_rk(j,l,ipoint)
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, big_array_dgemm, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 big_array_dgemm= 0.d0

 call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,erf_mu_r12_inv_r12_rk(1,1,1),ao_num*ao_num & 
                   ,ao_prod_on_grid(1,1,1),n_points_final_grid,1.d0,big_array_dgemm,ao_num*ao_num)

! call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,erf_mu_r12_inv_r12_rk(1,1,1),ao_num*ao_num & 
!                     ,ao_prod_on_grid(1,1,1,m),n_points_final_grid,1.d0,ac_mat,ao_num*ao_num)
END_PROVIDER 


subroutine test_big_array
 implicit none
 integer :: i,j,k,l
 double precision :: accu,contrib
 accu = 0.d0
  do l = 1, ao_num
   do j = 1, ao_num
    do k = 1, ao_num
     do i = 1, ao_num
      contrib = dabs(big_array_dgemm(i,k,j,l) - big_array_naive(i,k,j,l))
      if(contrib .gt. 1.d-10)then
       print*,'i,k,j,l',i,k,j,l
       print*,contrib,big_array_dgemm(i,k,j,l) , big_array_naive(i,k,j,l)
      endif
      accu += contrib
     enddo
    enddo
   enddo
  enddo
  print*,'accu = ',accu
  print*,'accu/n**4 = ',accu/dble(ao_num**4)

end
