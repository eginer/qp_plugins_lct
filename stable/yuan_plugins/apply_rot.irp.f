BEGIN_PROVIDER [double precision, new_mo_coefs_rot_mat, (ao_num, mo_num)]
 implicit none
 integer :: i,j,k,l
 double precision :: rotation_matrix(mo_num,mo_num)
 ! Read the rotation matrix 
 rotation_matrix = 0.d0
 open(1, file = 'rotation_matrix') 
 do i = 1, mo_num
  read(1,*)rotation_matrix(i,1:mo_num)
  write(*,'(100(F16.10,X))')rotation_matrix(i,1:mo_num)
 enddo
 close(1) 
 ! <AO_k| new_mo_j> = \sum_i U_ij <AO_k| old_mo_i>
 new_mo_coefs_rot_mat = 0.d0
 do j = 1, mo_num ! 
  do i = 1, mo_num
   do k = 1, ao_num
!    do l = 1, ao_num
!     new_mo_coefs_rot_mat(k,j) += mo_coef(k,i) * rotation_matrix(i,j) * ao_overlap(l,k)
     new_mo_coefs_rot_mat(k,j) += mo_coef(k,i) * rotation_matrix(i,j) 
!    enddo
   enddo
  enddo
 enddo
END_PROVIDER 


subroutine save_new_mo_coefs
 implicit none
 character*(64)                 :: label
 mo_coef = new_mo_coefs_rot_mat
 label = "Canonical"
 mo_label = label 
 soft_touch mo_coef mo_label 
 call save_mos
end
