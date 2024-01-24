BEGIN_PROVIDER [double precision, mo_bi_ortho_tc_two_e_chemist, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! mo_bi_ortho_tc_two_e_chemist(k,i,l,j) = <k l|V(r_12)|i j> where i,j are right MOs and k,l are left MOs
 END_DOC
 integer :: i,j,k,l,m,n,p,q
 double precision, allocatable :: mo_tmp_1(:,:,:,:),mo_tmp_2(:,:,:,:),mo_tmp_3(:,:,:,:)
 
 allocate(mo_tmp_1(mo_num,ao_num,ao_num,ao_num))
 mo_tmp_1 = 0.d0
 do m = 1, ao_num
  do p = 1, ao_num
   do n = 1, ao_num
    do q = 1, ao_num
     do k = 1, mo_num
      !       (k n|p m)    = sum_q c_qk * (q n|p m)
      mo_tmp_1(k,n,p,m) += mo_l_coef_transp(k,q) * ao_two_e_tc_tot(q,n,p,m)
     enddo
    enddo
   enddo
  enddo
 enddo
 allocate(mo_tmp_2(mo_num,mo_num,ao_num,ao_num))
 mo_tmp_2 = 0.d0
 do m = 1, ao_num
  do p = 1, ao_num
   do n = 1, ao_num
    do i = 1, mo_num
     do k = 1, mo_num
      !       (k i|p m) = sum_n c_ni * (k n|p m)
      mo_tmp_2(k,i,p,m) += mo_r_coef_transp(i,n) * mo_tmp_1(k,n,p,m)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(mo_tmp_1)
 allocate(mo_tmp_1(mo_num,mo_num,mo_num,ao_num))
 mo_tmp_1 = 0.d0
 do m = 1, ao_num
  do p = 1, ao_num
   do l = 1, mo_num
    do i = 1, mo_num
     do k = 1, mo_num
      mo_tmp_1(k,i,l,m) += mo_l_coef_transp(l,p) * mo_tmp_2(k,i,p,m)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(mo_tmp_2)
 mo_bi_ortho_tc_two_e_chemist = 0.d0 
 do m = 1, ao_num
  do j = 1, mo_num
   do l = 1, mo_num
    do i = 1, mo_num
     do k = 1, mo_num
      mo_bi_ortho_tc_two_e_chemist(k,i,l,j) += mo_r_coef_transp(j,m)  * mo_tmp_1(k,i,l,m)
     enddo
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, mo_bi_ortho_tc_two_e, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! mo_bi_ortho_tc_two_e(k,l,i,j) = <k l| V(r_12) |i j> where i,j are right MOs and k,l are left MOs
!
! the potential V(r_12) contains ALL TWO-E CONTRIBUTION OF THE TC-HAMILTONIAN
 END_DOC
 integer :: i,j,k,l
 do j = 1, mo_num
  do i = 1, mo_num
   do l = 1, mo_num
    do k = 1, mo_num
     !                                                         (k i|l j) = <k l|V(r_12)|i j>
     mo_bi_ortho_tc_two_e(k,l,i,j) = mo_bi_ortho_tc_two_e_chemist(k,i,l,j)
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 
