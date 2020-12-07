BEGIN_PROVIDER [double precision, mo_ten_no_dr12_pot_chemist, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
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
      mo_tmp_1(k,n,p,m) += mo_coef_transp(k,q) * ao_ten_no_dr12_pot(q,n,p,m)
     enddo
    enddo
   enddo
  enddo
 enddo
 free ao_ten_no_dr12_pot 
 allocate(mo_tmp_2(mo_num,mo_num,ao_num,ao_num))
 mo_tmp_2 = 0.d0
 do m = 1, ao_num
  do p = 1, ao_num
   do n = 1, ao_num
    do i = 1, mo_num
     do k = 1, mo_num
      !       (k i|p m) = sum_n c_ni * (k n|p m)
      mo_tmp_2(k,i,p,m) += mo_coef_transp(i,n) * mo_tmp_1(k,n,p,m)
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
      mo_tmp_1(k,i,l,m) += mo_coef_transp(l,p) * mo_tmp_2(k,i,p,m)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(mo_tmp_2)
 mo_ten_no_dr12_pot_chemist = 0.d0
 do m = 1, ao_num
  do j = 1, mo_num
   do l = 1, mo_num
    do i = 1, mo_num
     do k = 1, mo_num
      mo_ten_no_dr12_pot_chemist(k,i,l,j) += mo_coef_transp(j,m)  * mo_tmp_1(k,i,l,m)
     enddo
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, mo_ten_no_dr12_pot_physicist, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 integer :: i,j,k,l
 do j = 1, mo_num
  do i = 1, mo_num
   do l = 1, mo_num
    do k = 1, mo_num
     mo_ten_no_dr12_pot_physicist(k,l,i,j) = mo_ten_no_dr12_pot_chemist(k,i,l,j)
    enddo
   enddo
  enddo
 enddo
 free mo_ten_no_dr12_pot_chemist
END_PROVIDER 
