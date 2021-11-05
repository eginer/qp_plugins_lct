
 BEGIN_PROVIDER [ double precision, direct_normal_ordered, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, exch_normal_ordered, (mo_num, mo_num)]
 implicit none
 integer :: i,j
 do i =1, mo_num
  do j = 1, mo_num
   direct_normal_ordered(j,i) = normal_two_body(j,j,i,i)
   exch_normal_ordered(j,i) = normal_two_body(j,i,j,i)
  enddo 
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, three_indx_direct_normal_ordered, (mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_indx_exch_normal_ordered, (mo_num, mo_num, mo_num)]
 implicit none
 integer :: i,j,k
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    three_indx_direct_normal_ordered(k,j,i) = normal_two_body(k,i,k,j)
    three_indx_exch_normal_ordered(k,j,i) = normal_two_body(k,i,j,k)
   enddo
  enddo
 enddo
END_PROVIDER 

