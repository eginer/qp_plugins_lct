BEGIN_PROVIDER [integer, order_d, (3,3)]
 implicit none
 BEGIN_DOC
! needed for printing the AOs in the correct order 
 END_DOC
 order_d = -10000
 order_d(1,1) = 1 ! XX = first
 order_d(2,2) = 2 ! YY = second 
 order_d(3,3) = 3 ! ZZ = third 
 order_d(1,2) = 4 ! XY = fourth
 order_d(1,3) = 5 ! XZ = fifth
 order_d(2,3) = 6 ! YZ = sixth 

END_PROVIDER 

BEGIN_PROVIDER [integer, order_f, (10,10,10)]
 implicit none
 BEGIN_DOC
! needed for printing the AOs in the correct order 
 END_DOC
 order_f = -10000
 order_f(1,1,1) = 1 ! XXX = first
 order_f(2,2,2) = 2 ! YYY = second 
 order_f(3,3,3) = 3 ! ZZZ = third 
 order_f(1,1,2) = 4
 order_f(1,1,3) = 5
 order_f(2,2,1) = 6
 order_f(2,2,3) = 7
 order_f(3,3,1) = 8
 order_f(3,3,2) = 9
 order_f(1,2,3) = 10

END_PROVIDER 

 BEGIN_PROVIDER [integer, list_first_d, (ao_num)]
&BEGIN_PROVIDER [integer, n_first_d]
 implicit none
 integer :: i
 logical :: first
 n_first_d = 0
 list_first_d = -10000
 do i=1,ao_num
  if(ao_l(i) == 2 .and. ao_l_char_space(i) == 'XX  ')then
   n_first_d += 1
   list_first_d(n_first_d) = i
  endif
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [integer, new_order_d_AOS_for_molden, (6)]
&BEGIN_PROVIDER [integer, new_order_d_AOS_reverse_for_molden, (6)]
 implicit none
 integer :: i_ao,k,i
 new_order_d_AOS_for_molden = -1000
 i_ao = list_first_d(1) ! index of the first d
 do k = 1, 6 
  i = i_ao + k-1
  new_order_d_AOS_for_molden(k) = order_d( ao_l_powers(i,1) , ao_l_powers(i,2) )
  new_order_d_AOS_reverse_for_molden(order_d( ao_l_powers(i,1) , ao_l_powers(i,2) )) = k 
 enddo
 
END_PROVIDER 


 BEGIN_PROVIDER [integer, list_first_f, (ao_num)]
&BEGIN_PROVIDER [integer, n_first_f]
 implicit none
 integer :: i
 logical :: first
 n_first_f = 0
 list_first_f = -10000
 do i=1,ao_num
  if(ao_l(i) == 3 .and. ao_l_char_space(i) == 'XXX ')then
   n_first_f += 1
   list_first_f(n_first_f) = i
  endif
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [integer, new_order_f_AOS_for_molden, (10)]
&BEGIN_PROVIDER [integer, new_order_f_AOS_reverse_for_molden, (10)]
 implicit none
 integer :: i_ao,k,i
 new_order_f_AOS_for_molden = -1000
 i_ao = list_first_f(1) ! index of the first f
 do k = 1, 10
  i = i_ao + k-1
  new_order_f_AOS_for_molden(k) = order_f( ao_l_powers(i,1) , ao_l_powers(i,2) ,ao_l_powers(i,3) )
  new_order_f_AOS_reverse_for_molden(order_f( ao_l_powers(i,1) , ao_l_powers(i,2) ,ao_l_powers(i,3) )) = k
 enddo
 
END_PROVIDER 


 BEGIN_PROVIDER [integer, new_order_all_AOS_for_molden, (ao_num)]
 implicit none
 integer :: i,i_new,i_begin,k
 do i = 1, ao_num
  new_order_all_AOS_for_molden(i) = i
 enddo
 
 do i = 1, n_first_d
  i_begin = list_first_d(i)
  do k = 1, 6
   i_new = i_begin + new_order_d_AOS_reverse_for_molden(k) - 1
   new_order_all_AOS_for_molden(i_begin + k-1) = i_new
  enddo
 enddo

 do i = 1, n_first_f
  i_begin = list_first_f(i)
  do k = 1, 10
   i_new = i_begin + new_order_f_AOS_reverse_for_molden(k) - 1
   new_order_all_AOS_for_molden(i_begin + k-1) = i_new
  enddo
 enddo



!do i = 1, ao_num
! print*,'      '
! print*,i,ao_l(i),ao_l_char_space(i),new_order_all_AOS_for_molden(i)
! print*,'after'
! print*,i,ao_l(new_order_all_AOS_for_molden(i)),ao_l_char_space(new_order_all_AOS_for_molden(i))
!enddo



 END_PROVIDER 
