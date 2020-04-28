program coll_salv_fun
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  call test_int 
!  call test_fit
!  call print_mo
end

subroutine print_mo
 implicit none
 integer :: i,j
 double precision :: norm,alpha
 norm = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   norm += mo_coef(i,1) * mo_coef(j,1) * ao_overlap(j,i)
  enddo
 enddo
 print*,'norm = ',norm
 do i = 1, ao_num
  print*,' coef_fit_slat_gauss(',i,') = ',mo_coef(i,1)*ao_coef_normalized_ordered_transp(1,i)
 enddo
 do i = 1, ao_num
  do j=1,ao_prim_num(i)
   alpha = ao_expo_ordered_transp(j,i)     
   print*,' expo_fit_slat_gauss(',i,') = ',alpha
  enddo
 enddo
     

end

