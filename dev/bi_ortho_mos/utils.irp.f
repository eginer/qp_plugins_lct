
subroutine rotate_mo_coef(mo_coef_old,mo_mo_coef,mo_mo_overlap,mo_coef_new)
 implicit none
 BEGIN_DOC
 ! You have mo_coef_new which is based on a MO->MO transformation through mo_mo_coef
 END_DOC
 double precision, intent(in)  :: mo_coef_old(ao_num, mo_num)
 double precision, intent(in)  :: mo_mo_coef(mo_num, mo_num),mo_mo_overlap(mo_num, mo_num)
 double precision, intent(out) :: mo_coef_new(ao_num, mo_num)
! call dgemm('N','N',ao_num,mo_num,mo_num,1.d0,mo_coef_old,size(mo_coef_old,1),& 
!      mo_mo_coef,size(mo_mo_coef,1),&
!      0.d0,mo_coef_new,size(mo_coef_new,1))
 integer :: i,j,q,k
 mo_coef_new = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do q = 1, ao_num
     mo_coef_new(q,i) += mo_coef_old(q,j) * mo_mo_overlap(k,j) * mo_mo_coef(j,i) 
    enddo
   enddo
  enddo
 enddo

end
