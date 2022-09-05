BEGIN_PROVIDER [ double precision, mo_bi_ortho_tc_one_e, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! mo_bi_ortho_tc_one_e(k,i) = <MO^L_k | h_c | MO^R_i>
 END_DOC
 integer :: i,k,p,q
 
!double precision, allocatable :: mo_tmp(:,:),tmp(:,:)
! allocate(mo_tmp(ao_num, mo_num),tmp(mo_num, mo_num))
! mo_tmp = 0.d0
! do i = 1, mo_num
!  do p = 1, ao_num
!   do q = 1, ao_num
!    mo_tmp(p,i) +=  ao_one_e_integrals(p,q) * mo_r_coef(q,i) 
!   enddo
!  enddo
! enddo
! tmp = 0.d0
! do i = 1, mo_num
!  do k = 1, mo_num
!   do q = 1, ao_num
!    tmp(k,i) +=  mo_tmp(q,i) * mo_l_coef(q,k) 
!   enddo
!  enddo
! enddo
 call ao_to_mo_bi_ortho(ao_one_e_integrals,ao_num,mo_bi_ortho_tc_one_e,mo_num)
!do i = 1, mo_num
! do k = 1, mo_num
!   print*,k,i
!   print*,dabs(mo_bi_ortho_tc_one_e(k,i) - tmp(k,i)),mo_bi_ortho_tc_one_e(k,i),tmp(k,i)
!  if(dabs(mo_bi_ortho_tc_one_e(k,i) - tmp(k,i)).gt.1.d-12)then
!   print*,'problem in mo_bi_ortho_tc_one_e'
!   stop
!  endif
! enddo
!enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, mo_bi_ortho_tc_one_e_slow, (mo_num, mo_num)]
 implicit none
 integer :: i,k,p,q
 mo_bi_ortho_tc_one_e_slow = 0.D0
 do i = 1, mo_num
  do k = 1, mo_num
   do p = 1, ao_num
    do q = 1, ao_num
     mo_bi_ortho_tc_one_e_slow(k,i) +=   mo_l_coef(p,k) * ao_one_e_integrals(q,p) * mo_r_coef(q,i)
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 
