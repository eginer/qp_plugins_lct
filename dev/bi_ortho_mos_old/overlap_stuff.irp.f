 BEGIN_PROVIDER [ double precision, overlap_mo_r, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_mo_l, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! overlap_mo_r_mo(j,i) = <MO_i|MO_R_j>
 END_DOC
 integer :: i,j,p,q
 overlap_mo_r= 0.d0
 overlap_mo_l= 0.d0
 do i = 1, mo_num
   do j = 1, mo_num
    do p = 1, ao_num
     do q = 1, ao_num
      overlap_mo_r(j,i) += mo_r_coef(q,i) * mo_r_coef(p,j) * ao_overlap(q,p) 
      overlap_mo_l(j,i) += mo_l_coef(q,i) * mo_l_coef(p,j) * ao_overlap(q,p)
     enddo
    enddo
   enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, overlap_mo_r_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_mo_l_mo, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! overlap_mo_r_mo(j,i) = <MO_j|MO_R_i>
 END_DOC
 integer :: i,j,p,q
 overlap_mo_r_mo = 0.d0
 overlap_mo_l_mo = 0.d0
 do i = 1, mo_num
   do j = 1, mo_num
    do p = 1, ao_num
     do q = 1, ao_num
      overlap_mo_r_mo(j,i) += mo_coef(p,j) * mo_r_coef(q,i) * ao_overlap(q,p)
      overlap_mo_l_mo(j,i) += mo_coef(p,j) * mo_l_coef(q,i) * ao_overlap(q,p)
     enddo
    enddo
   enddo
 enddo
END_PROVIDER 
