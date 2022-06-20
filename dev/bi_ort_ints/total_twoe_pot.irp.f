BEGIN_PROVIDER [double precision, ao_two_e_tc_tot, (ao_num, ao_num, ao_num, ao_num) ]
 integer :: i,j,k,l
 BEGIN_DOC
! ao_two_e_tc_tot(k,i,l,j) = (ki|V^TC(r_12)|lj) = <lk| V^TC(r_12) |ji> where V^TC(r_12) is the total TC operator 
!
! including both hermitian and non hermitian parts
!
! WARNING :: non hermitian ! acts on "the right functions" (i,j)
 END_DOC
 double precision :: integral_sym, integral_nsym, get_ao_tc_sym_two_e_pot
 PROVIDE ao_tc_sym_two_e_pot_in_map
 do j = 1, ao_num
  do l = 1, ao_num
   do i = 1, ao_num
    do k = 1, ao_num
     integral_sym  = get_ao_tc_sym_two_e_pot(i,j,k,l,ao_tc_sym_two_e_pot_map)
     ! ao_non_hermit_term_chemist(k,i,l,j) = < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the AO basis
     integral_nsym = ao_non_hermit_term_chemist(k,i,l,j)
     ao_two_e_tc_tot(k,i,l,j) = integral_sym + integral_nsym 
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 


double precision function bi_ortho_mo_ints(l,k,j,i)
 implicit none
 BEGIN_DOC
! <mo^L_k mo^L_l | V^TC(r_12) | mo^R_i mo^R_j>
 END_DOC
 integer, intent(in) :: i,j,k,l
 integer :: m,n,p,q
 bi_ortho_mo_ints = 0.d0
 do m = 1, ao_num
  do p = 1, ao_num
   do n = 1, ao_num
    do q = 1, ao_num
     bi_ortho_mo_ints += ao_two_e_tc_tot(n,q,m,p) * mo_l_coef(m,l) * mo_l_coef(n,k) * mo_r_coef(p,j) * mo_r_coef(q,i)
    enddo
   enddo
  enddo
 enddo

end
