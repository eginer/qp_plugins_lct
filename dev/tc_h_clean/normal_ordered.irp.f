BEGIN_PROVIDER [ double precision, normal_two_body, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,j,h1,p1,h2,p2
 integer :: kkk,k,m,l,n,hh1,hh2,pp1,pp2,kk
 integer                        :: occ(N_int*bit_kind_size,2),Ne(2)
 integer(bit_kind), allocatable :: key_i_core(:,:)
 double precision :: integral
 allocate(key_i_core(N_int,2))
 if(core_tc_op)then
  do i = 1, N_int
   key_i_core(i,1) = xor(ref_bitmask(i,1),core_bitmask(i,1))
   key_i_core(i,2) = xor(ref_bitmask(i,2),core_bitmask(i,2))
  enddo
  call bitstring_to_list_ab(key_i_core,occ,Ne,N_int)
 else
  call bitstring_to_list_ab(ref_bitmask,occ,Ne,N_int)
 endif
 normal_two_body = 0.d0
 do hh1 = 1, n_act_orb
  h1 = list_act(hh1) 
  do pp1 = hh1+1, n_act_orb
   p1 = list_act(pp1)
   do hh2 = 1, n_act_orb
    h2 = list_act(hh2) 
    do pp2 = hh2+1, n_act_orb
     p2 = list_act(pp2)
      do kkk = 1, Ne(2) 
       ! double alpha / alpha 
       kk = occ(kkk,2)
       k = kk
       j = h1
       m = h2
       l = p1
       n = p2
       ! direct beta - alpha/alpha
       call give_integrals_3_body(j,m,k,l,n,k,integral)
       normal_two_body(p2,h2,p1,h1) += -2.d0 * integral 
       ! three_body_5_index_132
       call give_integrals_3_body(j,m,k,k,l,n,integral)
       normal_two_body(p2,h2,p1,h1) += -integral 
       ! three_body_5_index_312
       call give_integrals_3_body(j,m,k,n,k,l,integral)
       normal_two_body(p2,h2,p1,h1) += -integral 
       ! three_body_5_index_exch_12
       call give_integrals_3_body(j,m,k,n,l,k,integral)
       normal_two_body(p2,h2,p1,h1) += integral 
       ! three_body_5_index_exch_13 
       call give_integrals_3_body(j,m,k,k,l,n,integral)
       normal_two_body(p2,h2,p1,h1) += integral 
       ! three_body_5_index_exch_32
       call give_integrals_3_body(j,m,k,k,n,l,integral)
       normal_two_body(p2,h2,p1,h1) += integral 
      enddo
      do kkk = Ne(2)+1, Ne(1) ! alpha/alpha/alpha
       kk = occ(kkk,1)
       k = kk
       j = h1
       m = h2
       l = p1
       n = p2
       ! direct 
       call give_integrals_3_body(j,m,k,l,n,k,integral)
       normal_two_body(p2,h2,p1,h1) += -2.d0 * integral 
       ! three_body_5_index_132
       call give_integrals_3_body(j,m,k,k,l,n,integral)
       normal_two_body(p2,h2,p1,h1) += -integral 
       ! three_body_5_index_312
       call give_integrals_3_body(j,m,k,n,k,l,integral)
       normal_two_body(p2,h2,p1,h1) += -integral 
       ! three_body_5_index_exch_12
       call give_integrals_3_body(j,m,k,n,l,k,integral)
       normal_two_body(p2,h2,p1,h1) += integral 
       ! three_body_5_index_exch_13 
       call give_integrals_3_body(j,m,k,k,l,n,integral)
       normal_two_body(p2,h2,p1,h1) += integral 
       ! three_body_5_index_exch_32
       call give_integrals_3_body(j,m,k,k,n,l,integral)
       normal_two_body(p2,h2,p1,h1) += integral 
      enddo 
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 
