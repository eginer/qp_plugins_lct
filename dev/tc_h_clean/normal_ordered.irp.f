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


BEGIN_PROVIDER [ double precision, normal_two_body_ab, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! Normal ordered two-body sector of the three-body terms for opposite spin double excitations 
 END_DOC
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,ii,j,h1,p1,h2,p2
 integer :: kkk,k,m,l,n,hh1,hh2,pp1,pp2,kk
 integer                        :: occ(N_int*bit_kind_size,2),Ne(2)
 integer(bit_kind), allocatable :: key_i_core(:,:)
 double precision :: int_direct,int_exc_12,int_exc_13,integral
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
 normal_two_body_ab = 0.d0
 do hh1 = 1, n_act_orb
  h1 = list_act(hh1) 
  do pp1 = hh1, n_act_orb
   p1 = list_act(pp1)
   do hh2 = 1, n_act_orb
    h2 = list_act(hh2) 
    do pp2 = hh2, n_act_orb
     p2 = list_act(pp2)
      !!!! double alpha/beta
      do ii = 1, Ne(2) ! purely closed shell part 
       i = occ(ii,2)
       call give_integrals_3_body(i ,p2,p1,i,h2,h1,integral)
       int_direct = -1.d0 * integral
       call give_integrals_3_body(p1,p2, i,i,h2,h1,integral)
       int_exc_13 = -1.d0 * integral
       call give_integrals_3_body(p2, i,p1,i,h2,h1,integral)
       int_exc_12 = -1.d0 * integral
       normal_two_body_ab(p2,h2,p1,h1) += 2.d0 * int_direct - 1.d0 * ( int_exc_13 + int_exc_12)
      enddo
      do ii = Ne(2) + 1, Ne(1) ! purely open-shell part 
       i = occ(ii,1)
       call give_integrals_3_body(i ,p2,p1,i,h2,h1,integral)
       int_direct = -1.d0 * integral
       call give_integrals_3_body(p1,p2, i,i,h2,h1,integral)
       int_exc_13 = -1.d0 * integral
       call give_integrals_3_body(p2, i,p1,i,h2,h1,integral)
       int_exc_12 = -1.d0 * integral
       normal_two_body_ab(p2,h2,p1,h1) += 1.d0 * int_direct - 0.5d0* ( int_exc_13 + int_exc_12)
      enddo
    enddo
   enddo
  enddo
 enddo

 do hh1 = 1, n_act_orb
  h1 = list_act(hh1) 
  do pp1 = hh1+1, n_act_orb
   p1 = list_act(pp1)
   do hh2 = 1, n_act_orb
    h2 = list_act(hh2) 
    do pp2 = hh2+1, n_act_orb
     p2 = list_act(pp2)
     normal_two_body_ab(h2,p2,h1,p1) = normal_two_body_ab(p2,h2,p1,h1)
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, normal_two_body_aa_bb, (n_act_orb, n_act_orb, n_act_orb, n_act_orb)]
 implicit none
 BEGIN_DOC
! Normal ordered two-body sector of the three-body terms for opposite spin double excitations 
 END_DOC
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,ii,j,h1,p1,h2,p2
 integer :: kkk,k,m,l,n,hh1,hh2,pp1,pp2,kk
 integer                        :: occ(N_int*bit_kind_size,2),Ne(2)
 integer(bit_kind), allocatable :: key_i_core(:,:)
 double precision :: int_direct,int_exc_12,int_exc_13,int_exc_23,hthree
 double precision :: integral,int_exc_l,int_exc_ll
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
 normal_two_body_aa_bb = 0.d0
 do hh1 = 1, n_act_orb
! do hh1 = 1, 1
  h1 = list_act(hh1) 
  do pp1 = 1 , n_act_orb
!  do pp1 = 4,4
   p1 = list_act(pp1)
   do hh2 = 1, n_act_orb
!   do hh2 = 2,2
    h2 = list_act(hh2) 
    do pp2 = 1 , n_act_orb
!    do pp2 = 7,7
     p2 = list_act(pp2)
     call give_aab_contraction(h1,h2,p1,p2,Ne,occ,hthree)
     normal_two_body_aa_bb(p2,h2,p1,h1) = hthree
    enddo
   enddo
  enddo
 enddo

! do hh1 = 1, n_act_orb
!  h1 = list_act(hh1) 
!  do pp1 = hh1+1, n_act_orb
!   p1 = list_act(pp1)
!   do hh2 = 1, n_act_orb
!    h2 = list_act(hh2) 
!    do pp2 = hh2+1, n_act_orb
!     p2 = list_act(pp2)
!     normal_two_body_aa_bb(h2,p2,h1,p1) = normal_two_body_aa_bb(p2,h2,p1,h1)
!    enddo
!   enddo
!  enddo
! enddo

END_PROVIDER 


subroutine update_aa_contraction(h1,h2,p1,p2,Ne,occ,array,n_array)
 implicit none
 integer, intent(in) :: h1,h2,p1,p2
 integer, intent(in) :: Ne(2),occ(N_int*bit_kind_size,2),n_array
 double precision, intent(inout) :: array(n_array,n_array,n_array,n_array)
 integer :: ii,i
 double precision :: int_direct,int_exc_12,int_exc_13,int_exc_23
 double precision :: integral,int_exc_l,int_exc_ll
 do ii = 1, Ne(2) ! purely closed shell part 
  i = occ(ii,2)
  print*,'--'
  print*,i
  call give_integrals_3_body(i ,p2,p1,i,h2,h1,integral)
  int_direct = -1.d0 * integral
  call give_integrals_3_body(p2,p1,i ,i,h2,h1,integral)
  int_exc_l = -1.d0 * integral
  call give_integrals_3_body(p1,i ,p2,i,h2,h1,integral)
  int_exc_ll= -1.d0 * integral
  call give_integrals_3_body(p2,i ,p1,i,h2,h1,integral)
  int_exc_12= -1.d0 * integral
  call give_integrals_3_body(p1,p2, i,i,h2,h1,integral)
  int_exc_13= -1.d0 * integral
  call give_integrals_3_body(i ,p1,p2,i,h2,h1,integral)
  int_exc_23= -1.d0 * integral
  print*,'positive'
  print*,int_direct,int_exc_l,int_exc_ll
  print*,'negative'
  print*,int_exc_12,int_exc_13,int_exc_23
  print*,'--'

  array(p2,h2,p1,h1) +=  & 
  2.d0 * int_direct + int_exc_l + int_exc_ll -( int_exc_12+ int_exc_13+2.d0 * int_exc_23  )
 enddo
 do ii = Ne(2)+1,Ne(1) ! purely open-shell part 
  i = occ(ii,1)
  call give_integrals_3_body(i ,p2,p1,i,h2,h1,integral)
  int_direct = -1.d0 * integral
  call give_integrals_3_body(p2,p1,i ,i,h2,h1,integral)
  int_exc_l = -1.d0 * integral
  call give_integrals_3_body(p1,i ,p2,i,h2,h1,integral)
  int_exc_ll= -1.d0 * integral
  call give_integrals_3_body(p2,i ,p1,i,h2,h1,integral)
  int_exc_12= -1.d0 * integral
  call give_integrals_3_body(p1,p2, i,i,h2,h1,integral)
  int_exc_13= -1.d0 * integral
  call give_integrals_3_body(i ,p1,p2,i,h2,h1,integral)
  int_exc_23= -1.d0 * integral

  array(p2,h2,p1,h1) +=  & 
  1.d0 * int_direct + 0.5d0 * (int_exc_l + int_exc_ll -( int_exc_12+ int_exc_13+2.d0 * int_exc_23  ))
 enddo

end

subroutine give_aab_contraction(h1,h2,p1,p2,Ne,occ,hthree)
 implicit none
 integer, intent(in) :: h1,h2,p1,p2
 integer, intent(in) :: Ne(2),occ(N_int*bit_kind_size,2)
 double precision, intent(inout) :: hthree
 integer :: ii,i
 double precision :: int_direct,int_exc_12,int_exc_13,int_exc_23
 double precision :: integral,int_exc_l,int_exc_ll
 hthree = 0.d0
 do ii = 1, Ne(2) ! purely closed shell part 
  i = occ(ii,2)
  call give_integrals_3_body(p2,p1,i,h2,h1,i,integral)
  int_direct = -1.d0 * integral
  call give_integrals_3_body(p1,p2,i,h2,h1,i,integral)
  int_exc_23= -1.d0 * integral
  hthree  +=  1.d0 * int_direct - int_exc_23
 enddo

end
