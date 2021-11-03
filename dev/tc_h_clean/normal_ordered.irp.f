BEGIN_PROVIDER [ double precision, normal_two_body, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC 
! Normal ordering of the three body interaction on the HF density
 END_DOC 
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,h1,p1,h2,p2
 integer :: hh1,hh2,pp1,pp2
 integer                        :: occ(N_int*bit_kind_size,2),Ne(2)
 integer(bit_kind), allocatable :: key_i_core(:,:)
 double precision :: hthree_aba,hthree_aaa,hthree_aab
 double precision :: wall0,wall1
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
 print*,'Providing normal_two_body ...'
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (hh1,h1,hh2,h2,pp1,p1,pp2,p2,hthree_aba,hthree_aab,hthree_aaa) & 
 !$OMP SHARED (n_act_orb,list_act,Ne,occ,normal_two_body)
 !$OMP DO SCHEDULE (static) 
 do hh1 = 1, n_act_orb
  h1 = list_act(hh1) 
  do pp1 = hh1, n_act_orb
   p1 = list_act(pp1)
   do hh2 = 1, n_act_orb
    h2 = list_act(hh2) 
    do pp2 = hh2, n_act_orb
     p2 = list_act(pp2)
     ! opposite spin double excitations 
     call give_aba_contraction(h1,h2,p1,p2,Ne,occ,hthree_aba)
     ! same spin double excitations with opposite spin contributions 
     call give_aab_contraction(h1,h2,p1,p2,Ne,occ,hthree_aab)
     ! same spin double excitations with same spin contributions 
     if(Ne(2).ge.3)then
      call give_aaa_contraction(h1,h2,p1,p2,Ne,occ,hthree_aaa)
     else
      hthree_aaa = 0.d0
     endif
     normal_two_body(p2,h2,p1,h1) = 0.5d0*(hthree_aba + hthree_aab + hthree_aaa)
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 do hh1 = 1, n_act_orb
  h1 = list_act(hh1) 
  do pp1 = hh1+1, n_act_orb
   p1 = list_act(pp1)
   do hh2 = 1, n_act_orb
    h2 = list_act(hh2) 
    do pp2 = hh2+1, n_act_orb
     p2 = list_act(pp2)
     normal_two_body(h2,p2,h1,p1) = normal_two_body(p2,h2,p1,h1)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall1)
 print*,'Wall time for normal_two_body ',wall1-wall0

END_PROVIDER 


subroutine give_aba_contraction(h1,h2,p1,p2,Ne,occ,hthree)
 implicit none
 integer, intent(in) :: h1,h2,p1,p2
 integer, intent(in) :: Ne(2),occ(N_int*bit_kind_size,2)
 double precision, intent(out) :: hthree
 double precision :: int_direct,int_exc_12,int_exc_13,integral
 integer :: ii,i
 !!!! double alpha/beta
 hthree = 0.d0
 do ii = 1, Ne(2) ! purely closed shell part 
  i = occ(ii,2)
  call give_integrals_3_body(i ,p2,p1,i,h2,h1,integral)
  int_direct = -1.d0 * integral
  call give_integrals_3_body(p1,p2, i,i,h2,h1,integral)
  int_exc_13 = -1.d0 * integral
  call give_integrals_3_body(p2, i,p1,i,h2,h1,integral)
  int_exc_12 = -1.d0 * integral
  hthree += 2.d0 * int_direct - 1.d0 * ( int_exc_13 + int_exc_12)
 enddo
 do ii = Ne(2) + 1, Ne(1) ! purely open-shell part 
  i = occ(ii,1)
  call give_integrals_3_body(i ,p2,p1,i,h2,h1,integral)
  int_direct = -1.d0 * integral
  call give_integrals_3_body(p1,p2, i,i,h2,h1,integral)
  int_exc_13 = -1.d0 * integral
  call give_integrals_3_body(p2, i,p1,i,h2,h1,integral)
  int_exc_12 = -1.d0 * integral
  hthree += 1.d0 * int_direct - 0.5d0* ( int_exc_13 + int_exc_12)
 enddo

end

BEGIN_PROVIDER [ double precision, normal_two_body_ab, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! Normal ordered two-body sector of the three-body terms for opposite spin double excitations 
 END_DOC
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: h1,p1,h2,p2,i
 integer :: hh1,hh2,pp1,pp2
 integer                        :: occ(N_int*bit_kind_size,2),Ne(2)
 integer(bit_kind), allocatable :: key_i_core(:,:)
 double precision :: hthree
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
     call give_aba_contraction(h1,h2,p1,p2,Ne,occ,hthree)
     normal_two_body_ab(p2,h2,p1,h1) = hthree    
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
 integer :: hh1,hh2,pp1,pp2
 integer                        :: occ(N_int*bit_kind_size,2),Ne(2)
 integer(bit_kind), allocatable :: key_i_core(:,:)
 double precision :: hthree_aab,hthree_aaa
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
  h1 = list_act(hh1) 
  do pp1 = hh1 , n_act_orb
   p1 = list_act(pp1)
   do hh2 = 1, n_act_orb
    h2 = list_act(hh2) 
    do pp2 = hh2 , n_act_orb
     p2 = list_act(pp2)
     call give_aab_contraction(h1,h2,p1,p2,Ne,occ,hthree_aab)
     if(Ne(2).ge.3)then
      call give_aaa_contraction(h1,h2,p1,p2,Ne,occ,hthree_aaa)
     else
      hthree_aaa = 0.d0
     endif
     normal_two_body_aa_bb(p2,h2,p1,h1) = hthree_aab + hthree_aaa
    enddo
   enddo
  enddo
 enddo

 do hh1 = 1, n_act_orb
  h1 = list_act(hh1) 
  do pp1 = hh1, n_act_orb
   p1 = list_act(pp1)
   do hh2 = 1, n_act_orb
    h2 = list_act(hh2) 
    do pp2 = hh2, n_act_orb
     p2 = list_act(pp2)
     normal_two_body_aa_bb(h2,p2,h1,p1) = normal_two_body_aa_bb(p2,h2,p1,h1)
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 


subroutine give_aaa_contraction(h1,h2,p1,p2,Ne,occ,hthree)
 implicit none
 integer, intent(in) :: h1,h2,p1,p2
 integer, intent(in) :: Ne(2),occ(N_int*bit_kind_size,2)
 double precision, intent(out) :: hthree
 integer :: ii,i
 double precision :: int_direct,int_exc_12,int_exc_13,int_exc_23
 double precision :: integral,int_exc_l,int_exc_ll
 hthree = 0.d0
 do ii = 1, Ne(2) ! purely closed shell part 
  i = occ(ii,2)
!  print*,'--'
!  print*,i
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
!  print*,'positive'
!  print*,int_direct,int_exc_l,int_exc_ll
!  print*,'negative'
!  print*,int_exc_12,int_exc_13,int_exc_23
!  print*,'--'

  hthree +=  & 
  1.d0 * int_direct + int_exc_l + int_exc_ll -( int_exc_12+ int_exc_13+ int_exc_23  )
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

  hthree +=  & 
  1.d0 * int_direct + 0.5d0 * (int_exc_l + int_exc_ll -( int_exc_12+ int_exc_13+ int_exc_23  ))
 enddo

end

subroutine give_aab_contraction(h1,h2,p1,p2,Ne,occ,hthree)
 implicit none
 integer, intent(in) :: h1,h2,p1,p2
 integer, intent(in) :: Ne(2),occ(N_int*bit_kind_size,2)
 double precision, intent(out) :: hthree
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
