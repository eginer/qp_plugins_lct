subroutine provide_all_three_ints
 implicit none
PROVIDE three_body_3_index three_body_3_index_exch_23 three_body_3_index_exch_12 
PROVIDE three_body_3_index_exch_13 three_body_3_index_exch_231 
PROVIDE three_body_4_index three_body_4_index_exch_12_part three_body_4_index_exch_12
PROVIDE three_body_4_index_exch_231 three_body_4_index_exch_312
PROVIDE three_body_4_index_exch_12_part_bis
if(.not.double_normal_ord)then
 PROVIDE three_body_5_index three_body_5_index_132 three_body_5_index_312 
 PROVIDE three_body_5_index_exch_12 three_body_5_index_exch_13 three_body_5_index_exch_32
else
 PROVIDE normal_two_body
endif
end
subroutine diag_htilde_mu_mat_three_body(Nint, key_i, hthree)

  BEGIN_DOC
  !  diagonal element of htilde ONLY FOR THREE-BODY TERMS
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_i(Nint,2)
  double precision, intent(out) :: hthree
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  integer(bit_kind)             :: key_i_core(Nint,2)
  double precision              :: direct_int, exchange_int
  double precision              :: exchange_int_12, exchange_int_13, exchange_int_23
  double precision              :: exchange_int_231,exchange_int_312

  if(core_tc_op)then
   do i = 1, Nint
    key_i_core(i,1) = xor(key_i(i,1),core_bitmask(i,1))
    key_i_core(i,2) = xor(key_i(i,2),core_bitmask(i,2))
   enddo
   call bitstring_to_list_ab(key_i_core,occ,Ne,Nint)
  else
   call bitstring_to_list_ab(key_i,occ,Ne,Nint)
  endif
  hthree = 0.d0

  if(Ne(1)+Ne(2).ge.3)then
!!  ! alpha/alpha/beta three-body
   do i = 1, Ne(1)
    ii = occ(i,1) 
    do j = i+1, Ne(1)
     jj = occ(j,1) 
     do k = 1, Ne(2)
      kk = occ(k,2) 
      direct_int = three_body_3_index(kk,jj,ii)
      exchange_int_23 = three_body_3_index_exch_23(kk,jj,ii)
      hthree += direct_int - exchange_int_23
     enddo
    enddo
   enddo
  
   ! beta/beta/alpha three-body
   do i = 1, Ne(2)
    ii = occ(i,2) 
    do j = i+1, Ne(2)
     jj = occ(j,2) 
     do k = 1, Ne(1)
      kk = occ(k,1) 
      direct_int = three_body_3_index(kk,jj,ii)
      exchange_int_23 = three_body_3_index_exch_23(kk,jj,ii)
      hthree += direct_int - exchange_int_23
     enddo
    enddo
   enddo

   ! alpha/alpha/alpha three-body
   do i = 1, Ne(1)
    ii = occ(i,1) ! 1
    do j = i+1, Ne(1)
     jj = occ(j,1) ! 2 
     do k = j+1, Ne(1)
      kk = occ(k,1) ! 3 
      !               direct   :       3  2  1  
      direct_int = three_body_3_index(kk,jj,ii)
      exchange_int_231 = three_body_3_index_exch_231(kk,jj,ii)
      exchange_int_312 = three_body_3_index_exch_231(kk,jj,ii)
      exchange_int_12  = three_body_3_index_exch_12(kk,jj,ii)
      exchange_int_13  = three_body_3_index_exch_13(kk,jj,ii)
      exchange_int_23  = three_body_3_index_exch_23(kk,jj,ii)
      hthree +=  direct_int + exchange_int_231 + exchange_int_312 & 
               - exchange_int_12 - exchange_int_13 - exchange_int_23  
     enddo
    enddo
   enddo

   ! beta/beta/beta three-body
   do i = 1, Ne(2)
    ii = occ(i,2) ! 1
    do j = i+1, Ne(2)
     jj = occ(j,2) ! 2
     do k = j+1, Ne(2)
      kk = occ(k,2) ! 3
      !               direct   :       3  2  1  
      direct_int = three_body_3_index(kk,jj,ii)
      exchange_int_231 = three_body_3_index_exch_231(kk,jj,ii)
      exchange_int_312 = three_body_3_index_exch_231(kk,jj,ii)
      exchange_int_12 = three_body_3_index_exch_12(kk,jj,ii)
      exchange_int_13 = three_body_3_index_exch_13(kk,jj,ii)
      exchange_int_23 = three_body_3_index_exch_23(kk,jj,ii)
      hthree +=  direct_int + exchange_int_231 + exchange_int_312 & 
               - exchange_int_12 - exchange_int_13 - exchange_int_23  
     enddo
    enddo
   enddo
  endif

end



subroutine single_htilde_mu_mat_three_body(Nint, key_j, key_i, hthree)

  BEGIN_DOC
  ! <key_j | H_tilde | key_i> for single excitation ONLY FOR THREE-BODY TERMS
  !!
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2),key_i(Nint,2)
  double precision, intent(out) :: hthree
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  double precision              :: direct_int,exchange_int_12,exchange_int_23,exchange_int_13,phase
  double precision              :: exchange_int_231,exchange_int_312
  integer                       :: other_spin(2)
  integer(bit_kind)             :: key_j_core(Nint,2),key_i_core(Nint,2)

  other_spin(1) = 2
  other_spin(2) = 1


  hthree = 0.d0
  call get_excitation_degree(key_i,key_j,degree,Nint)
  if(degree.ne.1)then
   return
  endif
  if(core_tc_op)then
   do i = 1, Nint
    key_i_core(i,1) = xor(key_i(i,1),core_bitmask(i,1))
    key_i_core(i,2) = xor(key_i(i,2),core_bitmask(i,2))
    key_j_core(i,1) = xor(key_j(i,1),core_bitmask(i,1))
    key_j_core(i,2) = xor(key_j(i,2),core_bitmask(i,2))
   enddo
   call bitstring_to_list_ab(key_i_core, occ, Ne, Nint)
  else
   call bitstring_to_list_ab(key_i, occ, Ne, Nint)
  endif

  call get_single_excitation(key_i, key_j, exc, phase, Nint)
  call decode_exc(exc, 1, h1, p1, h2, p2, s1, s2)

   ! alpha/alpha/beta three-body
   if(Ne(1)+Ne(2).ge.3)then
    if(s1 == 2)then ! single beta 
     ! alpha-alpha + hole/particle beta 
     do i = 1, Ne(1)
      ii = occ(i,1) 
      do j = i+1, Ne(1)
       jj = occ(j,1) 
       !                         b  a a     b a a       b  a a   b a a
       !                       < h1 j  i | p1 j i > - < h1 j i | p1 i j >
       direct_int  = three_body_4_index(jj,ii,h1,p1)
       exchange_int_12 = three_body_4_index_exch_12(jj,ii,h1,p1)
       hthree += direct_int - exchange_int_12
      enddo
     enddo
  
     ! alpha-beta + hole/particle beta
     do i = 1, Ne(1)
      ii = occ(i,1) 
      do j = 1, Ne(2)
       jj = occ(j,2) 
       direct_int  = three_body_4_index(jj,ii,h1,p1)
       exchange_int_12 = three_body_4_index_exch_12_part(jj,ii,h1,p1)
       !                         b  b a   b b a         b  b a   b b a
       !                       < h1 j  i | p1 j i > - < h1 j i | j p1 i >
       hthree += direct_int - exchange_int_12
      enddo
     enddo

     ! beta-beta-beta
     do i = 1, Ne(2)
      ii = occ(i,2)
      do j = i+1, Ne(2)
       jj = occ(j,2)
       direct_int = three_body_4_index(jj,ii,h1,p1)                        ! < h1 jj ii | p1 jj ii >
       exchange_int_231 = three_body_4_index_exch_231(jj,ii,h1,p1)
       exchange_int_312 = three_body_4_index_exch_312(jj,ii,h1,p1)
       exchange_int_23 = three_body_4_index_exch_12(jj,ii,h1,p1)           ! < h1 jj ii | p1 ii jj >
       exchange_int_12 = three_body_4_index_exch_12_part(ii,jj,h1,p1)      ! < h1 jj ii | ii p1 jj >
       exchange_int_13 = three_body_4_index_exch_12_part_bis(ii,jj,h1,p1)  ! < h1 jj ii | ii p1 jj >
       hthree += direct_int + exchange_int_231 + exchange_int_312 & 
              -  exchange_int_23 & ! ii <-> jj
              -  exchange_int_12 & ! p1 <-> jj
              -  exchange_int_13   ! p1 <-> ii
      enddo
     enddo
  
    else ! single alpha 
     ! beta-beta + hole/particle alpha 
     do i = 1, Ne(2)
      ii = occ(i,2) 
      do j = i+1, Ne(2)
       jj = occ(j,2)
       direct_int  = three_body_4_index(jj,ii,h1,p1)
       exchange_int_12 = three_body_4_index_exch_12(jj,ii,h1,p1)
       !                         a  b b   a  b b       a  b b   a  b b
       !                       < h1 j i | p1 j i > - < h1 j i | p1 i j >
       hthree += direct_int - exchange_int_12
      enddo
     enddo
     ! alpha-beta + hole/particle alpha 
     do i = 1, Ne(2)
      ii = occ(i,2) 
      do j = 1, Ne(1)
       jj = occ(j,1)
       direct_int  = three_body_4_index(jj,ii,h1,p1)
       exchange_int_12 = three_body_4_index_exch_12_part(jj,ii,h1,p1)
       !                         a  a b   a  a b                       a  a b   a  a b   
       !                       < h1 j i | p1 j i >  -                < h1 j i | j p1 i >  
       hthree += direct_int - exchange_int_12
      enddo
     enddo

     ! alpha-alpha-alpha
     do i = 1, Ne(1)
      ii = occ(i,1)
      do j = i+1, Ne(1)
       jj = occ(j,1)
       direct_int = three_body_4_index(jj,ii,h1,p1)                       ! < h1 jj ii | p1 jj ii >
       exchange_int_231 = three_body_4_index_exch_231(jj,ii,h1,p1)        ! < h1 jj ii | ii p1 jj >
       exchange_int_312 = three_body_4_index_exch_312(jj,ii,h1,p1)        ! < h1 jj ii | jj ii p1 >
       exchange_int_23 = three_body_4_index_exch_12(jj,ii,h1,p1)          ! < h1 jj ii | p1 ii jj >
       exchange_int_12 = three_body_4_index_exch_12_part(ii,jj,h1,p1)     ! < h1 jj ii | ii jj p1 >
       exchange_int_13 = three_body_4_index_exch_12_part_bis(ii,jj,h1,p1) ! < h1 jj ii | jj p1 ii >
       hthree += direct_int + exchange_int_231 + exchange_int_312 & 
              -  exchange_int_23 & ! ii <-> jj
              -  exchange_int_12 & ! p1 <-> jj
              -  exchange_int_13   ! p1 <-> ii
      enddo
     enddo
  
    endif
   endif

  hthree  *= phase

end



subroutine double_htilde_mu_mat_three_body(Nint, key_j, key_i, hthree)

  BEGIN_DOC
  ! <key_j | H_tilde | key_i> for double excitation ONLY FOR THREE-BODY TERMS 
  !!
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2),key_i(Nint,2)
  double precision, intent(out) :: hthree
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  double precision              :: phase
  integer                       :: other_spin(2)
  integer(bit_kind)             :: key_i_core(Nint,2)
  double precision              :: integral,integral_exch

  other_spin(1) = 2
  other_spin(2) = 1

  call get_excitation_degree(key_i, key_j, degree, Nint)

  hthree = 0.d0

  if(degree.ne.2)then
   return
  endif

  if(core_tc_op)then
   do i = 1, Nint
    key_i_core(i,1) = xor(key_i(i,1),core_bitmask(i,1))
    key_i_core(i,2) = xor(key_i(i,2),core_bitmask(i,2))
   enddo
   call bitstring_to_list_ab(key_i_core, occ, Ne, Nint)
  else
   call bitstring_to_list_ab(key_i, occ, Ne, Nint)
  endif
  call get_double_excitation(key_i, key_j, exc, phase, Nint)
  call decode_exc(exc, 2, h1, p1, h2, p2, s1, s2)

    ! alpha/alpha/beta threee-body 
    if(Ne(1)+Ne(2).ge.3)then
     if(s1.eq.s2.and.s2.eq.1)then ! double alpha 
      do k = 1, Ne(2) ! beta - alpha/alpha
       kk = occ(k,2)
       hthree += three_body_5_index(kk,h1,h2,p1,p2) - three_body_5_index_exch_12(kk,h1,h2,p1,p2)
      enddo
      if(Ne(1).ge.3)then
       do k = 1, Ne(1) ! alpha/alpha/alpha
        kk = occ(k,1)
        hthree +=  three_body_5_index(kk,h1,h2,p1,p2)
        hthree +=  three_body_5_index_132(kk,h1,h2,p1,p2)
        hthree +=  three_body_5_index_312(kk,h1,h2,p1,p2)
        hthree -=  three_body_5_index_exch_12(kk,h1,h2,p1,p2)
        hthree -=  three_body_5_index_exch_13(kk,h1,h2,p1,p2)
        hthree -=  three_body_5_index_exch_32(kk,h1,h2,p1,p2)
       enddo 
      endif
     else if(s1.eq.s2.and.s2.eq.2)then ! double beta 
      do k = 1, Ne(1) ! alpha- beta/beta
       kk = occ(k,1)
       hthree += three_body_5_index(kk,h1,h2,p1,p2) - three_body_5_index_exch_12(kk,h1,h2,p1,p2)
      enddo
      if(Ne(2).ge.3)then
       do k = 1, Ne(2) ! beta/beta/beta
        kk = occ(k,2)
        hthree +=  three_body_5_index(kk,h1,h2,p1,p2)
        hthree +=  three_body_5_index_132(kk,h1,h2,p1,p2)
        hthree +=  three_body_5_index_312(kk,h1,h2,p1,p2)
        hthree -=  three_body_5_index_exch_12(kk,h1,h2,p1,p2)
        hthree -=  three_body_5_index_exch_13(kk,h1,h2,p1,p2)
        hthree -=  three_body_5_index_exch_32(kk,h1,h2,p1,p2)
       enddo 
      endif
     else ! double alpha/beta 
      if(s1.eq.1.and.s2.eq.2)then ! s1 == alpha , s2 == beta 
       do k = 1, Ne(1)
        kk = occ(k,1) ! direct - exchange in alpha 
        hthree += three_body_5_index(kk,h1,h2,p1,p2) - three_body_5_index_exch_13(kk,h1,h2,p1,p2)
       enddo
       do k = 1, Ne(2)
        kk = occ(k,2)! direct - exchange in beta 
        hthree +=  three_body_5_index(kk,h1,h2,p1,p2) - three_body_5_index_exch_32(kk,h1,h2,p1,p2)
       enddo
      else if(s1.eq.2.and.s2.eq.1)then  ! s1 == beta, s2 == alpha 
       do k = 1, Ne(2)
        kk = occ(k,2) ! direct - exchange in beta 
        hthree +=  three_body_5_index(kk,h1,h2,p1,p2) - three_body_5_index_exch_13(kk,h1,h2,p2,p1)
       enddo
       do k = 1, Ne(1)
        kk = occ(k,1)! direct - exchange in alpha 
        hthree +=  three_body_5_index(kk,h1,h2,p1,p2) - three_body_5_index_exch_13(kk,h1,h2,p2,p1)
       enddo
      endif 
     endif
    endif
  hthree  *= phase
 end


