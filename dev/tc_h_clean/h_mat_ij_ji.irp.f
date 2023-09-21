
subroutine double_hij_ji_mu_mat_scal_map(Nint, key_j, key_i, hmono, heff, hderiv_ji, hderiv_ij, htot_ij, htot_ji)

  BEGIN_DOC
  ! <key_j | H_tilde | key_i> for double excitation  ONLY FOR ONE- AND TWO-BODY TERMS 
  !!
  ! returns both htot_ji =  <key_j | H_tilde | key_i>
  !
  ! AND          htot_ij =  <key_i | H_tilde | key_j>
  !
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint 
  integer(bit_kind), intent(in) :: key_j(Nint,2), key_i(Nint,2)
  double precision, intent(out) :: hmono, heff, hderiv_ij, hderiv_ji, htot_ij, htot_ji
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2), i, j, ii, jj, ispin, jspin, k, kk
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  integer                       :: other_spin(2)
  integer(bit_kind)             :: key_i_core(Nint,2)
  double precision              :: get_mo_two_e_integral_tc_int,phase

  PROVIDE mo_two_e_integrals_tc_int_in_map mo_non_hermit_term

  other_spin(1) = 2
  other_spin(2) = 1

  call get_excitation_degree(key_i, key_j, degree, Nint)

  hmono = 0.d0
  heff  = 0.d0
  hderiv_ji= 0.d0
  hderiv_ij= 0.d0
  htot_ji = 0.d0
  htot_ij = 0.d0

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

  if(s1.ne.s2)then
   ! opposite spin two-body 
   if(s1==1)then
    heff    = get_mo_two_e_integral_tc_int(p2,p1,h2,h1,mo_integrals_tc_int_map) 
    if(double_normal_ord.and.+Ne(1).gt.2)then
     heff += normal_two_body(h1,p1,h2,p2)
    endif
    hderiv_ji  = mo_non_hermit_term(p2,p1,h2,h1) 
    hderiv_ij  = mo_non_hermit_term(h2,h1,p2,p1) 
   else
    heff    = get_mo_two_e_integral_tc_int(p1,p2,h1,h2,mo_integrals_tc_int_map) 
    if(double_normal_ord.and.Ne(2)+Ne(1).gt.2)then
     heff += normal_two_body(h1,p1,h2,p2)
    endif
    hderiv_ji  = mo_non_hermit_term(p1,p2,h1,h2) 
    hderiv_ij  = mo_non_hermit_term(h1,h2,p1,p2) 
   endif
  else
   ! same spin two-body 
   ! direct terms 
   heff    = get_mo_two_e_integral_tc_int(p2,p1,h2,h1,mo_integrals_tc_int_map) 
   if(double_normal_ord.and.+Ne(1).gt.2)then
    heff += normal_two_body(h1,p1,h2,p2)
   endif
   hderiv_ji  = mo_non_hermit_term(p2,p1,h2,h1)  
   hderiv_ij  = mo_non_hermit_term(h2,h1,p2,p1)  
   ! exchange terms 
   heff   -= get_mo_two_e_integral_tc_int(p1,p2,h2,h1,mo_integrals_tc_int_map) 
   if(double_normal_ord.and.+Ne(1).gt.2)then
    heff -= normal_two_body(h2,p1,h1,p2)
   endif
   hderiv_ji -= mo_non_hermit_term(p1,p2,h2,h1) 
   hderiv_ij -= mo_non_hermit_term(h1,h2,p2,p1) 
  endif
  heff   *= phase
  hderiv_ji *= phase
  hderiv_ij *= phase
  htot_ji = heff + hderiv_ji 
  htot_ij = heff + hderiv_ij 

end


subroutine single_hij_ji_mu_mat_scal_map(Nint, key_j, key_i, hmono_ji,hmono_ij, heff, hderiv_ji, hderiv_ij, htot_ij, htot_ji)

  BEGIN_DOC
  ! <key_j | H_tilde | key_i> for single excitation ONLY FOR ONE- AND TWO-BODY TERMS 
  !!
  ! returns both htot_ji =  <key_j | H_tilde | key_i>
  !
  ! AND          htot_ij =  <key_i | H_tilde | key_j>
  !
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2), key_i(Nint,2)
  double precision, intent(out) :: hmono_ij, hmono_ji,heff, hderiv_ij, hderiv_ji, htot_ij, htot_ji
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2), i, j, ii, jj, ispin, jspin, k, kk
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  double precision              :: get_mo_two_e_integral_tc_int, phase
  double precision              :: direct_int, exchange_int_12, exchange_int_23, exchange_int_13
  integer                       :: other_spin(2)
  integer(bit_kind)             :: key_j_core(Nint,2), key_i_core(Nint,2)

  PROVIDE mo_two_e_integrals_tc_int_in_map mo_non_hermit_term

  PROVIDE core_bitmask core_fock_operator mo_integrals_erf_map

  PROVIDE j1b_gauss

  other_spin(1) = 2
  other_spin(2) = 1

  hmono_ji = 0.d0
  hmono_ij = 0.d0
  heff  = 0.d0
  hderiv_ji= 0.d0
  hderiv_ij= 0.d0
  htot_ji = 0.d0
  call get_excitation_degree(key_i, key_j, degree, Nint)
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
  call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)

  hmono_ji = mo_one_e_integrals(h1,p1) * phase
  hmono_ij = mo_one_e_integrals(p1,h1) * phase

  if( j1b_gauss .eq. 1 ) then
    hmono_ji += ( mo_j1b_gauss_hermI  (h1,p1) &
             +     mo_j1b_gauss_hermII(h1,p1) &
             +    mo_j1b_gauss_nonherm(h1,p1) ) * phase
    hmono_ij += ( mo_j1b_gauss_hermI  (p1,h1) &
             +     mo_j1b_gauss_hermII(p1,h1) &
             +    mo_j1b_gauss_nonherm(p1,h1) ) * phase
  endif

  if(core_tc_op)then
   hmono_ji += phase * core_fock_operator(h1,p1)
   hmono_ij += phase * core_fock_operator(p1,h1)
  endif
  
   ! alpha/beta two-body 
   ispin = other_spin(s1)
   if(s1==1)then
    ! single alpha 
    do i = 1, Ne(ispin) ! electron 2 
     ii = occ(i,ispin) 
     heff   += get_mo_two_e_integral_tc_int(ii,p1,ii,h1,mo_integrals_tc_int_map) 
     hderiv_ji += mo_non_hermit_term(ii,p1,ii,h1) 
     hderiv_ij += mo_non_hermit_term(ii,h1,ii,p1) 
    enddo
   else
    ! single beta 
    do i = 1, Ne(ispin) ! electron 1 
     ii = occ(i,ispin) 
     heff   += get_mo_two_e_integral_tc_int(p1,ii,h1,ii,mo_integrals_tc_int_map) 
     hderiv_ji += mo_non_hermit_term(p1,ii,h1,ii) 
     hderiv_ij += mo_non_hermit_term(h1,ii,p1,ii) 
    enddo
   endif
!   ! same spin two-body 
   do i = 1, Ne(s1)
    ii = occ(i,s1) 
    ! (h1p1|ii ii) - (h1 ii| p1 ii)
    heff   += get_mo_two_e_integral_tc_int(ii,p1,ii,h1,mo_integrals_tc_int_map) & 
             -get_mo_two_e_integral_tc_int(p1,ii,ii,h1,mo_integrals_tc_int_map)
    hderiv_ji += mo_non_hermit_term(ii,p1,ii,h1) - mo_non_hermit_term(p1,ii,ii,h1) 
    hderiv_ij += mo_non_hermit_term(ii,h1,ii,p1) - mo_non_hermit_term(h1,ii,ii,p1) 
   enddo
   
  heff    *= phase
  hderiv_ji  *= phase
  hderiv_ij  *= phase
  htot_ji = hmono_ji + heff + hderiv_ji 
  htot_ij = hmono_ij + heff + hderiv_ij 

end



