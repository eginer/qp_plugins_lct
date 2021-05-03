
subroutine diag_htilde_mu_mat_3_index(key_i,hmono,herf,heff,hderiv,hthree,htot)
  use bitmasks
  BEGIN_DOC
!  diagonal element of htilde 
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv,htot,hthree
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  double precision :: get_mo_two_e_integral_erf,mo_two_e_integral_eff_pot
  double precision :: get_two_e_integral
  PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals 
  PROVIDE mo_two_e_integrals_eff_pot_in_map mo_two_e_integrals_erf_in_map
  integer(bit_kind) :: key_i_core(N_int,2)
  if(core_tc_op)then
   do i = 1, N_int
    key_i_core(i,1) = xor(key_i(i,1),core_bitmask(i,1))
    key_i_core(i,2) = xor(key_i(i,2),core_bitmask(i,2))
   enddo
   call bitstring_to_list_ab(key_i_core,occ,Ne,N_int)
   hmono = core_energy - nuclear_repulsion
  else
   call bitstring_to_list_ab(key_i,occ,Ne,N_int)
   hmono = 0.d0
  endif
  herf  = 0.d0
  heff  = 0.d0
  hderiv= 0.d0
  hthree = 0.d0
  htot = 0.d0

  do ispin = 1, 2 
   do i = 1, Ne(ispin) ! 
    ii = occ(i,ispin) 
    hmono += mo_one_e_integrals(ii,ii)
    if(core_tc_op)then
     hmono += core_fock_operator(ii,ii) ! add the usual Coulomb - Exchange from the core 
    endif
   enddo
  enddo


  if(.not.adjoint_tc_h)then ! Usual transcorrelated Hamiltonian 
   ! alpha/beta two-body
   ispin = 1
   jspin = 2 
   do i = 1, Ne(ispin) 
    ii = occ(i,ispin) 
    do j = 1, Ne(jspin)
     jj = occ(j,jspin) 
     herf += get_mo_two_e_integral_erf(ii,jj,ii,jj,mo_integrals_erf_map)
     heff += mo_two_e_integral_eff_pot(ii,jj,ii,jj) 
     hderiv += mo_two_e_eff_dr12_pot_array_physicist(ii,jj,ii,jj) 
    enddo
   enddo
 
   ! alpha/alpha two-body
   do i = 1, Ne(ispin)
    ii = occ(i,ispin) 
    do j = i+1, Ne(ispin)
     jj = occ(j,ispin) 
     herf += get_mo_two_e_integral_erf(ii,jj,ii,jj,mo_integrals_erf_map) - get_mo_two_e_integral_erf(ii,jj,jj,ii,mo_integrals_erf_map)
     heff += mo_two_e_integral_eff_pot(ii,jj,ii,jj) - mo_two_e_integral_eff_pot(ii,jj,jj,ii)
     hderiv += mo_two_e_eff_dr12_pot_array_physicist(ii,jj,ii,jj) & 
     - 0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(ii,jj,jj,ii) - 0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(jj,ii,ii,jj) 
    enddo
   enddo
 
   ! beta/beta two-body
   do i = 1, Ne(jspin)
    ii = occ(i,jspin) 
    do j = i+1, Ne(jspin)
     jj = occ(j,jspin) 
     herf += get_mo_two_e_integral_erf(ii,jj,ii,jj,mo_integrals_erf_map) - get_mo_two_e_integral_erf(ii,jj,jj,ii,mo_integrals_erf_map)
     heff += mo_two_e_integral_eff_pot(ii,jj,ii,jj) - mo_two_e_integral_eff_pot(ii,jj,jj,ii)
     hderiv += mo_two_e_eff_dr12_pot_array_physicist(ii,jj,ii,jj) & 
     - 0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(ii,jj,jj,ii) & 
     - 0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(jj,ii,ii,jj) 
    enddo
   enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!! ADJOINT OF THE TC HAMILTONIAN !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!! ADJOINT OF THE TC HAMILTONIAN !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!! ADJOINT OF THE TC HAMILTONIAN !!!!!!!!!!!!!!!!!!!!!!!
  else 
   ! alpha/beta two-body
   ispin = 1
   jspin = 2 
   do i = 1, Ne(ispin) 
    ii = occ(i,ispin) 
    do j = 1, Ne(jspin)
     jj = occ(j,jspin) 
     herf += 2.d0 * get_two_e_integral(ii,jj,ii,jj,mo_integrals_map) ! 2 / r12
     herf -= get_mo_two_e_integral_erf(ii,jj,ii,jj,mo_integrals_erf_map) ! - erf(mu r12)/r12
     heff += mo_two_e_integral_eff_pot(ii,jj,ii,jj) ! thanks to "adjoint_tc_h" keyword, proper eff two e pot
     hderiv -= mo_two_e_eff_dr12_pot_array_physicist(ii,jj,ii,jj)  !MARK THE MINUS SIGN HERE 
    enddo
   enddo
 
   ! alpha/alpha two-body
   do i = 1, Ne(ispin)
    ii = occ(i,ispin) 
    do j = i+1, Ne(ispin)
     jj = occ(j,ispin) 
     herf += 2.d0 * (get_two_e_integral(ii,jj,ii,jj,mo_integrals_map) - get_two_e_integral(ii,jj,jj,ii,mo_integrals_map))
     herf -= get_mo_two_e_integral_erf(ii,jj,ii,jj,mo_integrals_erf_map) - get_mo_two_e_integral_erf(ii,jj,jj,ii,mo_integrals_erf_map)
     heff += mo_two_e_integral_eff_pot(ii,jj,ii,jj) - mo_two_e_integral_eff_pot(ii,jj,jj,ii)
     hderiv -= mo_two_e_eff_dr12_pot_array_physicist(ii,jj,ii,jj) & 
     - 0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(ii,jj,jj,ii) - 0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(jj,ii,ii,jj) 
    enddo
   enddo
 
   ! beta/beta two-body
   do i = 1, Ne(jspin)
    ii = occ(i,jspin) 
    do j = i+1, Ne(jspin)
     jj = occ(j,jspin) 
     herf -= get_mo_two_e_integral_erf(ii,jj,ii,jj,mo_integrals_erf_map) - get_mo_two_e_integral_erf(ii,jj,jj,ii,mo_integrals_erf_map)
     herf += 2.d0 * (get_two_e_integral(ii,jj,ii,jj,mo_integrals_map) - get_two_e_integral(ii,jj,jj,ii,mo_integrals_map))
     heff += mo_two_e_integral_eff_pot(ii,jj,ii,jj) - mo_two_e_integral_eff_pot(ii,jj,jj,ii)
     hderiv -= mo_two_e_eff_dr12_pot_array_physicist(ii,jj,ii,jj) & 
     - 0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(ii,jj,jj,ii) - 0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(jj,ii,ii,jj) 
    enddo
   enddo
  endif

  if(three_body_h_tc)then
   if(Ne(1)+Ne(2).ge.3)then
    double precision :: direct_int, exchange_int
    double precision :: exchange_int_12, exchange_int_13, exchange_int_23
!!!  ! alpha/alpha/beta three-body
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
       exchange_int_12 = three_body_3_index_exch_12(kk,jj,ii)
       exchange_int_13 = three_body_3_index_exch_13(kk,jj,ii)
       exchange_int_23 = three_body_3_index_exch_23(kk,jj,ii)
       hthree +=  direct_int - exchange_int_12 - exchange_int_13 - exchange_int_23 
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
       exchange_int_12 = three_body_3_index_exch_12(kk,jj,ii)
       exchange_int_13 = three_body_3_index_exch_13(kk,jj,ii)
       exchange_int_23 = three_body_3_index_exch_23(kk,jj,ii)
       hthree +=  direct_int - exchange_int_12 - exchange_int_13 - exchange_int_23 
      enddo
     enddo
    enddo


   endif
  endif
  htot = hmono + herf + heff + hderiv + hthree

end

subroutine double_htilde_mu_mat_5_index(key_j,key_i,hmono,herf,heff,hderiv,hthree,htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> for double excitation  
!!
!! WARNING !!
! 
! Non hermitian !!
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv,hthree,htot
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  integer                        :: degree,exc(0:2,2,2)
  integer                        :: h1, p1, h2, p2, s1, s2
  double precision :: get_mo_two_e_integral_erf,mo_two_e_integral_eff_pot,phase
  double precision :: get_two_e_integral
  integer :: other_spin(2)
  PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals 
  PROVIDE mo_two_e_integrals_eff_pot_in_map mo_two_e_integrals_erf_in_map
  other_spin(1) = 2
  other_spin(2) = 1

  integer(bit_kind) :: key_i_core(N_int,2)
  call get_excitation_degree(key_i,key_j,degree,N_int)

  hmono = 0.d0
  herf  = 0.d0
  heff  = 0.d0
  hderiv= 0.d0
  hthree = 0.d0
  htot = 0.d0

  if(degree.ne.2)then
   return
  endif

  if(core_tc_op)then
   do i = 1, N_int
    key_i_core(i,1) = xor(key_i(i,1),core_bitmask(i,1))
    key_i_core(i,2) = xor(key_i(i,2),core_bitmask(i,2))
   enddo
   call bitstring_to_list_ab(key_i_core,occ,Ne,N_int)
  else
   call bitstring_to_list_ab(key_i,occ,Ne,N_int)
  endif
  call get_double_excitation(key_i,key_j,exc,phase,N_int)
  call decode_exc(exc,2,h1,p1,h2,p2,s1,s2)

  if(.not.adjoint_tc_h)then ! Usual transcorrelated Hamiltonian 
   ! opposite spin two-body 
   herf    = get_mo_two_e_integral_erf(p1,p2,h1,h2,mo_integrals_erf_map)   
   heff    = mo_two_e_integral_eff_pot(p1,p2,h1,h2) 
   hderiv  = mo_two_e_eff_dr12_pot_array_physicist(p1,p2,h1,h2) 
   ! same spin two-body 
   if(s1.eq.s2)then
    herf   -= get_mo_two_e_integral_erf(p1,p2,h2,h1,mo_integrals_erf_map)   
    heff   -= mo_two_e_integral_eff_pot(p1,p2,h2,h1) 
    hderiv -= 0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(p1,p2,h2,h1) & 
             +0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(p2,p1,h1,h2)
   endif
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!! ADJOINT OF THE TC HAMILTONIAN !!!!!!!!!!!!!!!!!!!!!!!
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!! ADJOINT OF THE TC HAMILTONIAN !!!!!!!!!!!!!!!!!!!!!!!
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!! ADJOINT OF THE TC HAMILTONIAN !!!!!!!!!!!!!!!!!!!!!!!
  else 
   ! opposite spin two-body 
   herf    = 2.d0 * get_two_e_integral(p1,p2,h1,h2,mo_integrals_map)   
   herf   += -get_mo_two_e_integral_erf(p1,p2,h1,h2,mo_integrals_erf_map)   
   heff    = mo_two_e_integral_eff_pot(p1,p2,h1,h2) 
   hderiv  = -mo_two_e_eff_dr12_pot_array_physicist(p1,p2,h1,h2) 
   ! same spin two-body 
   if(s1.eq.s2)then
    herf   -= 2.d0 * get_two_e_integral(p1,p2,h2,h1,mo_integrals_map)   
    herf   += get_mo_two_e_integral_erf(p1,p2,h2,h1,mo_integrals_erf_map)   
    heff   -= mo_two_e_integral_eff_pot(p1,p2,h2,h1) 
    hderiv += 0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(p1,p2,h2,h1) & 
             +0.5d0 * mo_two_e_eff_dr12_pot_array_physicist(p2,p1,h1,h2)
   endif
  endif
  if(three_body_h_tc)then
   if(double_3_body_tc)then
    ! alpha/alpha/beta threee-body 
    if(Ne(1)+Ne(2).ge.3)then
     if(s1.eq.s2.and.s2.eq.1)then ! double alpha 
      do k = 1, Ne(2)
       kk = occ(k,2)
       hthree += three_body_5_index(kk,h1,h2,p1,p2)
      enddo
     else if(s1.eq.s2.and.s2.eq.2)then ! double beta 
      do k = 1, Ne(1)
       kk = occ(k,1)
       hthree += three_body_5_index(kk,h1,h2,p1,p2)
      enddo
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
   endif
  endif
  herf   *= phase
  heff   *= phase
  hderiv *= phase
  hthree  *= phase
  htot = herf + heff + hderiv + hthree
 end
