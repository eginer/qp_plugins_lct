
! For a derivation of the transcorrelated Hamiltonian matrix elements, see Luo-JCP-10 
! http://dx.doi.org/10.1063/1.3505037
subroutine direct_diag_htilde_mu_mat(key_i,hmono,herf,heff,hderiv,hthree,htot)
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
   double precision :: direct_int, exchange_int
   if(Ne(1)+Ne(2).ge.3)then
!!!  ! alpha/alpha/beta three-body
    do i = 1, Ne(1)
     ii = occ(i,1) 
     do j = i+1, Ne(1)
      jj = occ(j,1) 
      do k = 1, Ne(2)
       kk = occ(k,2) 
       call give_integrals_3_body(kk,jj,ii,kk,jj,ii,direct_int)
       call give_integrals_3_body(kk,jj,ii,kk,ii,jj,exchange_int)
       hthree += -(direct_int - exchange_int)
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
       call give_integrals_3_body(kk,jj,ii,kk,jj,ii,direct_int)
       call give_integrals_3_body(kk,ii,jj,kk,jj,ii,exchange_int)
       hthree += -(direct_int - exchange_int)
      enddo
     enddo
    enddo

    double precision :: exchange_int_12, exchange_int_13, exchange_int_23
    ! alpha/alpha/alpha three-body
    do i = 1, Ne(1)
     ii = occ(i,1)
     do j = i+1, Ne(1)
      jj = occ(j,1)
      do k = j+1, Ne(1)
       kk = occ(k,1)
       !               direct   :  1  2  3  1  2  3
       call give_integrals_3_body(ii,jj,kk,ii,jj,kk,direct_int)
       !        exchange 1<->2  :  2  1  3  1  2  3
       call give_integrals_3_body(jj,ii,kk,ii,jj,kk,exchange_int_12)
       !        exchange 1<->3  :  3  2  1  1  2  3
       call give_integrals_3_body(kk,jj,ii,ii,jj,kk,exchange_int_13)
       !        exchange 2<->3  :  1  3  2  1  2  3
       call give_integrals_3_body(ii,kk,jj,ii,jj,kk,exchange_int_23)
       hthree += -( direct_int - exchange_int_12 - exchange_int_13 - exchange_int_23 )
      enddo
     enddo
    enddo

    ! beta/beta/beta three-body
    do i = 1, Ne(2)
     ii = occ(i,2)
     do j = i+1, Ne(2)
      jj = occ(j,2)
      do k = j+1, Ne(2)
       kk = occ(k,2)
       !               direct   :  1  2  3  1  2  3
       call give_integrals_3_body(ii,jj,kk,ii,jj,kk,direct_int)
       !        exchange 1<->2  :  2  1  3  1  2  3
       call give_integrals_3_body(jj,ii,kk,ii,jj,kk,exchange_int_12)
       !        exchange 1<->3  :  3  2  1  1  2  3
       call give_integrals_3_body(kk,jj,ii,ii,jj,kk,exchange_int_13)
       !        exchange 2<->3  :  1  3  2  1  2  3
       call give_integrals_3_body(ii,kk,jj,ii,jj,kk,exchange_int_23)
       hthree += -( direct_int - exchange_int_12 - exchange_int_13 - exchange_int_23 )
      enddo
     enddo
    enddo

   endif
  endif
  htot = hmono + herf + heff + hderiv + hthree

end
