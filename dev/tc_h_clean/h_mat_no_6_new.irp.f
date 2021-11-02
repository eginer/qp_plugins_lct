!!!!!! 

subroutine diag_htilde_mu_mat_scal_map(key_i,hmono,heff,hderiv,htot)
  use bitmasks
  BEGIN_DOC
!  diagonal element of htilde ONLY FOR ONE- AND TWO-BODY TERMS 
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_i(N_int,2)
  double precision, intent(out)  :: hmono,heff,hderiv,htot
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  double precision :: get_mo_two_e_integral_tc_int
  PROVIDE mo_two_e_integrals_tc_int_in_map mo_non_hermit_term
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
  heff  = 0.d0
  hderiv= 0.d0
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


   ! alpha/beta two-body
   ispin = 1
   jspin = 2 
   do i = 1, Ne(ispin) ! electron 1 (so it can be associated to mu(r1))
    ii = occ(i,ispin) 
    do j = 1, Ne(jspin) ! electron 2 
     jj = occ(j,jspin) 
     heff   += get_mo_two_e_integral_tc_int(jj,ii,jj,ii,mo_integrals_tc_int_map) 
     hderiv += mo_non_hermit_term(jj,ii,jj,ii) 
    enddo
   enddo
 
   ! alpha/alpha two-body
   do i = 1, Ne(ispin)
    ii = occ(i,ispin) 
    do j = i+1, Ne(ispin)
!    do j = 1, Ne(ispin)
     jj = occ(j,ispin) 
     heff += get_mo_two_e_integral_tc_int(ii,jj,ii,jj,mo_integrals_tc_int_map)  & 
            -get_mo_two_e_integral_tc_int(jj,ii,ii,jj,mo_integrals_tc_int_map)
     hderiv += mo_non_hermit_term(ii,jj,ii,jj) - mo_non_hermit_term(ii,jj,jj,ii) 
    enddo
   enddo
 
   ! beta/beta two-body
   do i = 1, Ne(jspin)
    ii = occ(i,jspin) 
    do j = i+1, Ne(jspin)
     jj = occ(j,jspin) 
     heff += get_mo_two_e_integral_tc_int(ii,jj,ii,jj,mo_integrals_tc_int_map) & 
            -get_mo_two_e_integral_tc_int(jj,ii,ii,jj,mo_integrals_tc_int_map) 
     hderiv += mo_non_hermit_term(ii,jj,ii,jj) - mo_non_hermit_term(ii,jj,jj,ii) 
              
    enddo
   enddo
  htot = hmono + heff + hderiv 

end

subroutine double_htilde_mu_mat_scal_map(key_j,key_i,hmono,heff,hderiv,htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> for double excitation  ONLY FOR ONE- AND TWO-BODY TERMS 
!!
!! WARNING !!
! 
! Non hermitian !!
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  double precision, intent(out)  :: hmono,heff,hderiv,htot
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  integer                        :: degree,exc(0:2,2,2)
  integer                        :: h1, p1, h2, p2, s1, s2
  double precision :: get_mo_two_e_integral_tc_int,phase
  integer :: other_spin(2)
  other_spin(1) = 2
  other_spin(2) = 1
  PROVIDE mo_two_e_integrals_tc_int_in_map mo_non_hermit_term

  integer(bit_kind) :: key_i_core(N_int,2)
  call get_excitation_degree(key_i,key_j,degree,N_int)

  hmono = 0.d0
  heff  = 0.d0
  hderiv= 0.d0
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

  if(s1.ne.s2)then
   ! opposite spin two-body 
   if(s1==1)then
    heff    = get_mo_two_e_integral_tc_int(p2,p1,h2,h1,mo_integrals_tc_int_map) 
    if(double_normal_ord)then
     heff += normal_two_body(h1,p1,h2,p2)
    endif
    hderiv  = mo_non_hermit_term(p2,p1,h2,h1) 
   else
    heff    = get_mo_two_e_integral_tc_int(p1,p2,h1,h2,mo_integrals_tc_int_map) 
    if(double_normal_ord)then
     heff += normal_two_body(h1,p1,h2,p2)
    endif
    hderiv  = mo_non_hermit_term(p1,p2,h1,h2) 
   endif
  else
   ! same spin two-body 
   ! direct terms 
   heff    = get_mo_two_e_integral_tc_int(p2,p1,h2,h1,mo_integrals_tc_int_map) 
   if(double_normal_ord)then
    heff += normal_two_body(h1,p1,h2,p2)
   endif
   hderiv  = mo_non_hermit_term(p2,p1,h2,h1)  
   ! exchange terms 
   heff   -= get_mo_two_e_integral_tc_int(p1,p2,h2,h1,mo_integrals_tc_int_map) 
   if(double_normal_ord)then
    heff -= normal_two_body(h2,p1,h1,p2)
   endif
   hderiv -= mo_non_hermit_term(p1,p2,h2,h1) 
  endif
  heff   *= phase
  hderiv *= phase
  htot = heff + hderiv 
 end


subroutine single_htilde_mu_mat_scal_map(key_j,key_i,hmono,heff,hderiv,htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> for single excitation ONLY FOR ONE- AND TWO-BODY TERMS 
!!
!! WARNING !!
! 
! Non hermitian !!
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  double precision, intent(out)  :: hmono,heff,hderiv, htot
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  integer                        :: degree,exc(0:2,2,2)
  integer                        :: h1, p1, h2, p2, s1, s2
  double precision :: get_mo_two_e_integral_tc_int,phase
  double precision :: direct_int,exchange_int_12,exchange_int_23,exchange_int_13
  integer :: other_spin(2)
  PROVIDE mo_two_e_integrals_tc_int_in_map mo_non_hermit_term
  other_spin(1) = 2
  other_spin(2) = 1

  integer(bit_kind) :: key_j_core(N_int,2),key_i_core(N_int,2)

  hmono = 0.d0
  heff  = 0.d0
  hderiv= 0.d0
  htot = 0.d0
  call get_excitation_degree(key_i,key_j,degree,N_int)
  if(degree.ne.1)then
   return
  endif
  if(core_tc_op)then
   do i = 1, N_int
    key_i_core(i,1) = xor(key_i(i,1),core_bitmask(i,1))
    key_i_core(i,2) = xor(key_i(i,2),core_bitmask(i,2))
    key_j_core(i,1) = xor(key_j(i,1),core_bitmask(i,1))
    key_j_core(i,2) = xor(key_j(i,2),core_bitmask(i,2))
   enddo
   call bitstring_to_list_ab(key_i_core,occ,Ne,N_int)
  else
   call bitstring_to_list_ab(key_i,occ,Ne,N_int)
  endif


  call get_single_excitation(key_i,key_j,exc,phase,N_int)
  call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)

  hmono = mo_one_e_integrals(h1,p1) * phase
  if(core_tc_op)then
   hmono += phase * core_fock_operator(h1,p1)
  endif
  
   ! alpha/beta two-body 
   ispin = other_spin(s1)
   if(s1==1)then
    ! single alpha 
    do i = 1, Ne(ispin) ! electron 2 
     ii = occ(i,ispin) 
     heff   += get_mo_two_e_integral_tc_int(ii,p1,ii,h1,mo_integrals_tc_int_map) 
     hderiv += mo_non_hermit_term(ii,p1,ii,h1) 
    enddo
   else
    ! single beta 
    do i = 1, Ne(ispin) ! electron 1 
     ii = occ(i,ispin) 
     heff   += get_mo_two_e_integral_tc_int(p1,ii,h1,ii,mo_integrals_tc_int_map) 
     hderiv += mo_non_hermit_term(p1,ii,h1,ii) 
    enddo
   endif
!   ! same spin two-body 
   do i = 1, Ne(s1)
    ii = occ(i,s1) 
    ! (h1p1|ii ii) - (h1 ii| p1 ii)
    heff   += get_mo_two_e_integral_tc_int(ii,p1,ii,h1,mo_integrals_tc_int_map) & 
             -get_mo_two_e_integral_tc_int(p1,ii,ii,h1,mo_integrals_tc_int_map)
    hderiv += mo_non_hermit_term(ii,p1,ii,h1) - mo_non_hermit_term(p1,ii,ii,h1) 
   enddo
   
  heff    *= phase
  hderiv  *= phase
  htot = hmono + heff + hderiv 
end

