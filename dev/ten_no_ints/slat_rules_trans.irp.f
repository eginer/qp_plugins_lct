subroutine diag_htilde_ten_no_mat(key_i,hmono,herf,heff,hderiv,htot)
  use bitmasks
  BEGIN_DOC
!  diagonal element of htilde 
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv,htot
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: ne(2),i,j,ii,jj,ispin,jspin
  double precision :: get_two_e_integral 
  provide mo_ten_no_eff_sq_lpl_pot_physicist mo_ten_no_dr12_pot_physicist
  provide mo_two_e_integrals_in_map mo_integrals_map

  call bitstring_to_list_ab(key_i,occ,Ne,N_int)

  hmono = 0.d0
  do ispin = 1, 2 
   do i = 1, Ne(ispin) ! 
    ii = occ(i,ispin) 
    hmono += mo_one_e_integrals(ii,ii)
   enddo
  enddo

  herf  = 0.d0
  heff  = 0.d0
  hderiv= 0.d0

  ispin = 1
  jspin = 2 
  do i = 1, Ne(ispin) 
   ii = occ(i,ispin) 
   do j = 1, Ne(jspin)
    jj = occ(j,jspin) 
    herf   += get_two_e_integral(ii,jj,ii,jj,mo_integrals_map)
    heff   += mo_ten_no_eff_sq_lpl_pot_physicist(ii,jj,ii,jj) 
    hderiv += mo_ten_no_dr12_pot_physicist(ii,jj,ii,jj) 
   enddo
  enddo
  htot = hmono + herf + heff + hderiv

end

subroutine single_htilde_ten_no_mat(key_j,key_i,hmono,herf,heff,hderiv,htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> for single excitation  
!!
!! WARNING !!
! 
! Non hermitian !!
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv,htot
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: ne(2),i,j,ii,jj,ispin,jspin
  integer                        :: degree,exc(0:2,2,2)
  integer                        :: h1, p1, h2, p2, s1, s2
  double precision :: get_two_e_integral,mo_two_e_integral_eff_pot,phase
  integer :: other_spin(2)
  provide mo_ten_no_eff_sq_lpl_pot_physicist mo_ten_no_dr12_pot_physicist
  provide mo_two_e_integrals_in_map mo_integrals_map

  other_spin(1) = 2
  other_spin(2) = 1

  call get_excitation_degree(key_i,key_j,degree,N_int)
  hmono = 0.d0
  herf = 0.d0
  heff = 0.d0
  hderiv = 0.d0
  if(degree.ne.1)then
   return
  endif
  call bitstring_to_list_ab(key_i,occ,Ne,N_int)
  call get_single_excitation(key_i,key_j,exc,phase,N_int)
  call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)

  hmono = mo_one_e_integrals(h1,p1) * phase
  
  ! Coulomb 
  ispin = other_spin(s1)
  do i = 1, Ne(ispin)
   ii = occ(i,ispin) 
   herf   += get_two_e_integral(ii,p1,ii,h1,mo_integrals_map)
   heff   += mo_ten_no_eff_sq_lpl_pot_physicist(ii,p1,ii,h1) 
   hderiv += mo_ten_no_dr12_pot_physicist(ii,p1,ii,h1) 
  enddo
  herf    *= phase
  heff    *= phase
  hderiv  *= phase
  htot = hmono + herf + heff + hderiv
end

subroutine double_htilde_ten_no_mat(key_j,key_i,hmono,herf,heff,hderiv,htot)
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
  double precision, intent(out)  :: hmono,herf,heff,hderiv,htot
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: ne(2),i,j,ii,jj,ispin,jspin
  integer                        :: degree,exc(0:2,2,2)
  integer                        :: h1, p1, h2, p2, s1, s2
  double precision :: get_two_e_integral,mo_two_e_integral_eff_pot,phase
  integer :: other_spin(2)
  provide mo_ten_no_eff_sq_lpl_pot_physicist mo_ten_no_dr12_pot_physicist
  provide mo_two_e_integrals_in_map mo_integrals_map
  other_spin(1) = 2
  other_spin(2) = 1

  call get_excitation_degree(key_i,key_j,degree,N_int)
  hmono = 0.d0
  herf = 0.d0
  heff = 0.d0
  hderiv = 0.d0
  if(degree.ne.2)then
   return
  endif
  call get_double_excitation(key_i,key_j,exc,phase,N_int)
  call decode_exc(exc,2,h1,p1,h2,p2,s1,s2)
  herf    = get_two_e_integral(p1,p2,h1,h2,mo_integrals_map)   
  heff    = mo_ten_no_eff_sq_lpl_pot_physicist(p1,p2,h1,h2) 
  hderiv  = mo_ten_no_dr12_pot_physicist(p1,p2,h1,h2) 
  herf   *= phase
  heff   *= phase
  hderiv *= phase
  htot = herf + heff + hderiv
 end


subroutine htilde_ten_no_mat(key_j,key_i,hmono,herf,heff,hderiv,htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> 
!!
!! WARNING !!
! 
! Non hermitian !!
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv,htot
  integer                        :: degree
   call get_excitation_degree(key_j,key_i,degree,N_int)
   hmono = 0.d0
   herf = 0.d0
   heff = 0.d0
   hderiv = 0.d0
   if(degree.gt.2)then
    return
   else if(degree == 2)then
    call double_htilde_ten_no_mat(key_j,key_i,hmono,herf,heff,hderiv,htot)
   else if(degree == 1)then
    call single_htilde_ten_no_mat(key_j,key_i,hmono,herf,heff,hderiv,htot)
   else if(degree == 0)then
    call diag_htilde_ten_no_mat(key_i,hmono,herf,heff,hderiv,htot)
   endif
   
end
