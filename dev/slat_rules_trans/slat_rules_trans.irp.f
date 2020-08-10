subroutine diag_htilde_mat(key_i,hmono,herf,heff,hderiv)
  use bitmasks
  implicit none
  integer(bit_kind), intent(in)  :: key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: ne(2),i,j,ii,jj,ispin,jspin
  double precision :: get_mo_two_e_integral_erf,mo_two_e_integral_eff_pot

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

  do ispin = 1,2
   do i = 1, Ne(ispin) 
    ii = occ(i,ispin) 
    do j = i+1, Ne(ispin)
     jj = occ(j,ispin) 
     herf += get_mo_two_e_integral_erf(ii,jj,ii,jj,mo_integrals_erf_map) & 
            -get_mo_two_e_integral_erf(ii,jj,jj,ii,mo_integrals_erf_map)
     heff += mo_two_e_integral_eff_pot(ii,jj,ii,jj) - mo_two_e_integral_eff_pot(ii,jj,jj,ii) 
     hderiv += mo_two_e_eff_dr12_pot_array(ii,jj,ii,jj) - mo_two_e_eff_dr12_pot_array(ii,jj,jj,ii)
    enddo
   enddo
  enddo
  ispin = 1
  jspin = 2 
  do i = 1, Ne(ispin) 
   ii = occ(i,ispin) 
   do j = 1, Ne(jspin)
    jj = occ(j,jspin) 
    herf += get_mo_two_e_integral_erf(ii,jj,ii,jj,mo_integrals_erf_map)
    heff += mo_two_e_integral_eff_pot(ii,jj,ii,jj) 
    hderiv += mo_two_e_eff_dr12_pot_array(ii,jj,ii,jj) 
    print*,'mo_two_e_eff_dr12_pot_array(ii,jj,ii,jj',mo_two_e_eff_dr12_pot_array(ii,jj,ii,jj)
   enddo
  enddo

end
