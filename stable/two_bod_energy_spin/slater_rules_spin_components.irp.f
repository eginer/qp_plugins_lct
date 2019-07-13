subroutine i_wee_j_single_alpha_beta(key1,key2,Nint,hij)
 use bitmasks
 implicit none
 BEGIN_DOC
 ! computes the <key1| w_ee(alpha/beta) |key2> where key1 and key2 are connected by a single excitation
 END_DOC
 integer, intent(in)            :: Nint
 integer(bit_kind), intent(in)  :: key1(Nint,2), key2(Nint,2)
 double precision, intent(out) :: hij
 integer :: i,j,h,p
 integer :: n_occ_ab
 integer           :: exc(0:2,2,2)
 double precision  :: phase
 integer           :: occ(Nint*bit_kind_size,2)
 call get_single_excitation(key1,key2,exc,phase,Nint)
 call bitstring_to_list_ab(key1, occ, n_occ_ab, Nint)
 hij = 0.d0
 if (exc(0,1,1) == 1) then
  ! Mono alpha
  h = exc(1,1,1)
  p = exc(1,2,1)
  do j = 1, elec_beta_num
   hij += big_array_coulomb_integrals(occ(j,2),p,h)
  enddo
 else 
  ! Mono beta
  h = exc(1,1,2)
  p = exc(1,2,2)
! print*,'h,p',h,p
  do j = 1, elec_alpha_num
   hij += big_array_coulomb_integrals(occ(j,1),p,h)
  enddo
 endif
 hij = hij * phase
end 

subroutine i_wee_j_single_alpha_alpha(key1,key2,Nint,hij)
 use bitmasks
 implicit none
 BEGIN_DOC
 ! computes the <key1| w_ee(alpha/alpha) |key2> where key1 and key2 are connected by a single excitation
 END_DOC
 integer, intent(in)            :: Nint
 integer(bit_kind), intent(in)  :: key1(Nint,2), key2(Nint,2)
 double precision, intent(out) :: hij
 integer :: i,j,h,p
 integer :: n_occ_ab
 integer           :: exc(0:2,2,2)
 double precision  :: phase
 integer           :: occ(Nint*bit_kind_size,2)
 call get_single_excitation(key1,key2,exc,phase,Nint)
 call bitstring_to_list_ab(key1, occ, n_occ_ab, Nint)
 hij = 0.d0
 if (exc(0,1,1) == 1) then
  ! Mono alpha
  h = exc(1,1,1)
  p = exc(1,2,1)
  do j = 1, elec_alpha_num
   hij += big_array_exchange_integrals(occ(j,1),p,h) 
  enddo
 else 
  ! Mono beta
  hij = 0.d0
 endif
 hij = hij * phase
end 

subroutine i_wee_j_single_beta_beta(key1,key2,Nint,hij)
 use bitmasks
 implicit none
 BEGIN_DOC
 ! computes the <key1| w_ee(beta/beta) |key2> where key1 and key2 are connected by a single excitation
 END_DOC
 integer, intent(in)            :: Nint
 integer(bit_kind), intent(in)  :: key1(Nint,2), key2(Nint,2)
 double precision, intent(out) :: hij
 integer :: i,j,h,p
 integer :: n_occ_ab
 integer           :: exc(0:2,2,2)
 double precision  :: phase
 integer           :: occ(Nint*bit_kind_size,2)
 call get_single_excitation(key1,key2,exc,phase,Nint)
 call bitstring_to_list_ab(key1, occ, n_occ_ab, Nint)
 hij = 0.d0
 if (exc(0,1,1) == 1) then
  ! Mono alpha
  hij = 0.d0
 else 
  ! Mono beta
  h = exc(1,1,2)
  p = exc(1,2,2)
  do j = 1, elec_beta_num
   hij += big_array_exchange_integrals(occ(j,2),p,h)
  enddo
 endif
 hij = hij * phase
end 


subroutine i_wee_i_diag_alpha_beta(key1,Nint,hij)
 use bitmasks
 implicit none
 BEGIN_DOC
 ! computes the <key1| w_ee(alpha/beta) |key1> 
 END_DOC
 integer, intent(in)            :: Nint
 integer(bit_kind), intent(in)  :: key1(Nint,2)
 double precision, intent(out) :: hij
 integer :: i,j,k,l
 integer :: n_occ_ab(2)
 integer           :: occ(Nint*bit_kind_size,2)
 call bitstring_to_list_ab(key1, occ, n_occ_ab, Nint)
 hij = 0.d0
 do i = 1, n_occ_ab(1)
  k = occ(i,1)
  do j = 1, n_occ_ab(2)
   l = occ(j,2)
   hij += mo_two_e_integrals_jj(k,l)
  enddo
 enddo 
end

subroutine i_wee_i_diag_alpha_alpha(key1,Nint,hij)
 use bitmasks
 implicit none
 BEGIN_DOC
 ! computes the <key1| w_ee(alpha/alpha) |key1> 
 END_DOC
 integer, intent(in)            :: Nint
 integer(bit_kind), intent(in)  :: key1(Nint,2)
 double precision, intent(out) :: hij
 integer :: i,j,k,l
 integer :: n_occ_ab
 integer           :: occ(Nint*bit_kind_size,2)
 call bitstring_to_list_ab(key1, occ, n_occ_ab, Nint)
 hij = 0.d0
 do i = 1, elec_alpha_num
  k = occ(i,1)
  do j = 1, elec_alpha_num
   l = occ(j,1)
   hij += mo_two_e_integrals_jj_anti(k,l)
  enddo
 enddo 
 hij = 0.5d0 * hij

end

subroutine i_wee_i_diag_beta_beta(key1,Nint,hij)
 use bitmasks
 implicit none
 BEGIN_DOC
 ! computes the <key1| w_ee(beta/beta) |key1> 
 END_DOC
 integer, intent(in)            :: Nint
 integer(bit_kind), intent(in)  :: key1(Nint,2)
 double precision, intent(out) :: hij
 integer :: i,j,k,l
 integer :: n_occ_ab
 integer           :: occ(Nint*bit_kind_size,2)
 call bitstring_to_list_ab(key1, occ, n_occ_ab, Nint)
 hij = 0.d0
 do i = 1, elec_beta_num
  k = occ(i,2)
  do j = 1, elec_beta_num
   l = occ(j,2)
   hij += mo_two_e_integrals_jj_anti(k,l)
  enddo
 enddo 
 hij = 0.5d0 * hij

end
