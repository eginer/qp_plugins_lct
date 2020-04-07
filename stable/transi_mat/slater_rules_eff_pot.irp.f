subroutine i_H_j_eff_pot(key_i,key_j,pot_a,pot_b,m,Nint,hij)
 implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
 integer, intent(in)            :: Nint,m
 integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
 double precision, intent(in)   :: pot_a(m,m),pot_b(m,m)
 double precision, intent(out)  :: hij
 integer :: degree
 call get_excitation_degree(key_i,key_j,degree,Nint)
 hij = 0.d0
 if(degree.gt.1)then
  return
 else if (degree == 0)then
  call i_H_j_eff_pot_diag(key_i,pot_a,pot_b,m,Nint,hij) 
 else
  call i_H_j_eff_pot_off_diag(key_i,key_j,pot_a,pot_b,m,Nint,hij) 
 endif

end

subroutine i_H_j_eff_pot_diag(key_i,pot_a,pot_b,m,Nint,hij)
 implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
 integer, intent(in)            :: Nint,m
 integer(bit_kind), intent(in)  :: key_i(Nint,2)
 double precision, intent(in)   :: pot_a(m,m),pot_b(m,m)
 double precision, intent(out)  :: hij
 integer, allocatable           :: occ(:,:)
 integer                        :: n_occ_ab(2),i,j
 allocate(occ(Nint*bit_kind_size,2))
 call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
 hij = 0.d0
 do i = 1, n_occ_ab(1) ! browsing the alpha electrons
  j = occ(i,1)
  hij += pot_a(j,j)
 enddo
 do i = 1, n_occ_ab(2) ! browsing the beta electrons
  j = occ(i,2)
  hij += pot_b(j,j)
 enddo

end

subroutine i_H_j_eff_pot_off_diag(key_i,key_j,pot_a,pot_b,m,Nint,hij)
 implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
 integer, intent(in)            :: Nint,m
 integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
 double precision, intent(in)   :: pot_a(m,m),pot_b(m,m)
 double precision               :: hij
 integer                        :: exc(0:2,2,2)
 double precision :: phase
 integer :: i,j
 call get_single_excitation(key_i,key_j,exc,phase,Nint)
 if (exc(0,1,1) == 1) then
   ! Single alpha
   i = exc(1,1,1)
   j = exc(1,2,1)
   hij = phase * pot_a(i,j)
 else
   ! Single beta
   i = exc(1,1,2)
   j = exc(1,2,2)
   hij = phase * pot_b(i,j)
 endif

end
