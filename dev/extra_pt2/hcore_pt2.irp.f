
subroutine i_H_j_hcore(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|h_{core}|j \rangle$ where $i$ and $j$ are determinants.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij

  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: get_two_e_integral
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem, phase
  integer                        :: n_occ_ab(2)
  PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)


  hij = 0.d0
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  integer :: spin
  select case (degree)
    case (2)
      hij = 0.d0
    case (1)
      call get_single_excitation(key_i,key_j,exc,phase,Nint)
      !DIR$ FORCEINLINE
      call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        spin = 1
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        spin = 2
      endif
      hij = mo_one_e_integrals(m,p)

    case (0)
      hij = 0.d0
      do m = 1, n_occ_ab(1)
       hij += mo_one_e_integrals(occ(m,1),occ(m,1))
      enddo
      do m = 1, n_occ_ab(2)
       hij += mo_one_e_integrals(occ(m,2),occ(m,2))
      enddo
  end select
end


