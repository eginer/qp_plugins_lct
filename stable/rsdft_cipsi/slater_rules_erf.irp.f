!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THIS FILE CONTAINS EVERYTHING YOU NEED TO COMPUTE THE LONG RANGE PART OF THE INTERACTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine i_H_j_erf(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|W_{ee}^{lr}|j> where i and j are determinants 
  ! and the W_{ee}^{lr} is the long range two-body interaction
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij
  
  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: get_mo_two_e_integral_erf
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem_erf, phase,phase_2
  integer                        :: n_occ_ab(2)
  PROVIDE mo_two_e_integrals_erf_in_map mo_integrals_erf_map int_erf_3_index_exc
  
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
      call get_double_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha, mono beta
        if(exc(1,1,1) == exc(1,2,2) )then
         hij = phase * int_erf_3_index_exc(exc(1,1,1),exc(1,1,2),exc(1,2,1))
        else if (exc(1,2,1) ==exc(1,1,2))then
         hij = phase * int_erf_3_index_exc(exc(1,2,1),exc(1,1,1),exc(1,2,2))
        else
         hij = phase*get_mo_two_e_integral_erf(                          &
             exc(1,1,1),                                              &
             exc(1,1,2),                                              &
             exc(1,2,1),                                              &
             exc(1,2,2) ,mo_integrals_erf_map)
        endif
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij = phase*(get_mo_two_e_integral_erf(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_erf_map) -                          &
            get_mo_two_e_integral_erf(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_erf_map) )
      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij = phase*(get_mo_two_e_integral_erf(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_erf_map) -                          &
            get_mo_two_e_integral_erf(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_erf_map) )
      endif
    case (1)
      call get_single_excitation(key_i,key_j,exc,phase,Nint)
      !DIR$ FORCEINLINE
      call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        spin = 1
        do i = 1, n_occ_ab(1)
         hij += -int_erf_3_index_exc(occ(i,1),m,p) + int_erf_3_index(occ(i,1),m,p)
        enddo
        do i = 1, n_occ_ab(2)
         hij += int_erf_3_index(occ(i,2),m,p)
        enddo
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        spin = 2
        do i = 1, n_occ_ab(2)
         hij += -int_erf_3_index_exc(occ(i,2),m,p) + int_erf_3_index(occ(i,2),m,p)
        enddo
        do i = 1, n_occ_ab(1)
         hij += int_erf_3_index(occ(i,1),m,p)
        enddo
      endif
      hij = hij * phase
    case (0)
      hij = diag_H_mat_elem_erf(key_i,Nint)
  end select
end


double precision function diag_H_mat_elem_erf(key_i,Nint)
 BEGIN_DOC 
! returns <i|W_{ee}^{lr}|i> where |i> is a determinant and 
! W_{ee}^{lr} is the two body long-range interaction
 END_DOC
 implicit none
 integer(bit_kind), intent(in) :: key_i(N_int,2)
 integer, intent(in)  :: Nint
 integer :: i,j
 integer                        :: occ(Nint*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
 diag_H_mat_elem_erf = 0.d0
 ! alpha - alpha
 do i = 1, n_occ_ab(1)
  do j = i+1, n_occ_ab(1)
   diag_H_mat_elem_erf += mo_two_e_int_erf_jj_anti(occ(i,1),occ(j,1))
  enddo
 enddo

 ! beta - beta 
 do i = 1, n_occ_ab(2)
  do j = i+1, n_occ_ab(2)
   diag_H_mat_elem_erf += mo_two_e_int_erf_jj_anti(occ(i,2),occ(j,2))
  enddo
 enddo

 ! alpha - beta 
 do i = 1, n_occ_ab(1)
  do j = 1, n_occ_ab(2)
   diag_H_mat_elem_erf += mo_two_e_int_erf_jj(occ(i,1),occ(j,2))
  enddo
 enddo
end
subroutine i_H_j_mono_spin_erf(key_i,key_j,Nint,spin,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants differing by a single excitation
  END_DOC
  integer, intent(in)            :: Nint, spin
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij
  
  integer                        :: exc(0:2,2)
  double precision               :: phase

  PROVIDE int_erf_3_index_exc mo_two_e_integrals_erf_in_map

  call i_H_j_erf(key_i,key_j,Nint,hij)
end



subroutine i_H_j_double_spin_erf(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants differing by a same-spin double excitation
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint), key_j(Nint)
  double precision, intent(out)  :: hij
  
  integer                        :: exc(0:2,2)
  double precision               :: phase
  double precision, external     :: get_mo_two_e_integral_erf

  PROVIDE int_erf_3_index_exc mo_two_e_integrals_erf_in_map

  call get_double_excitation_spin(key_i,key_j,exc,phase,Nint)
  hij = phase*(get_mo_two_e_integral_erf(                               &
      exc(1,1),                                                      &
      exc(2,1),                                                      &
      exc(1,2),                                                      &
      exc(2,2), mo_integrals_erf_map) -                                  &
      get_mo_two_e_integral_erf(                                        &
      exc(1,1),                                                      &
      exc(2,1),                                                      &
      exc(2,2),                                                      &
      exc(1,2), mo_integrals_erf_map) )
end

subroutine i_H_j_double_alpha_beta_erf(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants differing by an opposite-spin double excitation
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij
  
  integer                        :: exc(0:2,2,2)
  double precision               :: phase, phase2
  double precision, external     :: get_mo_two_e_integral_erf

  PROVIDE int_erf_3_index_exc mo_two_e_integrals_erf_in_map

  call get_single_excitation_spin(key_i(1,1),key_j(1,1),exc(0,1,1),phase,Nint)
  call get_single_excitation_spin(key_i(1,2),key_j(1,2),exc(0,1,2),phase2,Nint)
  phase = phase*phase2
  if (exc(1,1,1) == exc(1,2,2)) then
    hij = phase * int_erf_3_index_exc(exc(1,1,1),exc(1,1,2),exc(1,2,1))
  else if (exc(1,2,1) == exc(1,1,2)) then
    hij = phase * int_erf_3_index_exc(exc(1,2,1),exc(1,1,1),exc(1,2,2))
  else
    hij = phase*get_mo_two_e_integral_erf(                              &
        exc(1,1,1),                                                  &
        exc(1,1,2),                                                  &
        exc(1,2,1),                                                  &
        exc(1,2,2) ,mo_integrals_erf_map)
  endif
end


