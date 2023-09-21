

subroutine compute_mo_integrals_tc_int_jl(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for MO integrals
  END_DOC

  integer, intent(in)            :: j,l
  integer,intent(out)            :: n_integrals
  integer(key_kind),intent(out)  :: buffer_i(ao_num*ao_num)
  real(integral_kind),intent(out) :: buffer_value(ao_num*ao_num)

  integer                        :: i,k
  double precision               :: cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  double precision               :: thr
  integer                        :: kk, m, j1, i1
  double precision               :: get_two_e_integral,get_mo_tc_sym_two_e_pot

  thr = ao_integrals_threshold

  n_integrals = 0

  j1 = j+ishft(l*l-l,-1)
  do k = 1, ao_num           ! r1
    i1 = ishft(k*k-k,-1)
    if (i1 > j1) then
      exit
    endif
    do i = 1, k
      i1 += 1
      if (i1 > j1) then
        exit
      endif
      !DIR$ FORCEINLINE
      if(  trim(mo_class(i)) == "Core" .or. trim(mo_class(j)) == "Core" .or. trim(mo_class(k)) == "Core" .or. trim(mo_class(l)) == "Core")then
       integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
      else
       integral = get_mo_tc_sym_two_e_pot(i,j,k,l,mo_tc_sym_two_e_pot_map)
      endif
      if (abs(integral) < thr) then
        cycle
      endif
      n_integrals += 1
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
      buffer_value(n_integrals) = integral
    enddo
  enddo

end

subroutine compute_mo_integrals_tc_int_jl_zero(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for MO integrals
  END_DOC

  integer, intent(in)            :: j,l
  integer,intent(out)            :: n_integrals
  integer(key_kind),intent(out)  :: buffer_i(ao_num*ao_num)
  real(integral_kind),intent(out) :: buffer_value(ao_num*ao_num)

  integer                        :: i,k
  double precision               :: cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  double precision               :: thr
  integer                        :: kk, m, j1, i1
  double precision               :: get_two_e_integral,get_mo_tc_sym_two_e_pot

  thr = ao_integrals_threshold

  n_integrals = 0

  j1 = j+ishft(l*l-l,-1)
  do k = 1, ao_num           ! r1
    i1 = ishft(k*k-k,-1)
    if (i1 > j1) then
      exit
    endif
    do i = 1, k
      i1 += 1
      if (i1 > j1) then
        exit
      endif
      !DIR$ FORCEINLINE
      if(  trim(mo_class(i)) == "Core" .or. trim(mo_class(j)) == "Core" .or. trim(mo_class(k)) == "Core" .or. trim(mo_class(l)) == "Core")then
       integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
      else
       integral = get_mo_tc_sym_two_e_pot(i,j,k,l,mo_tc_sym_two_e_pot_map)  & 
           + 0.5d0 * (mo_non_hermit_term(i,j,k,l) + mo_non_hermit_term(k,l,i,j))
      endif
      if (abs(integral) < thr) then
        cycle
      endif
      n_integrals += 1
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
      buffer_value(n_integrals) = integral
    enddo
  enddo

end
