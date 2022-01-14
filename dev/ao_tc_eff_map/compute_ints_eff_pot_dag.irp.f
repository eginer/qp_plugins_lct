
subroutine compute_ao_integrals_gauss_dag_jl(j, l, n_integrals, buffer_i, buffer_value)

  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC

  use map_module

  implicit none

  integer, intent(in)             :: j, l
  integer,            intent(out) :: n_integrals
  integer(key_kind),  intent(out) :: buffer_i(ao_num*ao_num)
  real(integral_kind),intent(out) :: buffer_value(ao_num*ao_num)

  integer                         :: i, k, j1, i1
  double precision                :: thr, integral

  double precision                :: ao_tc_sym_two_e_pot_dag

  thr = ao_integrals_threshold

  n_integrals = 0

  j1 = j+ishft(l*l-l,-1)
  do k = 1, ao_num           ! r1
    i1 = ishft(k*k-k,-1)
    if(i1 > j1) exit

    do i = 1, k
      i1 += 1
      if(i1 > j1) exit

      if(ao_two_e_integral_erf_schwartz(i,k)*ao_two_e_integral_erf_schwartz(j,l) < thr ) cycle

      !DIR$ FORCEINLINE
      integral = ao_tc_sym_two_e_pot_dag(i,k,j,l)  ! i,k : r1    j,l : r2
      if(abs(integral) < thr) cycle
      n_integrals += 1
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
      buffer_value(n_integrals) = integral

    enddo
  enddo

end subroutine compute_ao_integrals_gauss_dag_jl



subroutine compute_ao_tc_sym_two_e_pot_dag_jl(j, l, n_integrals, buffer_i, buffer_value)

  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC

  use map_module

  implicit none

  integer, intent(in)              :: j, l
  integer,             intent(out) :: n_integrals
  integer(key_kind),   intent(out) :: buffer_i(ao_num*ao_num)
  real(integral_kind), intent(out) :: buffer_value(ao_num*ao_num)

  integer                        :: i, k, j1, i1
  double precision               :: integral, integral_pot, integral_erf, thr

  double precision               :: ao_tc_sym_two_e_pot_dag, ao_two_e_integral_erf

  thr = ao_integrals_threshold

  n_integrals = 0

  j1 = j+ishft(l*l-l,-1)
  do k = 1, ao_num           ! r1
    i1 = ishft(k*k-k,-1)
    if(i1 > j1) exit

    do i = 1, k
      i1 += 1
      if(i1 > j1) exit

      if( ao_two_e_integral_erf_schwartz(i,k)*ao_two_e_integral_erf_schwartz(j,l) < thr ) cycle

      !DIR$ FORCEINLINE
      integral_pot = ao_tc_sym_two_e_pot_dag(i,k,j,l)  ! i,k : r1    j,l : r2
      integral_erf = ao_two_e_integral_erf(i,k,j,l)
      integral = integral_erf - integral_pot
      if( abs(integral) < thr ) cycle
      n_integrals += 1
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i, j, k, l, buffer_i(n_integrals))
      buffer_value(n_integrals) = integral

    enddo
  enddo

end subroutine compute_ao_tc_sym_two_e_pot_dag_jl
