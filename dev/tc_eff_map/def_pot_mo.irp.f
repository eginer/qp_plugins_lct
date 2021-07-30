

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
  double precision               :: get_two_e_integral

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
!     if (ao_two_e_integral_zero(i,j,k,l)) then
      if (.False.) then
        cycle
      endif
      !DIR$ FORCEINLINE
      if(  trim(mo_class(i)) == "Core" .or. trim(mo_class(j)) == "Core" .or. trim(mo_class(k)) == "Core" .or. trim(mo_class(l)) == "Core")then
       integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
      else
       integral = scalar_mu_r_pot_physicist_mo(i,j,k,l) & 
                + 0.5d0 * (deriv_mu_r_pot_physicist_mo(i,j,k,l) + deriv_mu_r_pot_physicist_mo(k,l,i,j))
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
