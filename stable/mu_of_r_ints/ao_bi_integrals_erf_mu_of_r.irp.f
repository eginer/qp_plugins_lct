subroutine compute_ao_integrals_erf_mu_of_r_jl_old(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC
  
  integer, intent(in)            :: j,l
  integer,intent(out)            :: n_integrals
  integer(key_kind),intent(out)  :: buffer_i(ao_num*ao_num)
  real(integral_kind),intent(out) :: buffer_value(ao_num*ao_num)

  integer                        :: i,k
  double precision               :: ao_bielec_integral_erf,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  double precision               :: thr,erf_mu_of_r_ao
  integer                        :: kk, m, j1, i1

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
      if (ao_overlap_abs(i,k)*ao_overlap_abs(j,l) < thr) then
        cycle
      endif
!     if (ao_bielec_integral_erf_schwartz(i,k)*ao_bielec_integral_erf_schwartz(j,l) < thr ) then
!       cycle
!     endif
      !DIR$ FORCEINLINE
      integral = erf_mu_of_r_ao(i,k,j,l)  ! i,k : r1    j,l : r2
      if (abs(integral) < thr) then
        cycle
      endif
      n_integrals += 1
      !DIR$ FORCEINLINE
      call bielec_integrals_index(i,j,k,l,buffer_i(n_integrals))
      buffer_value(n_integrals) = integral
    enddo
  enddo
    
end

subroutine compute_ao_integrals_erf_mu_of_r_jl(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC
  
  integer, intent(in)            :: j,l
  integer,intent(out)            :: n_integrals
  integer(key_kind),intent(out)  :: buffer_i(ao_num*ao_num)
  real(integral_kind),intent(out) :: buffer_value(ao_num*ao_num)

  integer                        :: i,k
  double precision               :: ao_bielec_integral_erf,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  double precision               :: thr,erf_mu_of_r_ao
  integer                        :: kk, m, j1, i1

  thr = ao_integrals_threshold
  
  n_integrals = 0
  
  double precision :: integrals_matrix(ao_num,ao_num)
  call give_all_erf_mu_of_r_kl(j,l,integrals_matrix)
  j1 = j+ishft(l*l-l,-1)
  do k = 1, ao_num           ! r1
    i1 = ishft(k*k-k,-1)
    do i = 1, ao_num
      i1 += 1
      integral = integrals_matrix(k,i)  ! i,k : r1    j,l : r2
      if (abs(integral) < thr) then
        cycle
      endif
      n_integrals += 1
      !DIR$ FORCEINLINE
      call bielec_integrals_index_no_sym(i,j,k,l,ao_num,buffer_i(n_integrals))
      buffer_value(n_integrals) = integral
    enddo
  enddo
    
end


