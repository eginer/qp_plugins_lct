subroutine comp_mo_int_mu_of_r_jl_old(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC
  
  integer, intent(in)            :: j,l ! r1
  integer,intent(out)            :: n_integrals
  integer(key_kind),intent(out)  :: buffer_i(mo_num*mo_num)
  real(integral_kind),intent(out) :: buffer_value(mo_num*mo_num)

  integer                        :: i,k
  double precision               :: cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0,get_two_e_integral
  double precision               :: thr
  integer                        :: kk, m, j1, i1

  thr = mo_integrals_threshold
  
  n_integrals = 0
  
  j1 = j+ishft(l*l-l,-1)
  do k = 1, mo_num           ! r1
    i1 = ishft(k*k-k,-1)
    if (i1 > j1) then
      exit
    endif
    do i = 1, k
      i1 += 1
      if (i1 > j1) then
        exit
      endif
!      if (mo_overlap_abs(i,k)*mo_overlap_abs(j,l) < thr) then
!        cycle
!      endif
      !DIR$ FORCEINLINE
      integral = get_two_e_integral(i,j,k,l,mo_integrals_map) ! i,k : r1    j,l : r2
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

subroutine comp_mo_int_mu_of_r_jl(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC
  
  integer, intent(in)            :: j,l
  integer,intent(out)            :: n_integrals
  integer(key_kind),intent(out)  :: buffer_i(mo_num*mo_num)
  real(integral_kind),intent(out) :: buffer_value(mo_num*mo_num)

  integer                        :: i,k
  double precision               :: cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  double precision               :: thr
  integer                        :: kk, m, j1, i1

  thr = mo_integrals_threshold
  
  n_integrals = 0
  
  double precision :: integrals_matrix(mo_num,mo_num)
  call get_mo_two_e_integrals_i1j1(j,l,mo_num,integrals_matrix,mo_integrals_map)
  j1 = j+ishft(l*l-l,-1)
  do k = 1, mo_num           ! r1
    i1 = ishft(k*k-k,-1)
    do i = 1, mo_num
      i1 += 1
      integral = integrals_matrix(k,i)  + 1.d0 ! i,k : r1    j,l : r2
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


