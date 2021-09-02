
subroutine compute_all_ijkl_for_jl_mu_of_r_int(j,l,ao_integrals)
 implicit none
 integer, intent(in) :: j,l
 double precision, intent(out) :: ao_integrals(ao_num,ao_num)
 provide d_dn2_e_cmd_basis aos_in_r_array 
 integer :: ipoint,istate,i,k
 
 double precision, allocatable :: ao_tmp(:,:)
 allocate(ao_tmp(ao_num, n_points_final_grid))
 istate = 1
 do ipoint = 1, n_points_final_grid
  do k = 1, ao_num
   ao_tmp(k,ipoint) = phi_ij_eff_pot_in_r(ipoint,j,l) * aos_in_r_array(k,ipoint)
  enddo
 enddo
 
 ao_integrals = 0.d0
 call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                     &
                   ao_tmp,size(ao_tmp,1), &
                   aos_in_r_array,size(aos_in_r_array,1),1.d0,                                  &
                    ao_integrals,size(ao_integrals,1))
! do ipoint = 1, n_points_final_grid
!  do k = 1, ao_num
!   do i = 1, ao_num
!    ao_integrals(i,k) += ao_tmp(k,ipoint) * aos_in_r_array(i,ipoint)
!   enddo
!  enddo
! enddo
end

subroutine compute_ao_mu_of_r_integrals_jl(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC

  integer, intent(in)             :: j,l
  integer,intent(out)             :: n_integrals
  integer(key_kind),intent(out)   :: buffer_i(ao_num*ao_num)
  real(integral_kind),intent(out) :: buffer_value(ao_num*ao_num)
  logical, external               :: ao_two_e_integral_zero

  integer                         :: i,k
  double precision                :: ao_two_e_integral,cpu_1,cpu_2, wall_1, wall_2
  double precision                :: integral, wall_0
  double precision                :: thr
  integer                         :: kk, m, j1, i1
  double precision, allocatable :: ao_integrals(:,:)
  allocate( ao_integrals(ao_num,ao_num) )

  call compute_all_ijkl_for_jl_mu_of_r_int(j,l,ao_integrals)
  thr = ao_integrals_threshold

  n_integrals = 0

  j1 = j+shiftr(l*l-l,1)
  do k = 1, ao_num           ! r1
    i1 = shiftr(k*k-k,1)
    if (i1 > j1) then
      exit
    endif
    do i = 1, k
      i1 += 1
      if (i1 > j1) then
        exit
      endif
      if (ao_two_e_integral_zero(i,j,k,l)) then
        cycle
      endif
      !DIR$ FORCEINLINE
      integral = ao_two_e_integral(i,k,j,l)  ! i,k : r1    j,l : r2
      if (abs(integral) < thr) then
        cycle
      endif
      n_integrals += 1
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
      buffer_value(n_integrals) = integral + ao_integrals(i,k)
    enddo
  enddo

end
