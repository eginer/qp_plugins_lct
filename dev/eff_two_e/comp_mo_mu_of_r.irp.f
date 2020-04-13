subroutine comp_mo_int_mu_of_r_jl(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  j,l ==> compute all <ij| 1/r12 |kl> --> <ij|  weeb(r12)  |kl> 
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
  
  double precision, allocatable :: integrals_matrix(:,:),eff_int_mat(:,:)
  allocate(integrals_matrix(mo_num,mo_num),eff_int_mat(mo_num,mo_num))
  call get_mo_two_e_integrals_i1j1(j,l,mo_num,integrals_matrix,mo_integrals_map)
  call compute_all_ijkl_for_jl_mu_of_r_int(j,l,eff_int_mat)
  j1 = j+ishft(l*l-l,-1)
  do k = 1, mo_num           ! r1
    i1 = ishft(k*k-k,-1)
    do i = 1, mo_num
      i1 += 1
      integral = integrals_matrix(k,i)  + + eff_int_mat(k,i) ! i,k : r1    j,l : r2
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

subroutine compute_all_ijkl_for_jl_mu_of_r_int(j,l,mo_integrals)
 implicit none
 integer, intent(in) :: j,l
 double precision, intent(out) :: mo_integrals(mo_num,mo_num)
 provide d_dn2_e_cmd_sr_pbe_n2 mos_in_r_array 
 integer :: ipoint,istate,i,k
 
 double precision, allocatable :: mo_tmp(:,:)
 allocate(mo_tmp(mo_num, n_points_final_grid))
 istate = 1
 do ipoint = 1, n_points_final_grid
  do k = 1, mo_num
   mo_tmp(k,ipoint) = phi_ij_eff_pot_in_r(ipoint,j,l) * mos_in_r_array(k,ipoint)
  enddo
 enddo
 
 mo_integrals = 0.d0
 call dgemm('N','T',mo_num,mo_num,n_points_final_grid,1.d0,                                     &
                   mo_tmp,size(mo_tmp,1), &
                   mos_in_r_array,size(mos_in_r_array,1),1.d0,                                  &
                    mo_integrals,size(mo_integrals,1))
! do ipoint = 1, n_points_final_grid
!  do k = 1, mo_num
!   do i = 1, mo_num
!    mo_integrals(i,k) += mo_tmp(k,ipoint) * mos_in_r_array(i,ipoint)
!   enddo
!  enddo
! enddo
end
