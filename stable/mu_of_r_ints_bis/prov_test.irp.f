BEGIN_PROVIDER [ double precision, mo_mu_of_r_two_ints_phys, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 integer :: i,j,k,l
 double precision :: comp_mo_mu_of_r_two_ints_phys,integral
 integral = comp_mo_mu_of_r_two_ints_phys(1,1,1,1)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,k,l) & 
 !$OMP SHARED (mo_num,mo_mu_of_r_two_ints_phys)
 !$OMP DO SCHEDULE (dynamic)
 do l = 1, mo_num
  do k = 1, mo_num
   do j = 1, mo_num
    do i = 1, mo_num
     mo_mu_of_r_two_ints_phys(i,j,k,l) = comp_mo_mu_of_r_two_ints_phys(i,j,k,l)
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
END_PROVIDER 


double precision function comp_mo_mu_of_r_two_ints_phys(i,j,k,l)
 implicit none
 integer, intent(in) :: i,j,k,l
 integer :: ipoint
 double precision :: weight, contrib_jl, contrib_ik
 comp_mo_mu_of_r_two_ints_phys = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  contrib_ik = erf_mu_of_r_mo_mo_transp(ipoint,i,k) * mos_in_r_array_transp(ipoint,j) * mos_in_r_array_transp(ipoint,l)  
  contrib_jl = erf_mu_of_r_mo_mo_transp(ipoint,j,l) * mos_in_r_array_transp(ipoint,i) * mos_in_r_array_transp(ipoint,k) 
  comp_mo_mu_of_r_two_ints_phys +=  (contrib_jl + contrib_ik) * 0.5d0 * weight
 enddo
end


BEGIN_PROVIDER [ double precision, erf_mu_of_r_mo_mo, (mo_num, mo_num,n_points_final_grid)]
 implicit none
 integer :: ipoint
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (ipoint) & 
  !$OMP SHARED (ao_num,n_points_final_grid,erf_mu_of_r_mo_mo,erf_mu_of_r_ao_ao, mo_num)
  !$OMP DO SCHEDULE (dynamic)
  do ipoint = 1, n_points_final_grid
   call ao_to_mo(erf_mu_of_r_ao_ao(1,1,ipoint),ao_num,erf_mu_of_r_mo_mo(1,1,ipoint),mo_num)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

END_PROVIDER 

BEGIN_PROVIDER [ double precision, erf_mu_of_r_mo_mo_transp, (n_points_final_grid, mo_num, mo_num) ]
 implicit none
 integer :: i,j,ipoint
 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, ao_num
    erf_mu_of_r_mo_mo_transp(ipoint,j,i) = erf_mu_of_r_mo_mo(j,i,ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_mu_of_r_two_ints_phys, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 integer :: i,j,k,l
 double precision :: comp_ao_mu_of_r_two_ints_phys,integral
 integral = comp_ao_mu_of_r_two_ints_phys(1,1,1,1)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,k,l) & 
 !$OMP SHARED (ao_num,ao_mu_of_r_two_ints_phys)
 !$OMP DO SCHEDULE (dynamic)
 do l = 1, ao_num
  do k = 1, ao_num
   do j = 1, ao_num
    do i = 1, ao_num
     ao_mu_of_r_two_ints_phys(i,j,k,l) = comp_ao_mu_of_r_two_ints_phys(i,j,k,l)
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
END_PROVIDER 


double precision function comp_ao_mu_of_r_two_ints_phys(i,j,k,l)
 implicit none
 integer, intent(in) :: i,j,k,l
 integer :: ipoint
 double precision :: weight, contrib_jl, contrib_ik
 comp_ao_mu_of_r_two_ints_phys = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  contrib_ik = erf_mu_of_r_ao_ao_transp(ipoint,i,k) * aos_in_r_array_transp(ipoint,j) * aos_in_r_array_transp(ipoint,l)  
  contrib_jl = erf_mu_of_r_ao_ao_transp(ipoint,j,l) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,k) 
  comp_ao_mu_of_r_two_ints_phys +=  (contrib_jl + contrib_ik) * 0.5d0 * weight
 enddo
end

