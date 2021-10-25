program test
 implicit none
 my_grid_becke = .True.
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 call routine_read
end

subroutine routine_read
 implicit none
 character*(128) :: output
 integer :: i_unit_output
 integer :: getUnitAndOpen


 print*,'Reading the integrals '
 double precision, allocatable :: array_ints(:,:,:,:)
 allocate(array_ints(n_points_final_grid,3,mo_num,mo_num))
 output =trim(ezfio_filename)//'/x_w_ij_r'
 i_unit_output = getUnitAndOpen(output,'r')
 integer :: ipoint,m,i,j
 integer*8 ::itot,n_tot
 n_tot = n_points_final_grid * mo_num * mo_num * 3
 do itot = 1, n_tot
  read(i_unit_output,*)ipoint,m,j,i,array_ints(ipoint,m,j,i)
 enddo

 
 double precision :: accu
 accu = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do m = 1, 3
    do ipoint = 1, n_points_final_grid
    accu += dabs(array_ints(ipoint,m,j,i) - x_W_ij_erf_rk(ipoint,m,j,i))
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu
 print*,'accu/n ',accu/dble(n_tot)


 print*,'Reading the MOs '
 output =trim(ezfio_filename)//'/mos_in_r'
 i_unit_output = getUnitAndOpen(output,'r')
 double precision, allocatable :: array_mos(:,:)
 allocate(array_mos(n_points_final_grid,mo_num))
 n_tot = n_points_final_grid * mo_num 
 do itot = 1, n_tot
  read(i_unit_output,*)ipoint,i,array_mos(ipoint,i)
 enddo
 accu = 0.d0
 do i = 1, mo_num
  do ipoint = 1, n_points_final_grid
   accu += dabs(array_mos(ipoint,i)-mos_in_r_array_transp(ipoint,i) * sqrt_weight_at_r(ipoint))
  enddo
 enddo
 print*,'accu = ',accu


 print*,'Computing integrals '
 double precision :: integral,integral_ref
 integer :: k,l,n
 accu = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do m = 1, mo_num
    do k = 1, mo_num
     do l = 1, mo_num
      do n = j,j
    call compute_ijklmn(i,j,m,i,j,m,array_mos,array_ints,integral)
    call give_integrals_3_body(i,j,m,i,j,m,integral_ref)
    accu += dabs(integral_ref - integral)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 print*,'accu integrals = ',accu

end


subroutine compute_ijklmn(i,j,m,k,l,n,array_mos,array_ints,integral)
 implicit none
 double precision, intent(out) :: integral
 integer, intent(in) :: i,j,m,k,l,n
 double precision, intent(in) :: array_ints(n_points_final_grid,3,mo_num, mo_num)
 double precision, intent(in) :: array_mos(n_points_final_grid,mo_num)
 BEGIN_DOC
! <ijm|L|kln>
 END_DOC
 integer :: ipoint,mm
 integral = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   integral += array_mos(ipoint,i) * array_mos(ipoint,k) * array_ints(ipoint,mm,m,n) * array_ints(ipoint,mm,j,l) 
   integral += array_mos(ipoint,j) * array_mos(ipoint,l) * array_ints(ipoint,mm,m,n) * array_ints(ipoint,mm,i,k) 
   integral += array_mos(ipoint,m) * array_mos(ipoint,n) * array_ints(ipoint,mm,j,l) * array_ints(ipoint,mm,i,k) 
  enddo
 enddo
end
