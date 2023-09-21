
program fcidump_tc_h_3_body_5_idx
 implicit none
 my_grid_becke = .True. 
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid 
 call fcidump_3_tc
  call fcidump_2_tc_cst_mu
end

subroutine fcidump_3_tc
 implicit none
 use bitmasks
 integer :: i,j,k,l,m,n
 double precision :: integral 
 character*(128) :: output_physicist
 integer :: i_unit_output_physicist,inint,i_accu
 integer :: getUnitAndOpen
 integer(bit_kind), allocatable :: key(:)
 integer :: ii,jj,kk,ll,mm,nn
 character*(2048)                :: output(1)
 output_physicist =trim(ezfio_filename)//'/FCIDUMP_3_body_tc_5_idx'
 i_unit_output_physicist = getUnitAndOpen(output_physicist,'w')
 allocate(key(N_int))
!if(read_six_index_tensor)then
 do nn = 1, n_act_orb
  n = list_act(nn)
  do ll = 1, n_act_orb
   l = list_act(ll)
   do kk = 1, n_act_orb
    k = list_act(kk)
     do mm = 1, n_act_orb
      m = list_act(mm)
      do jj = 1, n_act_orb
       j = list_act(jj)
       do ii = 1, n_act_orb
        i = list_act(ii)
        key = 0_bit_kind 
        call set_bit_to_integer(i,key,N_int)
        call set_bit_to_integer(j,key,N_int)
        call set_bit_to_integer(m,key,N_int)
        call set_bit_to_integer(l,key,N_int)
        call set_bit_to_integer(n,key,N_int)
        call set_bit_to_integer(k,key,N_int)
        i_accu = 0
        do inint = 1, N_int
         i_accu += popcnt(key(inint))
        enddo
        if(i_accu .gt.5)cycle
        !                          1 2 3 1 2 3
        !                         <i j m|k l n>
        call give_integrals_3_body(i,j,m,k,l,n,integral)
                   
        integral = - integral * 1.d0/3.d0 !!!! For NECI convention 
        if(dabs(integral).lt.1.d-12)cycle
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, ii, jj, mm, kk, ll, nn 


       enddo
      enddo
     enddo
    enddo
   enddo
  enddo

end

