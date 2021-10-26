program fcidump_tc_h_3_body
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
 integer :: i,j,k,l,m,n
 double precision :: integral 
 character*(128) :: output_physicist
 integer :: i_unit_output_physicist
 integer :: getUnitAndOpen
 integer :: ii,jj,kk,ll,mm,nn
 output_physicist =trim(ezfio_filename)//'/FCIDUMP_3_body_tc'
 i_unit_output_physicist = getUnitAndOpen(output_physicist,'w')
!if(read_six_index_tensor)then
  provide three_body_ints
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
        !                          1 2 3 1 2 3
        !                         <i j m|k l n>
        !                         (ik|jl|mn)
        integral = three_body_ints(i,j,m,k,l,n)
        integral = integral * 1.d0/3.d0 !!!! For NECI convention 
        if(dabs(integral).lt.1.d-12)cycle
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, ii, jj, mm, kk, ll, nn
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo

end

