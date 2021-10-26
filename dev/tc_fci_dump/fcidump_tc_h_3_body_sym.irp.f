program fcidump_tc_h_3_body_sym
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
 integer :: ii,jj,kk,ll,mm,nn
 double precision :: integral 
 character*(128) :: output_physicist
 integer :: i_unit_output_physicist
 integer :: getUnitAndOpen
 output_physicist =trim(ezfio_filename)//'/FCIDUMP_3_body_tc_sym'
 i_unit_output_physicist = getUnitAndOpen(output_physicist,'w')
  do nn = 1, n_act_orb
   n = list_act(nn)
   do ll = 1, n_act_orb
    l = list_act(ll)
    do kk = 1, n_act_orb
     k = list_act(kk)
     do mm = n, n_act_orb
      m = list_act(mm)
      do jj = l, n_act_orb
       j = list_act(jj)
       do ii = k, n_act_orb
        i = list_act(ii)
        !                          1 2 3 1 2 3
        !                         <i j m|k l n>
        !                         (ik|jl|mn)
!        integral = three_body_ints(i,j,m,k,l,n)
        call give_integrals_3_body(i,j,m,k,l,n,integral)
        integral = -integral * 1.d0/3.d0 !!!! For NECI convention 
        if(dabs(integral).lt.1.d-12)cycle
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, ii,jj,mm,kk,ll,nn
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, kk,jj,mm,ii,ll,nn
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, kk,ll,mm,ii,jj,nn
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, kk,jj,nn,ii,ll,mm
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, kk,ll,nn,ii,jj,mm
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, ii,ll,mm,kk,jj,nn
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, kk,ll,mm,ii,jj,nn
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, ii,ll,nn,kk,jj,mm
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, ii,jj,nn,kk,ll,mm
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, kk,jj,nn,ii,ll,mm
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, ii,ll,nn,kk,jj,mm
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
end

