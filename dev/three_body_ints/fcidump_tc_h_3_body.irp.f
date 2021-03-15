program fcidump_tc_h_3_body
 implicit none
! call test_mo_erf
 integer :: i,j,k,l,m,n
 double precision :: integral 
 character*(128) :: output_chemist, output_physicist
 integer :: i_unit_output_chemist, i_unit_output_physicist
 integer :: getUnitAndOpen
! output_chemist =trim(ezfio_filename)//'.FCIDUMP_3_body_chemist'
! i_unit_output_chemist = getUnitAndOpen(output_chemist,'w')
 output_physicist =trim(ezfio_filename)//'.FCIDUMP_3_body_physicist'
 i_unit_output_physicist = getUnitAndOpen(output_physicist,'w')
 if(read_six_index_tensor)then
  provide three_body_ints
  do n = 1, mo_num
   do l = 1, mo_num
    do k = 1, mo_num
     do m = 1, mo_num
      do j = 1, mo_num
       do i = 1, mo_num
        !                          1 2 3 1 2 3
        !                         <i j m|k l n>
        !                         (ik|jl|mn)
        integral = three_body_ints(i,j,m,k,l,n)
        integral = integral * 1.d0/3.d0 !!!! For NECI convention 
        if(dabs(integral).lt.1.d-12)cycle
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, i, j, m, k, l, n
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 else

  do n = 1, mo_num
   do l = 1, mo_num
    do k = 1, mo_num
     do m = 1, mo_num
      do j = 1, mo_num
       do i = 1, mo_num
        !                          1 2 3 1 2 3
        !                         <i j m|k l n>
        !                         (ik|jl|mn)
        call give_integrals_3_body(i,j,m,k,l,n,integral)
        integral = -1.D0 * integral * 1.d0/3.d0 !!!! For NECI convention 
        if(dabs(integral).lt.1.d-12)cycle
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, i, j, m, k, l, n
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 endif

end
