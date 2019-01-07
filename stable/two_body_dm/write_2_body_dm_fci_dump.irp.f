program write_2_body_dm_fci_dump
 implicit none
 BEGIN_DOC
! This program prints the alpha/beta two-body density in the same format as an FCIDUMP
 END_DOC
 read_wf = .True.
 touch read_wf
 call routine_print

end

subroutine routine_print
  implicit none
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  output=trim(ezfio_filename)//'.2RDM'
  i_unit_output = getUnitAndOpen(output,'w')

  integer :: i,j,k,l
  integer :: i1,j1,k1,l1
  integer :: i2,j2,k2,l2
  integer*8 :: m
  character*(2), allocatable :: A(:)

  write(i_unit_output,*) '&FCI NORB=', n_act_orb, ', NELEC=', elec_num-n_core_orb*2, &
   ', MS2=', (elec_alpha_num-elec_beta_num), ','
  allocate (A(n_act_orb))
  A = '1,'
  write(i_unit_output,*) 'ORBSYM=', (A(i), i=1,n_act_orb) 
  write(i_unit_output,*) 'ISYM=0,'
  write(i_unit_output,*) '/'
  deallocate(A)
  

  double precision :: integral

  do l=1,mo_num
   do k=1,mo_num
    do j=1,mo_num
     do i=1,mo_num
      !                                                     1 2 1 2
      integral = two_bod_alpha_beta_mo_physicist(i,j,k,l,1)
      if(dabs(integral).lt.1.d-16)then
        integral = 0.d0
      endif
      write(i_unit_output,*) integral, i,k,j,l
     enddo
    enddo
   enddo
  enddo
end
