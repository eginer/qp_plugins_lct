program read_fci_dump
  implicit none
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  output=trim(ezfio_filename)//'.FCIDUMP_erf'
  i_unit_output = getUnitAndOpen(output,'r')

  integer :: i,j,k,l
  integer :: i1,j1,k1,l1
  integer :: i2,j2,k2,l2
  integer*8 :: m
  double precision, allocatable :: big_array(:,:,:,:),v_array(:,:)
  double precision :: core_e,integral
  allocate(big_array(n_act_orb,n_act_orb,n_act_orb,n_act_orb),v_array(n_act_orb,n_act_orb))
  do i = 1, n_lines
    read(i_unit_output,*) integral, i1,k1,j1,l1
    if(j1.ne.0)then
     big_array(i1,k1,j1,l1) = integral
     big_array(k1,i1,j1,l1) = integral
     big_array(i1,k1,l1,j1) = integral
     big_array(k1,i1,l1,j1) = integral
 
     big_array(j1,l1,i1,k1) = integral
     big_array(j1,l1,k1,i1) = integral
     big_array(l1,j1,i1,k1) = integral
     big_array(l1,j1,k1,i1) = integral
    else if(i1.ne.0)then
     v_array(i1,k1) = integral
     v_array(k1,i1) = integral
    else if(i1==0.and.k1==0)then
     core_e = integral
    endif
  enddo
  double precision :: energy
  energy = core_e
  do i = 1, 3
   energy += 2.d0 * v_array(i,i)
   do j = 1, 3
    energy += 2.d0 * big_array(i,i,j,j) - big_array(i,j,i,j)
   enddo
  enddo
  print*,'energy = ',energy
end
