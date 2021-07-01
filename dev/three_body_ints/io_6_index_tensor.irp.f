
subroutine write_array_6_index_tensor(n_orb,array_tmp,name_file)
 implicit none
 integer, intent(in) :: n_orb
 character*(128),  intent(in) :: name_file 
 double precision, intent(in) :: array_tmp(n_orb,n_orb,n_orb,n_orb,n_orb,n_orb)

 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 PROVIDE ezfio_filename                                                                                                  
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'W')
 write(i_unit_output)array_tmp
 close(unit=i_unit_output)
end

subroutine read_array_6_index_tensor(n_orb,array_tmp,name_file)
 implicit none
 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 integer, intent(in) :: n_orb
 character*(128),  intent(in)  :: name_file 
 double precision, intent(out) :: array_tmp(n_orb,n_orb,n_orb,n_orb,n_orb,n_orb)
 PROVIDE ezfio_filename                                                                                                  
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'R')
 read(i_unit_output)array_tmp
 close(unit=i_unit_output)
end

