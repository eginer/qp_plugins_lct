program density_mulliken
 implicit none
  read_wf = .True.
  touch read_wf
  call routine 

end 

subroutine routine
 implicit none
 double precision :: accu_tot, accu_charge
 integer :: i,istate
 do istate = 1, N_states
  print*,'Mulliken density analysis for state ',istate
  accu_charge = 0.d0
  accu_tot = 0.d0
  do i = 1, nucl_num
   print*,i,nucl_charge(i),nucl_charge(i)-mulliken_density_densities(i,istate),nucl_coord(i,3)
   accu_charge += nucl_charge(i)-mulliken_density_densities(i,istate)
   accu_tot+= mulliken_density_densities(i,istate)
  enddo
  print*,'net charge = ',accu_charge
  print*,'total elec = ',accu_tot
 enddo
end
