program print_md_sr_pbe
 implicit none
 read_wf = .True.
 touch read_wf 
 call routine

end

subroutine routine
 implicit none
 print*,'energy_x_md_sr_pbe = ',energy_x_md_sr_pbe
 print*,'energy_c_md_sr_pbe = ',energy_c_md_sr_pbe

end
