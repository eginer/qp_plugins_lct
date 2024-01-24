 program print_tc_dm
   implicit none
   BEGIN_DOC
 ! TODO : Put the documentation of the program here
   END_DOC
   print *, 'Hello world'
   my_grid_becke = .True.
   my_n_pt_r_grid = 30
   my_n_pt_a_grid = 50
   read_wf = .True.
   touch read_wf
   touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
   call print_dm
 end

subroutine print_dm
 implicit none
 integer :: i
 do i = 1, elec_alpha_num
   write(33,*) dabs(tc_transition_matrix(i,i,1,1)-2.d0)
 enddo
 do i = 1+ elec_alpha_num, mo_num
   write(33,*) dabs(tc_transition_matrix(i,i,1,1))
 enddo
end
