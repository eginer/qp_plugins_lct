program print_mat
 implicit none
 read_wf = .true.
 touch read_wf
 call routine_print
end

subroutine routine_print
 implicit none
 integer :: i,j
  print*,'Regular H matrix '
  do i = 1, N_det
   write(*,'(100(F14.7,X))')H_matrix_all_dets(i,:)
  enddo
  double precision :: eigvalues(N_det),eigvectors(N_det,N_det)
  print*,'Delta E =',H_matrix_all_dets(1,1) - H_matrix_all_dets(3,3)
  call lapack_diagd(eigvalues,eigvectors,H_matrix_all_dets,N_det,N_det)
  print*,'eigvalues = ',eigvalues(1) + nuclear_repulsion
  do i = 1, N_det
   print*,'ci = ',eigvectors(i,1)
  enddo



end
