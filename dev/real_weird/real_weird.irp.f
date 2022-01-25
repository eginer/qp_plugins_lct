program real_weird
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print*,' e0_effective_H  =',e0_effective_H
  integer :: i
  do i = 1, mo_num
   write(*,'(100(F16.10,X))'),effective_H(i,:)
  enddo
 
 character*(64)    :: label
 logical :: output
 output = .True.
 label = "Canonical"
 call mo_as_eigvectors_of_mo_matrix(effective_H,mo_num,mo_num,label,1.d0,output)
! call save_mos                                                                                                          
end

