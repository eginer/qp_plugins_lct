program cisd_sc2
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  touch read_wf 
  call routine
end

subroutine routine
 implicit none
 print*,'ecorr_tot = ',ecorr_tot
end
