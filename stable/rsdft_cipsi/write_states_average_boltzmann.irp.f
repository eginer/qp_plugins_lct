program toucher_les_poids
 implicit none
 BEGIN_DOC
 ! Documentation tres precise 
 END_DOC
 read_wf = .true.
 touch read_wf
 state_average_weight = state_average_weight_boltzmann
 touch state_average_weight 
 call pouetou
 print*,  '**************************************        **'
 write(*, '(A22,X,3(F16.10,X))') 'Weight boltz  1 = ',state_average_weight_boltzmann(1)
 write(*, '(A22,X,3(F16.10,X))') 'Weight boltz  2 = ',state_average_weight_boltzmann(2)
 write(*, '(A22,X,3(F16.10,X))') 'Weight boltz  3 = ',state_average_weight_boltzmann(3)
 print*,  '**************************************        **'
end
  
  subroutine pouetou
  implicit none
 call ezfio_set_determinants_state_average_weight(state_average_weight)
end
