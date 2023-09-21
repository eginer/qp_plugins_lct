program save_one_e_dm
  implicit none
 BEGIN_DOC
! Program that computes the one body density on the |MO| and |AO| basis
! for $\alpha$ and $\beta$ electrons from the wave function
! stored in the |EZFIO| directory, and then saves it into the
! :ref:`module_aux_quantities`.
!
! Then, the global variable :option:`aux_quantities data_one_e_dm_alpha_mo`
! and :option:`aux_quantities data_one_e_dm_beta_mo` (and the corresponding for |AO|)
! will automatically ! read this density in the next calculation. 
! This can be used to perform damping on the density in |RSDFT| calculations (see
! :ref:`module_density_for_dft`).
 END_DOC
  read_wf = .True.
  touch read_wf
  call routine_save_one_e_dm_fc

end

