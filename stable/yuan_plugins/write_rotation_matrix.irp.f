program write_rot_mat
 implicit none
 read_wf = .True.
 touch read_wf
! call routine_active_only_test(act_2_rdm_ab_mo)
 call write_rotation_matrix_total_one_e_rdm_uniq
end
