program test_e_two_rdm
 implicit none
 read_wf = .True.
 touch read_wf
 call routine_active_only_test(act_2_rdm_ab_mo)
!call routine_full_mos_test


end
