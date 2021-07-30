program write_2_rdm
 implicit none
 read_wf = .True.
 touch read_wf
 call write_two_rdm(act_2_rdm_spin_trace_mo)
end
