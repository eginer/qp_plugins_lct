BEGIN_PROVIDER [integer, n_soc_states]
 implicit none
 integer :: i,n_sm, n_sp, how_many_sm,how_many_sp
 double precision :: S,ms
 ms = dble(elec_alpha_num - elec_beta_num)
 n_soc_states = 0
 do i = 1, n_states
  S = s_values(i)
  n_soc_states += 1
  n_sm = how_many_sm(S,ms)
  n_sp = how_many_sp(S,ms)
  print*,'S,n_sm,n_sp',S,n_sm,n_sp
  n_soc_states += n_sm + n_sp
 enddo

END_PROVIDER 

