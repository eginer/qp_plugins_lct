BEGIN_PROVIDER[double precision, kernel_lr, (N_det, N_det, N_states)]
 implicit none
 BEGIN_DOC
 ! sum_I(1 -> n_det) sum_J(1 -> n_det) c_I c_J <phi_I | E_ij | phi_J>
 END_DOC 
 integer :: i_state
 integer :: det_I, det_J
 integer :: i, j, k, l
 
 kernel_lr = 0.d0

 do i_state=1, N_states
  do det_I = 1, N_det
   do det_J = 1, N_det

     do i=1, mo_num
      do j=1, mo_num
       do k=1, mo_num
        do l=1, mo_num
          kernel_lr(det_I, det_J, i_state) = kernel_lr(det_I, det_J, i_state) + &
          & delta_gamma_I_i_j(det_I, j, i, i_state) * delta_gamma_I_i_j(det_J, l, k, i_state)* eff_two_e_lr(l,k,j,i,i_state)
        enddo
       enddo
      enddo
     enddo

   enddo
  enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER[double precision, phase_from_hf, (N_det-1)]
&BEGIN_PROVIDER[integer, ex_from_hf, (2, N_det-1)]
 implicit none
  integer                        :: exc(0:2,2,2)
  double precision               :: phase
  integer                        :: h2,p2,s2,degree
  integer                        :: h1,p1,s1
  integer :: det_I, det_II
  
 ex_from_hf = -1
 det_II=0
 do det_I = 1, N_det
  if(det_I.eq.hf_index) cycle
  det_II+=1
  
  call get_excitation_degree(ref_bitmask,psi_det(1,1,det_I),degree,N_int)
  if (degree.eq.1) then
   call get_excitation(ref_bitmask,psi_det(1,1,det_I),exc,degree,phase,N_int)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2) 
   print*,"---------------"
   call debug_det(psi_det(1,1,det_I),N_int)
   print*, h1, p1, s1
   ex_from_hf(1, det_II) = h1 
   ex_from_hf(2, det_II) = p1 
   phase_from_hf(det_II) = phase
  endif
 enddo
 
END_PROVIDER

BEGIN_PROVIDER[double precision, K_test_cis, (N_det-1, N_det-1)]
 implicit none
 integer :: h1,p1,h2,p2
 integer :: det_II,det_JJ
 double precision :: phase1, phase2

 K_test_cis =0.d0

 do det_II = 1, N_det-1  
   if (ex_from_hf(1, det_II).eq.-1) cycle
   h1 = ex_from_hf(1, det_II)
   p1 = ex_from_hf(2, det_II)
   phase1 = phase_from_hf(det_II)
   
  do det_JJ = 1, N_det -1
   if (ex_from_hf(1, det_JJ).eq.-1) cycle
   h2 = ex_from_hf(1, det_JJ)
   p2 = ex_from_hf(2, det_JJ)
   phase2 = phase_from_hf(det_JJ)

   K_test_cis(det_JJ, det_II) = phase1*phase2*eff_two_e_lr(h1,p1,h2,p2,1)
  enddo

 enddo

END_PROVIDER
