subroutine s_plus_det_orb_j(det_in,orb_j,det_out,phase)
 implicit none
 BEGIN_DOC
! a^dagger_orb_j_up a_orb_j_down | det_in > = phase * | det_out >
 END_DOC
  use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), intent(in)  :: det_in(N_int,2)
 integer          , intent(in)  :: orb_j
 integer(bit_kind), intent(out) :: det_out(N_int,2)
 double precision, intent(out)  :: phase
 
 integer(bit_kind) :: open_shell(N_int),open_a_b(N_int,2)
 integer(bit_kind) :: det_orb_j(N_int)
 integer :: i,j,n_occ_ab(2),ispin,occ(N_int*bit_kind_size,2)

 det_out = det_in 

! open_shell = 0_bit_kind
 ! string with open shell orbitals 
 do i = 1, N_int
  open_shell(i) = xor(det_in(i,1),det_in(i,2))
 enddo
 ! string with open shell orbitals of alpha and beta electrons
 do i = 1, N_int
  open_a_b(i,1) = iand(open_shell(i),det_in(i,1))
  open_a_b(i,2) = iand(open_shell(i),det_in(i,2))
 enddo
 ! check if the excitation is possible 
 ! ==> is the orbital orb_j occupied in open_shell_beta 
 det_orb_j = 0_bit_kind
 call set_bit_to_integer(orb_j,det_orb_j,N_int)
 j = 0
 do i = 1, N_int
  j+= popcnt(iand(det_orb_j(i),open_a_b(i,2)))
 enddo
 if(j.ne.1)then
  return
  phase = 0.d0
 endif

 ! you do the spin flip : det_out is initialized to det_in
 ! set to 1 the bit corresponding to orb_j into the beta string
 call set_bit_to_integer(orb_j,det_out(1,1),N_int)
 ! clear the bit corresponding to orb_j into the alpha string
 call clear_bit_to_integer(orb_j,det_out(1,2),N_int)
 
 call bitstring_to_list_ab(open_a_b, occ, n_occ_ab, N_int)
 phase = 1.d0
 ! phase = (-1)**(number open shell orbitals before the orbital orb_j
 do ispin = 1, 2
  do i = 1, n_occ_ab(ispin)
   j = occ(i,ispin)
   if(j.lt.orb_j)then
    phase *= -1.d0
   endif
  enddo
 enddo
end


subroutine s_plus_det(det_in,det_out,phase,ndet_out)
 implicit none
 BEGIN_DOC
! a^dagger_orb_j_up a_orb_j_down | det_in > = phase * | det_out >
 END_DOC
  use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), intent(in)  :: det_in(N_int,2)
 integer(bit_kind), intent(out) :: det_out(N_int,2,elec_beta_num)
 integer, intent(out)           :: ndet_out
 double precision , intent(out) :: phase(elec_beta_num)
 
 integer(bit_kind) :: open_shell(N_int),open_a_b(N_int,2),det_tmp(N_int,2)
 integer(bit_kind) :: det_orb_j(N_int)
 double precision  :: phase_tmp
 integer :: i,j,n_occ_ab(2),ispin,occ(N_int*bit_kind_size,2),k
! call debug_det(det_in,N_int)
 ! string with open shell orbitals 
 do i = 1, N_int
  open_shell(i) = xor(det_in(i,1),det_in(i,2))
 enddo
 ! string with open shell orbitals of alpha and beta electrons
 do i = 1, N_int
  open_a_b(i,1) = iand(open_shell(i),det_in(i,1))
  open_a_b(i,2) = iand(open_shell(i),det_in(i,2))
 enddo
 call bitstring_to_list_ab(open_a_b, occ, n_occ_ab, N_int)
! call debug_det(open_a_b,N_int)
 ndet_out = 0
 do i = 1, n_occ_ab(2)
  j = occ(i,2)
  call s_plus_det_orb_j(det_in,j,det_tmp,phase_tmp)
  if(phase_tmp .ne. 0.d0)then
   ndet_out += 1
   det_out(:,:,ndet_out) = det_tmp(:,:) 
   phase(ndet_out) = phase_tmp
  endif
 enddo

end


double precision function factor_s_p(S,ms)
 implicit none
 double precision, intent(in) :: S,ms
 factor_s_p = dsqrt(S*(S+1.d0) - ms*(ms+1.d0))
end

double precision function factor_s_m(S,ms)
 implicit none
 double precision, intent(in) :: S,ms
 factor_s_m = dsqrt(S*(S+1.d0) - ms*(ms-1.d0))
end
