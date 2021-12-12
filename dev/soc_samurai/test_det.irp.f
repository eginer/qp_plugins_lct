program test_det
 implicit none
 read_wf = .True.
 touch read_wf 
 call routine
end

subroutine routine
 implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
 integer                        :: exc(0:2,2,2)
 integer                        :: h2,p2,s2,degree,i,j
 integer(bit_kind), allocatable :: det_out(:,:,:)
 double precision, allocatable  :: phase(:)
 double precision :: accu
 allocate(phase(elec_beta_num),det_out(N_int,2,elec_beta_num))
 accu = 0.d0
 do i = 1, N_det
  call s_plus_det(psi_det(1,1,i),det_out,phase)
  print*,''
  print*,'Applying S+ on '
  call debug_det(psi_det(1,1,i),N_int)
  do j = 1, elec_beta_num
   print*,'phase = ',phase
   call debug_det(det_out(1,1,j),N_int)
   accu += psi_coef(i,1) * phase(j)
  enddo
 enddo
 print*,'accu = ',accu
 print*,'sqrt(3)',dsqrt(3.d0)
end
