program print_mat
 implicit none
 read_wf = .true.
 touch read_wf
 ! total one-e integrals 
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  
 ! Vne integrals on the MO basis 
 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 ! kinetic integrals on the MO basis 
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 ! Vne integrals on the AO basis 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 ! kinetic integrals on the AO basis 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 

 ! regular 1/r12 integrals  on the MO basis
 io_mo_two_e_integrals = "None"
 touch io_mo_two_e_integrals
 ! regular 1/r12 integrals  on the AO basis
 io_ao_two_e_integrals = "None"
 touch io_ao_two_e_integrals
 ! integral of the effective potential 
 io_mo_int_mu_of_r = "None" 
 touch io_mo_int_mu_of_r
 call routine

end

subroutine routine
 implicit none
  use bitmasks
 integer :: i,j
 double precision              :: h_eff_one_e,h_eff_two_e

 double precision, allocatable :: h_dress_one(:,:),h_dress_two(:,:),h_tot(:,:)
 allocate(h_dress_one(N_det,N_det),h_dress_two(N_det,N_det),h_tot(N_det,N_det))
 provide tot_eff_two_e
 if(.True.)then
  do i = 1, N_det
   do j = 1, N_det
    call i_eff_H_j(psi_det(1,1,i),psi_det(1,1,j),h_eff_one_e,h_eff_two_e)
    h_dress_one(i,j) = h_eff_one_e 
    h_dress_two(i,j) = h_eff_two_e
   enddo
  enddo
  print*,'Regular H matrix '
  do i = 1, N_det
   write(*,'(100(F14.7,X))')H_matrix_all_dets(i,:)
  enddo
  print*,'Dresising matrix one e '
  do i = 1, N_det
   write(*,'(100(F14.7,X))')h_dress_one(i,:)
  enddo
  print*,'Dresising matrix two e '
  do i = 1, N_det
   write(*,'(100(F14.7,X))')h_dress_two(i,:)
  enddo
  h_tot = H_matrix_all_dets + h_dress_one + h_dress_two
  print*,'Total matrix two e '
  do i = 1, N_det
   write(*,'(100(F14.7,X))')h_tot(i,:)
  enddo
  double precision :: eigvalues(N_det),eigvectors(N_det,N_det)
  double precision :: de_0,de_0_one,de_0_two,de_tot
  de_0 = H_matrix_all_dets(3,3) - H_matrix_all_dets(1,1)
  de_0_one = h_dress_one(3,3) - h_dress_one(1,1)
  de_0_two = h_dress_two(3,3) -  h_dress_two(1,1) 
  de_tot   = h_tot(3,3) - h_tot(1,1)
  print*,'Delta E naked = ',de_0
  print*,'Delta E one   = ',de_0_one
  print*,'Delta E two   = ',de_0_two
  print*,'Delta E TOT   = ',de_tot
  call lapack_diagd(eigvalues,eigvectors,h_tot,N_det,N_det)
  print*,'eigvalues = ',eigvalues(1) + nuclear_repulsion
  do i = 1, N_det
   print*,'ci = ',eigvectors(i,1)
  enddo
  print*,'int_d_dn2_e_cmd_basis = ',int_d_dn2_e_cmd_basis
  call write_stuff
 endif


end

subroutine write_stuff
 implicit none
 integer :: i,j,istate,nx,m
 double precision :: r(3),dx,xmax,f_psi,rho2,mu_correction_of_on_top,d_dn2
 double precision :: mu_of_r,rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2
 double precision :: mos_array(mo_num)
 double precision  :: aos_array(ao_num)
 double precision  :: grad_aos_array(ao_num,3)
 istate = 1
 nx = 1000
 xmax = 20.d0
 dx = xmax/dble(nx)
 r(:) = nucl_coord_transp(:,1) 
 r(3) = 0.5d0 * (nucl_coord_transp(3,1) + nucl_coord_transp(3,2) ) 
 r(3) += - xmax * 0.5d0
 write(33,'((A400))')'#       r(3)        rho_a+rho_b          rho2          ec_srmuPBE       decdrho2          mu_of_r            mos_array(1)     mos_array(2)'
 do i = 1, nx
  call give_all_mos_at_r(r,mos_array)
  call give_mu_of_r_cas(r,istate,mu_of_r,f_psi,rho2)
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a , grad_rho_b, aos_array, grad_aos_array)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo
   
   rho2 = 2.d0*rho2
   rho2 = mu_correction_of_on_top(mu_of_r,rho2) ! extrapolation based on mu(r)
   call ecmdsrPBEn2(mu_of_r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)
   d_dn2 = 2.d0 * decdrho2
!   if(mu_of_r.gt.1.d+5)then
!    mu_of_r = 0.d0
!   endif
   write(33,'(100(F16.10,X))')r(3),rho_a+rho_b,rho2,ec_srmuPBE,d_dn2,mu_of_r,mos_array(1),mos_array(2), mos_array(3), mos_array(4)
  r(3) += dx
 enddo
end
