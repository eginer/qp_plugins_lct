subroutine local_energy_htilde(r1,r2,mu_in,nstates,psi_occ,coefs,ndet,e_loc,kin_e,pot_ee,pot_en,non_hermit_e,psi)
 implicit none
 use bitmasks
 include 'utils/constants.include.F'
 double precision, intent(out):: kin_e(nstates),pot_ee,pot_en,non_hermit_e(nstates),e_loc(nstates),psi(nstates)
 double precision, intent(in) :: r1(3), r2(3), mu_in
 integer, intent(in) :: ndet,nstates,psi_occ(2,ndet)
 double precision, intent(in) :: coefs(ndet,nstates)


 double precision :: r12,r,w_ee_eff
 integer :: i,istate
 double precision :: grad_psi(3,nstates), lapl_psi(3,nstates), d_d_r12(nstates)
 call get_two_e_general_grad_lapl_psi_at_r1r2(r1,r2,mu_in,coefs,nstates,psi_occ,ndet,psi,grad_psi,lapl_psi,d_d_r12)
 do istate = 1, nstates
  kin_e(istate) = -0.5d0 * ( lapl_psi(1,istate) + lapl_psi(2,istate) + lapl_psi(3,istate) ) 
  non_hermit_e(istate) = - d_d_r12(istate)
 enddo
 if(psi_occ(2,1).gt.0)then
  r12 = (r1(1) - r2(1))**2.D0 + (r1(2) - r2(2))**2.D0 + (r1(3) - r2(3))**2.D0
  r12 = dsqrt(r12)
  pot_ee = w_ee_eff(r12,mu_in)
  r = 0.d0
  pot_en = 0.d0
  do i = 1, nucl_num
   r = dsqrt((r1(1) - nucl_coord_transp(1,i))**2.D0 + (r1(2) - nucl_coord_transp(2,i))**2.D0 + (r1(3) - nucl_coord_transp(3,i))**2.D0)
   pot_en += -nucl_charge(i) / r
   r = dsqrt((r2(1) - nucl_coord_transp(1,i))**2.D0 + (r2(2) - nucl_coord_transp(2,i))**2.D0 + (r2(3) - nucl_coord_transp(3,i))**2.D0)
   pot_en += -nucl_charge(i) / r
  enddo
 else
  pot_ee = 0.d0
  pot_en = 0.d0
  do i = 1, nucl_num
   r = dsqrt((r1(1) - nucl_coord_transp(1,i))**2.D0 + (r1(2) - nucl_coord_transp(2,i))**2.D0 + (r1(3) - nucl_coord_transp(3,i))**2.D0)
   pot_en += -nucl_charge(i) / r
  enddo
 endif
 
 do istate = 1, nstates
  e_loc(istate) = pot_en + pot_ee + (non_hermit_e(istate) + kin_e(istate))/psi(istate)
 enddo
end


subroutine get_two_e_general_grad_lapl_psi_at_r1r2(r1,r2,mu_in,coefs,nstates,psi_occ,ndet,psi,grad_psi,lapl_psi,d_d_r12)
 implicit none
 integer, intent(in) :: ndet,nstates
 double precision, intent(in) :: r1(3),r2(3), coefs(ndet,nstates),mu_in
 integer, intent(in) :: psi_occ(2,ndet)
 double precision, intent(out):: psi(nstates),grad_psi(3,nstates),lapl_psi(3,nstates),d_d_r12(nstates)
 BEGIN_DOC
! provides psi(r1,r2), sum_{i=1,2} grad_{x,y,z} psi(r1,r2), sum_{i=1,2} lapl_{x,y,z} psi(r1,r2)
!
! and the non hermi term (1 - erf(mu_in r12)) d/d_r12 psi(r1,r2)
 END_DOC
 
 double precision :: mos_array_r1(mo_num),mos_grad_array_r1(3,mo_num),mos_lapl_array_r1(3,mo_num)
 double precision :: mos_array_r2(mo_num),mos_grad_array_r2(3,mo_num),mos_lapl_array_r2(3,mo_num)
 call give_all_mos_and_grad_and_lapl_at_r(r1,mos_array_r1,mos_grad_array_r1,mos_lapl_array_r1)
 call give_all_mos_and_grad_and_lapl_at_r(r2,mos_array_r2,mos_grad_array_r2,mos_lapl_array_r2)

 integer :: i,istate,i_up,i_down,m,j
 psi = 0.d0
 grad_psi = 0.d0
 lapl_psi = 0.d0
 do istate = 1, nstates
  do i = 1, ndet
   i_up   = psi_occ(1,i)
   i_down = psi_occ(2,i)
   if(i_down.gt.0)then
    psi(istate) += mos_array_r1(i_up) * mos_array_r2(i_down) * coefs(i,istate)
    do m = 1, 3
     grad_psi(m,istate) += coefs(i,istate) * (mos_grad_array_r1(m,i_up)   * mos_array_r2(i_down) & 
                                            + mos_grad_array_r2(m,i_down) * mos_array_r1(i_up)   )
     lapl_psi(m,istate) += coefs(i,istate) * (mos_lapl_array_r1(m,i_up)   * mos_array_r2(i_down) &  
                                            + mos_lapl_array_r2(m,i_down) * mos_array_r1(i_up)   )    
    enddo
   else
    psi(istate) += mos_array_r1(i_up) * coefs(i,istate)
    do m = 1, 3
     grad_psi(m,istate) += coefs(i,istate) * mos_grad_array_r1(m,i_up)
     lapl_psi(m,istate) += coefs(i,istate) * mos_lapl_array_r1(m,i_up)
    enddo
   endif
  enddo
 enddo

 double precision :: r12(3), dist_r12, dist_vec(3),poly(3)
 double precision :: erf_mu_r12,derf_mu_x,poly_tot(3)
 dist_r12 = 0.d0
 do m = 1, 3
  r12(m) = r1(m) - r2(m) 
  dist_r12 += r12(m)*r12(m)
 enddo
 dist_r12 = dsqrt(dist_r12)
 dist_vec(1) = dsqrt(r12(2)*r12(2) + r12(3)*r12(3))
 dist_vec(2) = dsqrt(r12(1)*r12(1) + r12(3)*r12(3))
 dist_vec(3) = dsqrt(r12(1)*r12(1) + r12(2)*r12(2))
 erf_mu_r12 = derf_mu_x(mu_in,dist_r12)
 call inv_r_times_poly(r12, dist_r12, dist_vec, poly)
 ! poly_tot(1) = (1 - erf(mu * r12))/(2 * r12) (x1 - x2)
 do m = 1, 3
  poly_tot(m) = 0.5d0 * (poly(m) - erf_mu_r12 * r12(m) )
 enddo

 d_d_r12 = 0.d0
 do istate = 1, nstates
  do i = 1, ndet
   i_up   = psi_occ(1,i)
   i_down = psi_occ(2,i)
   if(i_down.gt.0)then
    do m = 1, 3
     d_d_r12(istate) += poly_tot(m) * (mos_grad_array_r1(m,i_up)   * mos_array_r2(i_down)  &
                                     - mos_grad_array_r2(m,i_down) * mos_array_r1(i_up )) * coefs(i,istate)
    enddo
   else
    d_d_r12 = 0.d0
   endif
  enddo
 enddo

end

subroutine  give_occ_two_e_psi(dets,ndet,psi_occ)
 implicit none
 use bitmasks
 include 'utils/constants.include.F'
 integer, intent(in) :: ndet
 integer(bit_kind), intent(in)  :: dets(N_int,2,ndet)
 integer, intent(out) :: psi_occ(2, ndet)

 integer :: i,j,n_occ_ab(2)
 integer :: occ(N_int*bit_kind_size,2)
 if(elec_num.gt.2)then
  print*,'this provider works only for two electron WF'
  stop
 endif
 do i = 1, ndet
  call bitstring_to_list_ab(dets(1,1,i), occ, n_occ_ab, N_int)
  psi_occ(1,i) = occ(1,1)
  if(n_occ_ab(2)==0)then
   psi_occ(2,i) = -1
  else
   psi_occ(2,i) = occ(1,2)
  endif
 enddo
end

double precision function w_ee_eff(x,mu_in)
 implicit none
 double precision, intent(in) :: x,mu_in
 double precision :: derf_mu_in_x
 include 'utils/constants.include.F'
 w_ee_eff = derf_mu_in_x(x,mu_in)  + mu_in * inv_sq_pi * dexp(-(mu_in*x)**2.d0) - 0.25d0 * (1.d0 - derf(mu_in*x))**2.d0                  
end

double precision function derf_mu_in_x(x,mu_in)  
 implicit none
 double precision, intent(in) :: mu_in,x
 include 'utils/constants.include.F'
 if(dabs(x).gt.1.d-6)then
  derf_mu_in_x = derf(mu_in * x)/x
 else
  derf_mu_in_x =  inv_sq_pi * 2.d0 * mu_in * (1.d0 - mu_in*mu_in*x*x/3.d0)
 endif                                                                                                                   
end

