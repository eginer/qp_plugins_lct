
double precision function ao_bielec_integral_ijkl_r3(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !   \int dr^3  i(r) j(r) k(r) l(r)
  END_DOC
  integer,intent(in)             :: i,j,k,l
  integer                        :: p,q,r,s,px,py,pz,qx,qy,qz
  double precision               :: I_center(3),J_center(3),K_center(3),L_center(3)
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  double precision               :: integral,intx,inty,intz,overlap_gaussian_x
  include 'Utils/constants.include.F'
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision, allocatable  :: schwartz_kl(:,:)
  double precision               :: schwartz_ij
  double precision               :: coef2
  double precision               :: coef1
  double precision               :: coef3
  double precision               :: coef4
  double precision               :: p_inv,q_inv

  dim1 = n_pt_max_integrals
  
  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)
  ao_bielec_integral_ijkl_r3 = 0.d0
  double precision               :: thr
  
  allocate(schwartz_kl(0:ao_prim_num(l),0:ao_prim_num(k)))


    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
      I_center(p) = nucl_coord(num_i,p)
      J_center(p) = nucl_coord(num_j,p)
      K_center(p) = nucl_coord(num_k,p)
      L_center(p) = nucl_coord(num_l,p)
    enddo
    

    do p = 1, ao_prim_num(i)
       coef1 = ao_coef_normalized_ordered_transp(p,i)
      do q = 1, ao_prim_num(j)
         coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
        call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
            ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),                 &
            I_power,J_power,I_center,J_center,dim1)
        do r = 1, ao_prim_num(k)
           coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
          do s = 1, ao_prim_num(l)
             coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
                ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),             &
                K_power,L_power,K_center,L_center,dim1)
            
            intx = 0.d0
            inty = 0.d0
            intz = 0.d0
            do px = 0,iorder_p(1) 
             do qx = 0,iorder_q(1)
              intx +=  P_new(px,1) * Q_new(qx,1) * overlap_gaussian_x(P_center(1),Q_center(1),pp,qq,px,qx,max_dim) 
             enddo
            enddo
            do py = 0,iorder_p(2)
             do qy = 0,iorder_q(2)
              inty +=  P_new(py,2) * Q_new(qy,2) * overlap_gaussian_x(P_center(2),Q_center(2),pp,qq,py,qy,max_dim)
             enddo
            enddo
            do pz = 0,iorder_p(3)
             do qz = 0,iorder_q(3)
              intz +=  P_new(pz,3) * Q_new(qz,3) * overlap_gaussian_x(P_center(3),Q_center(3),pp,qq,pz,qz,max_dim)
             enddo
            enddo
            ao_bielec_integral_ijkl_r3 = ao_bielec_integral_ijkl_r3 + fact_p * fact_q * intx * inty * intz * coef4
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p
    
  deallocate (schwartz_kl)
  
end

