 BEGIN_PROVIDER [ double precision, dressing_column_h, (N_det,N_states) ]
&BEGIN_PROVIDER [ double precision, dressing_column_s, (N_det,N_states) ]
 implicit none
 BEGIN_DOC
 ! \Delta_{state-specific}. \Psi
 ! Diagonal element is divided by 2 because Delta = D + D^t
 END_DOC

 integer :: i,ii,k,j, l
 double precision :: f, tmp
 double precision, allocatable :: delta(:)

 allocate(delta(N_det))
 call get_delta_no_store_no_provide(psi_det,psi_coef,N_det,delta)

 dressing_column_h(:,:) = 0.d0
 dressing_column_s(:,:) = 0.d0

 l = 1
 do j = 1, n_det
   if (j == l) cycle
   dressing_column_h(j,1)  = delta(j) 
   dressing_column_h(l,1) -= psi_coef(j,1) * delta(j) / psi_coef(l,1)
 enddo
 dressing_column_h(l,1) += delta(l) 
 dressing_column_h(l,1) *= 0.5d0

END_PROVIDER

subroutine set_dressing_column_h_s(u_in)
 implicit none
 double precision, intent(in) :: u_in(N_det)
 integer :: i,ii,k,j, l
 double precision :: f, tmp
 double precision, allocatable :: delta(:)

 allocate(delta(N_det))
 call get_delta_no_store_no_provide(psi_det,u_in,N_det,delta)

 dressing_column_h(:,:) = 0.d0
 dressing_column_s(:,:) = 0.d0

 l = 1
 do j = 1, n_det
   if (j == l) cycle
   dressing_column_h(j,1)  = delta(j) 
   dressing_column_h(l,1) -= u_in(j) * delta(j) / u_in(l)
 enddo
 dressing_column_h(l,1) += delta(l) 
 dressing_column_h(l,1) *= 0.5d0
 soft_touch dressing_column_h 
end


