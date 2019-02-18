use bitmasks

 BEGIN_PROVIDER [ integer(bit_kind), psi_ref_fobo, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_ref_fobo_coef,  (psi_det_size,n_states) ]
&BEGIN_PROVIDER [ integer, idx_ref_fobo_cas_lmct_mlct, (psi_det_size,2) ]
&BEGIN_PROVIDER [ integer, idx_ref_fobo, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, idx_ref_cas_fobo, (psi_det_size) ]
&BEGIN_PROVIDER [ integer(bit_kind), psi_non_ref_fobo, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_non_ref_fobo_coef,  (psi_det_size,n_states) ]
&BEGIN_PROVIDER [ integer, idx_non_ref_fobo, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, N_det_non_ref_fobo ]
&BEGIN_PROVIDER [ integer, idx_lmct, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, idx_mlct, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, idx_lmct_orb, (n_list_lmct) ]
&BEGIN_PROVIDER [ integer, idx_mlct_orb, (n_list_mlct) ]
&BEGIN_PROVIDER [ integer, idx_lmct_orb_reverse, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, idx_mlct_orb_reverse, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, N_det_ref_fobo ]
&BEGIN_PROVIDER [ integer, N_det_ref_cas_fobo]
&BEGIN_PROVIDER [ integer, N_det_lmct ]
&BEGIN_PROVIDER [ integer, N_det_mlct ]
  implicit none
  BEGIN_DOC
  ! CAS wave function, defined from the application of the CAS bitmask on the 
  ! determinants. idx_cas gives the indice of the CAS determinant in psi_det.
  END_DOC
  integer :: i,j,k
  integer :: i_hole,i_particl,n_h,n_p,number_of_holes,number_of_particles
  logical :: is_the_hole_in_det,is_the_particl_in_det
  integer :: n_cas_tmp
  integer :: hole_particl_count(0:2,0:2)
  integer :: i_count_1h
  integer :: i_count_1p
  double precision :: coef_average

  N_det_ref_fobo = 0
  N_det_ref_cas_fobo = 0
  N_det_lmct = 0
  N_det_mlct = 0
  N_det_non_ref_fobo = 0
  n_cas_tmp = 0
  hole_particl_count = 0
  i_count_1h = 0
  i_count_1p = 0
  do j = 1, N_det
   n_h = number_of_holes(psi_det(1,1,j))
   n_p = number_of_particles(psi_det(1,1,j))
!  print*, n_h,n_p
   hole_particl_count(n_h,n_p) += 1
   
   if(n_h==1.and.n_p==0)then
    i_count_1h += 1
    logical :: is_good_hole
    is_good_hole = .False.
    do i = 1, n_list_lmct
     i_hole = list_lmct(i)
     if(is_the_hole_in_det(psi_det(1,1,j),1,i_hole).or. &
        is_the_hole_in_det(psi_det(1,1,j),2,i_hole))then
      coef_average = 0.d0
      do k = 1, N_states
       coef_average += dabs(psi_coef(j,k))
      enddo
      coef_average = coef_average / dble(N_states)
      if(coef_average.gt.thresh_coef_lmct)then
       N_det_ref_fobo +=1 
       N_det_lmct+= 1
       idx_ref_fobo(N_det_ref_fobo) = j
       idx_lmct(N_det_lmct) = j
       idx_lmct_orb(i) = j
       idx_lmct_orb_reverse(N_det_lmct) = i_hole
       idx_ref_fobo_cas_lmct_mlct(N_det_ref_fobo,1) = +1
       idx_ref_fobo_cas_lmct_mlct(N_det_ref_fobo,2) = i_hole
       do k=1,N_int
         psi_ref_fobo(k,1,N_det_ref_fobo) = psi_det(k,1,j)
         psi_ref_fobo(k,2,N_det_ref_fobo) = psi_det(k,2,j)
       enddo
       do k=1,N_states
         psi_ref_fobo_coef(N_det_ref_fobo,k) = psi_coef(j,k)
       enddo
       is_good_hole = .True.
       exit
      endif
     endif
    enddo
    if(.not.is_good_hole)then
     N_det_non_ref_fobo += 1
     idx_non_ref_fobo(N_det_non_ref_fobo) = j
     do k=1,N_int
       psi_non_ref_fobo(k,1,N_det_non_ref_fobo) = psi_det(k,1,j)
       psi_non_ref_fobo(k,2,N_det_non_ref_fobo) = psi_det(k,2,j)
     enddo
     do k=1,N_states
       psi_non_ref_fobo_coef(N_det_non_ref_fobo,k) = psi_coef(j,k)
     enddo
    endif
   else if (n_h==0.and.n_p==1)then
    i_count_1p += 1
    logical :: is_good_particl
    is_good_particl = .False.
    do i = 1, n_list_mlct
     i_particl = list_mlct(i)
     if(is_the_particl_in_det(psi_det(1,1,j),1,i_particl).or. &
        is_the_particl_in_det(psi_det(1,1,j),2,i_particl))then
      coef_average = 0.d0
      do k = 1, N_states
       coef_average += dabs(psi_coef(j,k))
      enddo
      coef_average = coef_average / dble(N_states)
      if(coef_average.gt.thresh_coef_mlct)then
       N_det_ref_fobo +=1 
       N_det_mlct+= 1
       idx_ref_fobo(N_det_ref_fobo) = j
       idx_mlct(N_det_mlct) = j
       idx_mlct_orb(i) = j
       idx_mlct_orb_reverse(N_det_mlct) = i_particl
       idx_ref_fobo_cas_lmct_mlct(N_det_ref_fobo,1) = -1
       idx_ref_fobo_cas_lmct_mlct(N_det_ref_fobo,2) = i_particl
       do k=1,N_int
         psi_ref_fobo(k,1,N_det_ref_fobo) = psi_det(k,1,j)
         psi_ref_fobo(k,2,N_det_ref_fobo) = psi_det(k,2,j)
       enddo
       do k=1,N_states
         psi_ref_fobo_coef(N_det_ref_fobo,k) = psi_coef(j,k)
       enddo
       is_good_particl = .True.
       exit
      endif
     endif
    enddo
    if(.not.is_good_particl)then
     N_det_non_ref_fobo += 1
     idx_non_ref_fobo(N_det_non_ref_fobo) = j
     do k=1,N_int
       psi_non_ref_fobo(k,1,N_det_non_ref_fobo) = psi_det(k,1,j)
       psi_non_ref_fobo(k,2,N_det_non_ref_fobo) = psi_det(k,2,j)
     enddo
     do k=1,N_states
       psi_non_ref_fobo_coef(N_det_non_ref_fobo,k) = psi_coef(j,k)
     enddo
    endif
   else if (n_h==0.and.n_p==0)then
    coef_average = 0.d0
    do k = 1, N_states
     coef_average += dabs(psi_coef(j,k))
    enddo
    coef_average = coef_average / dble(N_states)
    if(coef_average.gt.thresh_coef_cas)then
     N_det_ref_fobo +=1 
     N_det_ref_cas_fobo += 1
     idx_ref_fobo(N_det_ref_fobo) = j
     idx_ref_cas_fobo(N_det_ref_fobo) = j
     idx_ref_fobo_cas_lmct_mlct(N_det_ref_fobo,1) = 0 
     do k=1,N_int
       psi_ref_fobo(k,1,N_det_ref_fobo) = psi_det(k,1,j)
       psi_ref_fobo(k,2,N_det_ref_fobo) = psi_det(k,2,j)
     enddo
     do k=1,N_states
       psi_ref_fobo_coef(N_det_ref_fobo,k) = psi_coef(j,k)
     enddo
    else 
     N_det_non_ref_fobo += 1
     idx_non_ref_fobo(N_det_non_ref_fobo) = j
     do k=1,N_int
       psi_non_ref_fobo(k,1,N_det_non_ref_fobo) = psi_det(k,1,j)
       psi_non_ref_fobo(k,2,N_det_non_ref_fobo) = psi_det(k,2,j)
     enddo
     do k=1,N_states
       psi_non_ref_fobo_coef(N_det_non_ref_fobo,k) = psi_coef(j,k)
     enddo
    endif
   else 
    N_det_non_ref_fobo += 1
    idx_non_ref_fobo(N_det_non_ref_fobo) = j
    do k=1,N_int
      psi_non_ref_fobo(k,1,N_det_non_ref_fobo) = psi_det(k,1,j)
      psi_non_ref_fobo(k,2,N_det_non_ref_fobo) = psi_det(k,2,j)
    enddo
    do k=1,N_states
      psi_non_ref_fobo_coef(N_det_non_ref_fobo,k) = psi_coef(j,k)
    enddo
   endif
  enddo
  integer :: n_det_tmp
  n_det_tmp = 0
  do n_h = 0,2 
   do n_p = 0,2
   print*, 'n_h,n_p',n_h,n_p
   print*, hole_particl_count(n_h,n_p) 
   n_det_tmp += hole_particl_count(n_h,n_p)
   enddo
  enddo
  print*, 'N_det_cas           ',N_det_cas
  print*, 'N_det_ref_cas_fobo  ',N_det_ref_cas_fobo
  print*, 'N_det_lmct          ',N_det_lmct
  print*, 'N_det_mlct          ',N_det_mlct
  if(N_det_non_ref_fobo + N_det_ref_fobo.ne.N_det)then
   print*, 'pb !!'
   print*, 'N_det_non_ref_fobo + N_det_ref_fobo.ne.N_det'
   print*,  N_det_non_ref_fobo + N_det_ref_fobo, N_det 
   print*,  N_det_non_ref_fobo , N_det_ref_fobo, N_det 
   print*,  'N_det_cas,N_det_mlct,N_det_lmct'
   print*,  N_det_cas,N_det_mlct,N_det_lmct
   stop
  endif

END_PROVIDER


BEGIN_PROVIDER [ double precision, psi_ref_fobo_coef_inv, (psi_det_size,n_states) ]
 implicit none
 BEGIN_DOC
 ! 1/psi_ref_fobo_coef
 END_DOC
 integer :: i, i_state
 do i_state=1,N_states
  do i=1,N_det_ref_fobo
    psi_ref_fobo_coef_inv(i,i_state) = 1.d0/psi_ref_fobo_coef(i,i_state)
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), psi_ref_fobo_restart, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_ref_fobo_coef_restart,  (psi_det_size,n_states) ]
  implicit none
  BEGIN_DOC
  ! Projection of the CAS wave function on the restart wave function. 
  END_DOC
  integer :: i,j,k
  integer, save                  :: ifirst

  if(ifirst == 0)then
   ifirst = 1
   do i=1,N_det_ref_fobo
     do k=1,N_int
       psi_ref_fobo_restart(k,1,i) = psi_cas(k,1,i)
       psi_ref_fobo_restart(k,2,i) = psi_cas(k,2,i)
     enddo
   enddo
   do k=1,N_states
     do i=1,N_det_ref_fobo
       psi_ref_fobo_coef_restart(i,k) = psi_cas_coef(i,k)
     enddo
   enddo
  endif

END_PROVIDER

 BEGIN_PROVIDER [double precision, norm_psi_ref_fobo, (N_states)]
&BEGIN_PROVIDER [double precision, inv_norm_psi_ref_fobo, (N_states)]
  implicit none
  integer :: i,j
  norm_psi_ref_fobo = 0.d0
  do j = 1, N_states
   do i = 1, N_det_ref_fobo
    norm_psi_ref_fobo(j) += psi_ref_fobo_coef(i,j) * psi_ref_fobo_coef(i,j)
   enddo
   inv_norm_psi_ref_fobo(j) = 1.d0/(dsqrt(norm_psi_ref_fobo(j)))
   print *,  inv_norm_psi_ref_fobo(j)
  enddo

 END_PROVIDER

 BEGIN_PROVIDER [double precision, psi_ref_fobo_coef_interm_norm, (N_det_ref_fobo,N_states)]
  implicit none
  integer :: i,j
  do j = 1, N_states
   do i = 1, N_det_ref_fobo
    psi_ref_fobo_coef_interm_norm(i,j) = inv_norm_psi_ref_fobo(j) * psi_ref_fobo_coef(i,j)
   enddo
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [double precision, psi_non_ref_fobo_coef_interm_norm, (N_det_non_ref_fobo,N_states)]
  implicit none
  integer :: i,j
  do j = 1, N_states
   do i = 1, N_det_non_ref_fobo
    psi_non_ref_fobo_coef_interm_norm(i,j) = psi_non_ref_fobo_coef(i,j) * inv_norm_psi_ref_fobo(j)
   enddo
  enddo
 END_PROVIDER 

BEGIN_PROVIDER [integer, n_list_lmct]
 implicit none
 
 integer :: i,n_h,n_p,number_of_holes,number_of_particles
 n_list_lmct = 0
 do i = 1, N_det
  n_h = number_of_holes(psi_det(1,1,i))
  n_p = number_of_particles(psi_det(1,1,i))
  if(n_h==1.and.n_p==0.and.dabs(psi_coef(i,1)/psi_coef(1,1)).gt.thresh_coef_lmct)then
!  call debug_det(psi_det(1,1,i),N_int)
   n_list_lmct+=1 
   print*, n_h,n_p,dabs(psi_coef(i,1))
  endif
 enddo
 print*, 'n_list_lmct = ',n_list_lmct

END_PROVIDER 

BEGIN_PROVIDER [integer, list_lmct, (n_list_lmct)]
use bitmasks
 implicit none
 integer          :: exc(0:2,2,2), degree
 integer :: h1,p1,h2,p2,s1,s2
 double precision :: phase
 
 integer :: i,n_h,n_p,number_of_holes,number_of_particles,j
 j = 0
 do i = 1, N_det
  n_h = number_of_holes(psi_det(1,1,i))
  n_p = number_of_particles(psi_det(1,1,i))
  if(n_h==1.and.n_p==0.and.dabs(psi_coef(i,1)/psi_coef(1,1)).gt.thresh_coef_lmct)then
   j +=1 
   call get_excitation(psi_det(1,1,1),psi_det(1,1,i),exc,degree,phase,N_int)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
   if(list_inact_reverse(h1).ne.0)then
    list_lmct(j) = h1
   else 
    list_lmct(j) = h2
   endif
   print*, 'degree = ',degree
   print*, 'h1,p1,h2,p2',h1,p1,h2,p2
   print*, 'h1',list_lmct(j)
  endif
 enddo

END_PROVIDER 


BEGIN_PROVIDER [integer, n_list_mlct]
 implicit none
 
 integer :: i,n_h,n_p,number_of_holes,number_of_particles
 n_list_mlct = 0
 do i = 1, N_det
  n_h = number_of_holes(psi_det(1,1,i))
  n_p = number_of_particles(psi_det(1,1,i))
  if(n_h==0.and.n_p==1.and.dabs(psi_coef(i,1)/psi_coef(1,1)).gt.thresh_coef_mlct)then
   n_list_mlct+=1 
  endif
 enddo

END_PROVIDER 

BEGIN_PROVIDER [integer, list_mlct, (n_list_mlct)]
 implicit none
 integer          :: exc(0:2,2,2), degree
 integer :: h1,p1,h2,p2,s1,s2
 double precision :: phase
 
 integer :: i,n_h,n_p,number_of_holes,number_of_particles,j
 j = 0
 do i = 1, N_det
  n_h = number_of_holes(psi_det(1,1,i))
  n_p = number_of_particles(psi_det(1,1,i))
  if(n_h==0.and.n_p==1.and.dabs(psi_coef(i,1)/psi_coef(1,1)).gt.thresh_coef_mlct)then
   j +=1 
   call get_excitation(psi_det(1,1,1),psi_det(1,1,i),exc,degree,phase,N_int)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
   if(list_virt_reverse(p1).ne.0)then
    list_mlct(j) = p1
   else 
    list_mlct(j) = p2
   endif
  endif
 enddo

END_PROVIDER 
