
BEGIN_PROVIDER [double precision, rho_alpha_hf_ao_in_r_per_atom, (n_pts_max_per_atom,ao_num,nucl_num) ]
 implicit none
 integer :: ipoint, j,k,inucl
 rho_alpha_hf_ao_in_r_per_atom = 0.d0
 do inucl = 1, nucl_num
  do j = 1, ao_num
   do ipoint = 1, n_pts_per_atom(inucl)
    do k = 1, ao_num
     rho_alpha_hf_ao_in_r_per_atom(ipoint,j,inucl) += rho_alpha_hf_ao(k,j) * aos_in_r_array_per_atom(k,ipoint,inucl)
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 



BEGIN_PROVIDER [double precision, rho_beta_hf_ao_in_r_per_atom, (n_pts_max_per_atom,ao_num,nucl_num) ]
 implicit none
 integer :: ipoint, j,k,inucl
 rho_beta_hf_ao_in_r_per_atom = 0.d0
 do inucl = 1, nucl_num
  do j = 1, ao_num
   do ipoint = 1, n_pts_per_atom(inucl)
    do k = 1, ao_num
     rho_beta_hf_ao_in_r_per_atom(ipoint,j,inucl) += rho_beta_hf_ao(k,j) * aos_in_r_array_per_atom(k,ipoint,inucl)
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, full_dens_ao_in_r_per_atom, (n_pts_max_per_atom,ao_num,nucl_num) ]
 implicit none
 integer :: ipoint, j,k,inucl
 full_dens_ao_in_r_per_atom = 0.d0
 do inucl = 1, nucl_num
  do j = 1, ao_num
   do ipoint = 1, n_pts_per_atom(inucl)
    do k = 1, ao_num
     full_dens_ao_in_r_per_atom(ipoint,j,inucl) += full_dens_ao(k,j) * aos_in_r_array_per_atom(k,ipoint,inucl)
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 


 BEGIN_PROVIDER [integer, list_good_pairs_density_per_atom, (2,ao_num**2,nucl_num)]
&BEGIN_PROVIDER [integer, n_good_pairs_density_per_atom, (nucl_num)]
 implicit none
 integer :: ipoint, j,k,inucl
 double precision :: accu , thresh,norm(ao_num)
 double precision, allocatable :: ovmat(:,:)
 allocate(ovmat(ao_num,ao_num))
 thresh = thresh_pairs_mu_of_r
 n_good_pairs_density_per_atom = 0
 ovmat = 0.d0
 do inucl = 1, nucl_num
 !print*,''
 !print*,'inucl = ',inucl
 !print*,''
  do j = 1, ao_num
   do k = 1, ao_num
    accu = 0.d0
    if (ao_overlap_abs(k,j) < thresh) cycle
    do ipoint = 1, n_pts_per_atom(inucl)
     accu += dabs(full_dens_ao_in_r_per_atom(ipoint,j,inucl) * rho_alpha_hf_ao_in_r_per_atom(ipoint,k,inucl)) * final_weight_at_r_vector_per_atom(ipoint,inucl)
    enddo
    ovmat(k,j) = accu
   enddo
  enddo
  do j = 1, ao_num
   norm(j) = 1.d0/dsqrt(ovmat(j,j))
  enddo
  integer :: igood_at
  igood_at = 0
  do j = 1, ao_num
   integer :: igood
   igood = 0
   do k = 1, ao_num
!   ovmat(k,j) *= norm(j) * norm(k)
    if(ovmat(k,j).lt.thresh)cycle
    igood += 1
    n_good_pairs_density_per_atom(inucl) += 1
    list_good_pairs_density_per_atom(1,n_good_pairs_density_per_atom(inucl),inucl) = j
    list_good_pairs_density_per_atom(2,n_good_pairs_density_per_atom(inucl),inucl) = k
   enddo
  !igood_at += igood
  !print*,'igood = ',igood
  !write(*,'(1000(F10.5,X))')ovmat(j,:)
  !print*,'??????'
  enddo
 !print*,'n_good_pairs_density_per_atom(inucl) = ',n_good_pairs_density_per_atom(inucl)
 !print*,'igood_at                             = ',igood_at
 enddo


END_PROVIDER 

BEGIN_PROVIDER [double precision, f_hf_ab_ao_per_atom, (n_pts_max_per_atom,nucl_num)]
 implicit none
 integer :: alph,bet,gam,delt,i,ipoint,jpoint,jlist,ilist,inucl
 double precision :: thresh, integral 
 integer :: sze,sze_max,non_zero_int,m
 integer, allocatable :: out_val_index(:,:),iorder(:)
 double precision, allocatable :: out_val(:),alpha_dens(:),vcd_mat(:,:),v_eig_vectors(:,:), v_eig_values(:),sort_eigval(:)
 double precision :: wall1, wall2
 sze = ao_num
 sze_max = ao_num*ao_num
 thresh = thresh_int_mu_of_r

 f_hf_ab_ao_per_atom = 0.d0
 provide n_good_pairs_density_per_atom rho_alpha_hf_ao_in_r_per_atom rho_beta_hf_ao_in_r_per_atom full_dens_ao_in_r_per_atom
 print*,'all small stuff provided !'

 double precision :: tmp1,tmp2,tmp3,accu1,accu2,accu3
 
 double precision :: npairs_av
 call wall_time(wall1)
 allocate(out_val(sze_max),out_val_index(2,sze_max),alpha_dens(n_pts_max_per_atom),iorder(ao_num))
 allocate(vcd_mat(ao_num,ao_num), v_eig_vectors(ao_num,ao_num), v_eig_values(ao_num),sort_eigval(ao_num))
 !  First alpha pair of AO
 do inucl = 1, nucl_num
  print*,'inucl = ',inucl
  print*,'n good pairs = ',n_good_pairs_density_per_atom(inucl)
  print*,'n_pts_per_atom(inucl) = ',n_pts_per_atom(inucl)
  npairs_av = 0.d0
  do ilist = 1, n_good_pairs_density_per_atom(inucl)
   print*,'ilist = ',ilist
   alph = list_good_pairs_density_per_atom(1,ilist,inucl)
   gam = list_good_pairs_density_per_atom(2,ilist,inucl)
!  call get_ao_two_e_integrals_non_zero_jl_from_list(alph,gam,thresh,list_good_pairs_density_per_atom(1,1,inucl),n_good_pairs_density_per_atom(inucl),sze_max,out_val,out_val_index,non_zero_int)
   call get_ao_two_e_integrals_non_zero_jl(alph,gam,thresh,sze_max,sze,out_val,out_val_index,non_zero_int)
   vcd_mat = 0.d0
   do jlist = 1, non_zero_int
    bet  = out_val_index(1,jlist)
    delt = out_val_index(2,jlist)
    vcd_mat(bet,delt) = out_val(jlist)
   enddo
  !do jlist = 1, ao_num
  ! vcd_mat(jlist,jlist) = 0.d0
  !enddo
   call lapack_diagd(v_eig_values,v_eig_vectors,vcd_mat,ao_num,ao_num)
   do jlist = 1, ao_num
    sort_eigval(jlist) = dabs(v_eig_values(jlist))
    iorder(jlist) = jlist
   enddo
   call dsort(sort_eigval,iorder,ao_num)
   print*,''
   do jlist = 1, ao_num
    if(sort_eigval(jlist).lt.thresh)cycle
    print*,'v_eig_values = ',jlist,sort_eigval(jlist),v_eig_values(iorder(jlist))
   enddo
   cycle
   npairs_av += dble(non_zero_int)
   print*,non_zero_int
   do jlist = 1, non_zero_int
    bet  = out_val_index(1,jlist)
    delt = out_val_index(2,jlist)
    do ipoint = 1, n_pts_per_atom(inucl)
     f_hf_ab_ao_per_atom(ipoint,inucl) += rho_alpha_hf_ao_in_r_per_atom(ipoint,alph,inucl) * rho_beta_hf_ao_in_r_per_atom(ipoint,bet,inucl) * full_dens_ao_in_r_per_atom(ipoint,gam,inucl) * full_dens_ao_in_r_per_atom(ipoint,delt,inucl) * out_val(jlist) 
    enddo
   enddo
  enddo
  print*,'npairs av = ', npairs_av/dble(n_good_pairs_density_per_atom(inucl))
 enddo
 call wall_time(wall2)


 print*,'f_hf_ab_ao_per_atom provided in ',wall2 - wall1

END_PROVIDER 


BEGIN_PROVIDER [double precision, f_hf_ab_ao_per_atom_vector, (n_points_final_grid)]
 implicit none
 integer :: i,j,k,index_vector,inucl,ipoint
 do inucl = 1, nucl_num
  do ipoint = 1, n_pts_per_atom(inucl)
   k = index_final_points_per_atom(1,ipoint,inucl)
   j = index_final_points_per_atom(2,ipoint,inucl)
   i = index_final_points_per_atom(3,ipoint,inucl)
   index_vector = index_final_points_reverse(k,j,i)
   f_hf_ab_ao_per_atom_vector(index_vector) = f_hf_ab_ao_per_atom(ipoint,inucl)
  enddo
 enddo
END_PROVIDER 
