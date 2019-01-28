BEGIN_PROVIDER [double precision, metal_metal_dist, (n_metal,n_metal)]
 implicit none
 BEGIN_DOC
! metal_metal_dist(j,i) = distance between metal atom i and metal atom j
 END_DOC
 integer :: i,j,index_i,index_j
 double precision :: rx,ry,rz,d
 metal_metal_dist = 0.d0
 do i = 1, n_metal
  index_i = index_metal(i)
  do j = 1, n_metal
   index_j = index_metal(j)
   if(index_i == index_j)cycle
   rx = nucl_coord_transp(1,index_i) - nucl_coord_transp(1,index_j)
   ry = nucl_coord_transp(2,index_i) - nucl_coord_transp(2,index_j)
   rz = nucl_coord_transp(3,index_i) - nucl_coord_transp(3,index_j)
   d = dsqrt(rz**2 + ry**2 + rx**2)
   metal_metal_dist(i,j) = d 
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, metal_2_distance]
 implicit none
 BEGIN_DOC
! distance between the two first metal atoms 
 END_DOC
 metal_2_distance = metal_metal_dist(1,2)
END_PROVIDER 

BEGIN_PROVIDER [double precision, metal_non_metal_distance, (n_metal, n_non_metal)]
 implicit none
 BEGIN_DOC
 ! metal_non_metal_distance(i,j) = distance between metal atom i  and non metal atom j
 END_DOC
 integer :: i,j,index_i,index_j
 double precision :: rx,ry,rz,d
 do i = 1, n_metal
  index_i = index_metal(i)
  do j = 1, n_non_metal
   index_j = index_non_metal(j)
   rx = nucl_coord_transp(1,index_i) - nucl_coord_transp(1,index_j)
   ry = nucl_coord_transp(2,index_i) - nucl_coord_transp(2,index_j)
   rz = nucl_coord_transp(3,index_i) - nucl_coord_transp(3,index_j)
   d = dsqrt(rz**2 + ry**2 + rx**2)
   metal_non_metal_distance(i,j) = d
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [integer, n_metal_coord_sphere, (n_metal)]
&BEGIN_PROVIDER [integer, n_max_metal_coord_sphere]
 implicit none
 BEGIN_DOC
 ! n_metal_coord_sphere(i) = number of coordination atoms for the metal atom  i
 END_DOC
 integer :: i,j
 n_metal_coord_sphere = 0
 do i = 1, n_metal
  do j = 1, n_non_metal 
   if(metal_non_metal_distance(i,j) .le. metal_2_distance  )then
    n_metal_coord_sphere(i) += 1
   endif
  enddo
 enddo
 n_max_metal_coord_sphere = maxval(n_metal_coord_sphere) 
END_PROVIDER 

BEGIN_PROVIDER [integer, index_metal_coord_sphere, (n_metal,n_max_metal_coord_sphere)]
 implicit none
 BEGIN_DOC
 ! n_metal_coord_sphere(i,j) = index of coordinated non metal atom j to the metal atom i
 END_DOC
 integer :: i,j
 integer :: n_metal_coord_sphere_tmp(n_metal,n_max_metal_coord_sphere)
 n_metal_coord_sphere_tmp = 0
 do i = 1, n_metal
  do j = 1, n_non_metal 
   if(metal_non_metal_distance(i,j) .le. metal_2_distance  )then
    n_metal_coord_sphere_tmp(i,j) += 1
    index_metal_coord_sphere(i,n_metal_coord_sphere_tmp(i,j)) = j 
   endif
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [integer, n_coord_non_metal, (n_non_metal)]
&BEGIN_PROVIDER [integer, coord_non_metal, (n_non_metal,n_metal)]
 implicit none
 BEGIN_DOC
 ! n_coord_non_metal(i) = number of coordinated metal atoms for a non metal atom i
 END_DOC
 integer :: i,j,imin
 double precision :: dmin
 n_coord_non_metal = 0
 do i = 1, n_metal
  do j = 1, n_non_metal 
   if(metal_non_metal_distance(i,j) .le. metal_2_distance  )then
    n_coord_non_metal(j) += 1
    coord_non_metal(j,n_coord_non_metal(j)) = i
   endif
  enddo
 enddo
 do j = 1, n_non_metal
  if(n_coord_non_metal(j) .ne. 1 .and. n_coord_non_metal(j) .ne. 2)then
   print*,'PB !!!!!!!!'
   print*,'j = ',j
   print*,'n_coord_non_metal(j)',n_coord_non_metal(j)
   stop
  endif
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [integer, index_other_ext, (n_external_ligand)]
 implicit none
 BEGIN_DOC
 ! gives the other external ligand associated for a given external ligand 
 END_DOC
 integer :: i,index_i,index_final_i,metal_atom_i
 integer :: j,index_j,index_final_j,metal_atom_j
 index_other_ext = 0
 do i = 1, n_external_ligand 
  index_i       = index_external_ligand(i)    ! index of non metal atom 
  metal_atom_i  = coord_non_metal(index_i,1)  ! index of connected metal atom 
  do j = 1, n_external_ligand 
   if(i==j)cycle
   index_j       = index_external_ligand(i)    ! index of non metal atom 
   metal_atom_j  = coord_non_metal(index_j,1)  ! index of connected metal atom 
   if(metal_atom_j == metal_atom_i)then
    index_other_ext(i) = j 
   endif
  enddo 
 enddo
 

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, metal_barycentre, (3)]
&BEGIN_PROVIDER [double precision, metal_difference, (3)]
 implicit none
 integer :: index_1,index_2
 index_1 = index_metal(1)
 index_2 = index_metal(2)
 metal_barycentre(1) = nucl_coord_transp(1,index_1) + nucl_coord_transp(1,index_2)
 metal_barycentre(2) = nucl_coord_transp(2,index_1) + nucl_coord_transp(2,index_2)
 metal_barycentre(3) = nucl_coord_transp(3,index_1) + nucl_coord_transp(3,index_2)
 metal_barycentre = metal_barycentre * 0.5d0 

 metal_difference(1) = nucl_coord_transp(1,index_1) - nucl_coord_transp(1,index_2)
 metal_difference(2) = nucl_coord_transp(2,index_1) - nucl_coord_transp(2,index_2)
 metal_difference(3) = nucl_coord_transp(3,index_1) - nucl_coord_transp(3,index_2)

END_PROVIDER 
