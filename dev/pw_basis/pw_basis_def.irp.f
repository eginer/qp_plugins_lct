BEGIN_PROVIDER [double precision, vol_box_pw]
 BEGIN_DOC
! Volume of the box on which the plane waves are defined
 END_DOC
 implicit none
 vol_box_pw = pw_box_x * pw_box_y * pw_box_z
END_PROVIDER 

BEGIN_PROVIDER [integer, n_pw_array, (3)]
 implicit none
 BEGIN_DOC
! number of plane waves in each dimension stored in an array
 END_DOC
 n_pw_array(1) = pw_nx
 n_pw_array(2) = pw_ny
 n_pw_array(3) = pw_nz
END_PROVIDER 

BEGIN_PROVIDER [double precision, delta_k_pw, (3)]
 implicit none
 BEGIN_DOC
! increment of plane wave argument  in real space 
 END_DOC
 include 'constants.include.F'
 delta_k_pw(1) =  2.d0 * pi / pw_box_x
 delta_k_pw(2) =  2.d0 * pi / pw_box_y
 delta_k_pw(3) =  2.d0 * pi / pw_box_z
END_PROVIDER 


 BEGIN_PROVIDER [double precision, norm_box_pw, (3)]
&BEGIN_PROVIDER [double precision, norm_tot_pw ]
 BEGIN_DOC
! norm of plane wave in each direction and total norm of the plane wave 
 END_DOC
 implicit none
 norm_box_pw(1) = 1.d0/dsqrt(pw_box_x)
 norm_box_pw(2) = 1.d0/dsqrt(pw_box_y)
 norm_box_pw(3) = 1.d0/dsqrt(pw_box_z)
 norm_tot_pw    = norm_box_pw(1) * norm_box_pw(2) * norm_box_pw(3)
END_PROVIDER 

 BEGIN_PROVIDER [integer, n_pw_basis]
&BEGIN_PROVIDER [integer, n_max_pw_basis]
 implicit none
 BEGIN_DOC
! total number of plane waves, and maximum number of plane waves in all dimension
 END_DOC
 n_pw_basis = pw_nx * pw_ny * pw_nz
 print*,'n_pw_basis = ',n_pw_basis
 n_max_pw_basis = max(pw_nx , pw_ny)
 n_max_pw_basis = max(n_max_pw_basis,pw_nz)
END_PROVIDER

 BEGIN_PROVIDER [integer, list_pw, (3,n_pw_basis)]
&BEGIN_PROVIDER [integer, list_reverse_pw, (n_max_pw_basis,n_max_pw_basis,n_max_pw_basis)]
&BEGIN_PROVIDER [double precision, d_list_pw, (3,n_pw_basis)]
 integer :: i,j,k,ii
 BEGIN_DOC
! list of n_x,n_y,n_z for each plane wave, and reverse list
! By convention, we chose not to take the k_x = k_y = k_z = 0 plane wave as it is just a constant function
 END_DOC
 ii = 0
 do i = 1, pw_nx
  do j = 1, pw_ny
   do k = 1, pw_nz
    ii += 1
    list_pw(1,ii) = i
    list_pw(2,ii) = j
    list_pw(3,ii) = k
    list_reverse_pw(i,j,k) = ii
    d_list_pw(1,ii) = dble(i)
    d_list_pw(2,ii) = dble(j)
    d_list_pw(3,ii) = dble(k)
   enddo
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [integer, n_delta_k_uniq, (3) ]
&BEGIN_PROVIDER [integer, n_max_delta_k_uniq ]
 implicit none
 BEGIN_DOC
! Number of uniq |n_1 - n_2| for each dimension of plane waves
! This is useful for computing integrals of operators between plane waves
 END_DOC
 integer :: m
 do m = 1, 3
  n_delta_k_uniq(m) = n_pw_array(m)
 enddo
 n_max_delta_k_uniq = maxval(n_delta_k_uniq)

END_PROVIDER 

BEGIN_PROVIDER [integer, int_delta_k_uniq, (n_max_delta_k_uniq,3)]
 implicit none
 integer :: m,i
 do m = 1, 3
  do i = 1, n_pw_array(m)
   int_delta_k_uniq(i,m) = i - 1
  enddo
 enddo
END_PROVIDER 
