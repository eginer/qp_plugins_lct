program test
 read_wf = .True.
 touch read_wf
! call test_spherical_grid
 call routine_bis
end
subroutine routine
 implicit none
 integer :: i,j,k,i_point
 double precision :: r1(3),r12
 double precision :: two_dm_hf,two_dm_hf_laplacian,total_hf_dm
 double precision :: two_dm_psi,two_dm_psi_laplacian,total_psi_dm

 double precision :: spher_av_n2_hartree,num_sphe_hf
 double precision :: num_sphe_psi
 r12 = 0.00001d0
 do i_point = 1, n_points_final_grid
  r1(1) = final_grid_points(1,i_point)
  r1(2) = final_grid_points(2,i_point)
  r1(3) = final_grid_points(3,i_point)
  call spherical_averaged_two_dm_HF_at_second_order(r1,r12,two_dm_hf,two_dm_hf_laplacian,total_hf_dm)
  num_sphe_hf = spher_av_n2_hartree(r1,r12)
!  if(dabs(num_sphe_hf).gt.1.d-10)then
!   if(dabs( num_sphe_hf - total_hf_dm).gt.1.d-10)then
!    print*,num_sphe_hf,total_hf_dm,num_sphe_hf/total_hf_dm
!   endif
!  endif
  call spherical_averaged_two_dm_at_second_order(r1,r12,1,two_dm_psi,two_dm_psi_laplacian,total_psi_dm)
  if(two_dm_psi_laplacian/two_dm_psi-two_dm_hf_laplacian/two_dm_hf .lt. 0.d0)then
   print*,''
   print*,two_dm_psi_laplacian/two_dm_psi,two_dm_hf_laplacian/two_dm_hf,two_dm_psi_laplacian/two_dm_psi-two_dm_hf_laplacian/two_dm_hf
   print*,two_dm_psi,two_dm_hf
   print*,r1
  endif
!  num_sphe_psi = spher_av_n2_psi(r1,r12)
!  if(dabs(num_sphe_psi).gt.1.d-10)then
!   if(dabs( num_sphe_psi - total_psi_dm).gt.1.d-10)then
!    print*,num_sphe_psi,total_psi_dm,num_sphe_psi/total_psi_dm
!   endif
!  endif
 enddo

end

subroutine routine_bis
 implicit none
 double precision :: d2_hf,d2_psi
 integer :: i,j,k,i_point,nr12
 double precision :: r1(3),r12,dr12,r12max
 double precision :: two_dm_hf,two_dm_hf_laplacian,total_hf_dm
 double precision :: two_dm_psi,two_dm_psi_laplacian,total_psi_dm
 double precision :: j_r12,numerator

 double precision :: spher_av_n2_hartree,num_sphe_hf,spher_av_n2_hf
 double precision :: num_sphe_psi

 double precision :: num_sphe_hf_r,num_sphe_hf_l
 double precision :: num_sphe_psi_r,num_sphe_psi_l

 double precision :: f_HF_val_ab,two_bod_dens,w_hf,mu_hf
 call f_HF_valence_ab(r1,r1,f_HF_val_ab,two_bod_dens) 
 w_hf  = f_HF_val_ab/two_bod_dens

 r1 = 0.d0

 dr12 = 0.0001d0
 r12max = 2.d0
 r12 = 0.00001d0
 r1(1) = 1.2694538767039474d-2
 nr12 = int(r12max/dr12)
 call spherical_averaged_two_dm_HF_at_second_order(r1,r12  ,two_dm_hf ,two_dm_hf_laplacian ,total_hf_dm)
 call spherical_averaged_two_dm_at_second_order(   r1,r12,1,two_dm_psi,two_dm_psi_laplacian,total_psi_dm)
 numerator = two_dm_psi_laplacian * two_dm_hf - two_dm_hf_laplacian * two_dm_psi
 print*,'numerator = ',numerator
 print*,'lapl      = ',two_dm_psi_laplacian,two_dm_hf_laplacian
 print*,'ontop     = ',two_dm_psi,two_dm_hf
 print*,'num/n2hf  = ',numerator/(two_dm_hf)**2.d0

 print*,'nr12 = ',nr12
 print*,two_dm_psi_laplacian/two_dm_psi,two_dm_hf_laplacian/two_dm_hf,two_dm_psi_laplacian/two_dm_psi-two_dm_hf_laplacian/two_dm_hf
 print*,two_dm_psi,two_dm_hf
 print*,r1


 print*,''
 print*,'HF'
 print*,''
 do i = 2, 6
  dr12 = 10.d0**(-dble(i))
 
  r12 = dr12
  num_sphe_hf  = spher_av_n2_hf(r1,r12)
 
  r12 = 2.d0 * dr12
  num_sphe_hf_r  = spher_av_n2_hf(r1,r12)
 
  r12 = 0.d0
  num_sphe_hf_l  = spher_av_n2_hf(r1,r12)
  d2_hf = num_sphe_hf_r -2.d0 * num_sphe_hf + num_sphe_hf_l
  d2_hf = d2_hf/(dr12)**2.d0
  print*,'d2_hf  = ',d2_hf,two_dm_hf_laplacian,dr12
 enddo

 print*,''
 print*,'PSI'
 print*,''
 do i = 2, 6
  dr12 = 10.d0**(-dble(i))
 
  r12 = dr12
  call spher_av_n2_psi(r1,r12,num_sphe_psi,j_r12)
 
  r12 = 2.d0 * dr12
  call spher_av_n2_psi(r1,r12,num_sphe_psi_r,j_r12)
 
  r12 = 0.d0
  call spher_av_n2_psi(r1,r12,num_sphe_psi_l,j_r12)

  d2_psi = num_sphe_psi_r -2.d0 * num_sphe_psi + num_sphe_psi_l
  d2_psi = d2_psi/(dr12)**2.d0
  print*,'d2_psi = ',d2_psi,two_dm_psi_laplacian,dr12
 enddo
  print*,'d2_j   = ',(d2_psi/two_dm_psi - d2_hf/two_dm_hf)
 print*,'mu = ', 0.33333d0 * 1.5d0 * dsqrt(dacos(-1.d0)) * (d2_psi/two_dm_psi - d2_hf/two_dm_hf)
 mu_hf =  w_hf * dsqrt(dacos(-1.d0)) * 0.5d0
 print*,'mu_hf',mu_hf
 stop

 write(34,*)' # r12                n2hf               n2fci          n2fci/n2hf       j_r12            n2hf             n2fci'
 do i = 1, nr12
  num_sphe_hf  = spher_av_n2_hf(r1,r12)
  call spher_av_n2_psi(r1,r12,num_sphe_psi,j_r12)
  write(34,'(100(F16.10,X))')r12, num_sphe_hf , num_sphe_psi ,num_sphe_psi/num_sphe_hf, j_r12,two_dm_hf + 0.5d0 * two_dm_hf_laplacian * r12**2, two_dm_psi + 0.5d0 * r12**2
  r12 += dr12
 enddo


end
