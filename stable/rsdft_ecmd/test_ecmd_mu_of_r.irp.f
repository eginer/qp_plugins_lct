program two_body_dm
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
 read_wf = .True.
 touch read_wf
 
!call comparaison_mu_of_r
 call tracer_interaction
end


 subroutine comparaison_mu_of_r 
 implicit none
 integer :: i
 double precision :: mu_dif
 mu_dif = 0.d0 
 do i = 1, n_points_final_grid
  mu_dif += dabs(mu_of_r_hf_coal_vector_ecmd(i)- mu_of_r_hf_coal_vector(i)) * final_weight_at_r_vector(i) 
 enddo
 print*,'mu_dif int = ', mu_dif

 end



 subroutine tracer_interaction
 implicit none
 integer :: i,j,nr
 double precision :: r12, dr12,r12max,mu
 double precision :: r1(3),r2(3),local_potential,two_bod
 nr = 100 
 r12max = 2.d0
 dr12 = r12max/dble(nr)
 r1(1) = 0.0d0
 r1(2) = 0.0d0
 r1(3) = 0.0d0
 print*,'r1 = ',r1
!!!!!!!!! definition of mu
 r2 = r1 
 call f_HF_ab_ecmd(r1,r2,local_potential,two_bod)
 local_potential = local_potential /  two_bod
 mu =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 print*,'muuuuuuuuuuuuuuuuuuuu = ',mu
!!!!!! tracer interaction

 r12 = 0.d0
 do i = 1, nr
  r12 += dr12 
  r2(1) += dr12
  call f_HF_ab_ecmd(r1,r2,local_potential,two_bod)
  write(33,'(10(F16.10,X))')r12,local_potential/two_bod,erf(mu * r12)/r12,local_potential,two_bod
 enddo
 end
