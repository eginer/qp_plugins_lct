program pouet 
 implicit none
 integer :: ntheta,i
 double precision, allocatable  :: theta(:), r12_array(:), psi_ex(:)
 double precision :: r1,r2,r12,psimorgan,psi,psi2, a0,dtheta
 double precision, allocatable :: r1_array(:), r2_array(:,:)
 ntheta = 100
 allocate(theta(ntheta), r12_array(ntheta), psi_ex(ntheta),r1_array(3),r2_array(3,ntheta))
 dtheta = 2.d0 * dacos(-1.d0)/dble(ntheta)
 theta(1) = -dacos(-1.d0)
 do i = 2, ntheta
  theta(i) = theta(i-1) + dtheta
 enddo
 r1 = 2.00d0
 r2 = r1
 r1_array = 0.d0 
 r1_array(1) = r1
 do i = 1, ntheta
  r2_array(3,i) = r1_array(3)
  r2_array(1,i) = r1 * dcos(theta(i))
  r2_array(2,i) = r1 * dsin(theta(i))
  r12_array(i) = (r1_array(1) - r2_array(1,i))**2.d0 
  r12_array(i) = r12_array(i) + (r1_array(2) - r2_array(2,i))**2.d0 
  r12_array(i) = r12_array(i) + (r1_array(3) - r2_array(3,i))**2.d0 
  r12_array(i) = dsqrt(r12_array(i))
 enddo

 do i = 1, ntheta
  r12 = r12_array(i) 
  psi = psi2(r1,r2,r12)
  write(33,'(100(F16.10,X))')theta(i),r12,dsqrt(psi),r2_array(1,i),r2_array(2,i),r2_array(3,i)
  write(34,*)'theta_array_psi_ex(',i,',7) = dble(',theta(i),')'
  write(34,*)'array_psi_ex(',i,',7) = dble(',dsqrt(psi),')'
 enddo

 close(2)
end


      double precision function psi2(r1,r2,r12)
      implicit none
      integer j,l,m,maxci,maxcin,n,paramsread
      double precision chi, r1,r12,r2, psimorgan, scale, z
      parameter (maxci=418)
      common/frankoi/j(maxci),l(maxci),m(maxci),n(maxci),paramsread
      common/frankor/chi(maxci),scale,z
!     paramsread assumed to be initialized to !=1
      if(paramsread.ne.1) then
        call getparams
        paramsread=1
      endif

      psi2=psimorgan(r1,r2,r12)**2

      end

      subroutine getparams
      implicit none
      integer i,j,l,m,maxci,maxcin,n,paramsread
      double precision chi, scale, z
      parameter (maxci=418)
      common/frankoi/j(maxci),l(maxci),m(maxci),n(maxci),paramsread
      common/frankor/chi(maxci),scale,z

        open( 14,file='morgan.418')
        read (14,*) maxcin,z,scale
        if(maxcin.ne.maxci) stop 'maxci'
        do 108 i=1,maxci
          read(14,*) chi(i),n(i),l(i),m(i),j(i)
108     continue
        close(14)
      end

      double precision function psimorgan(r1,r2,r12)
      parameter( maxci=418 )   
      double precision chi,r1,r2,r12,dpsi,s,t,u,scale,vl(maxci),vm(maxci),vj(maxci),z
! old version without chi     double precision r1,r2,r12,dpsi,s,t,u,
!     >                 
      integer n,l,m,j,paramsread   

      common/frankoi/j(maxci),l(maxci),m(maxci),n(maxci),paramsread
      common/frankor/chi(maxci),scale,z
!      print*,'*****'
!      print*,r1,r2,r12
!      print*,z,scale
      s=2.d0*z*scale*(r1+r2)
      t=2.d0*z*scale*(r1-r2)
      u=2.d0*z*scale*r12
!      print*,s,t,u
!      print*,'*****'
!      pause
      if(s.eq.0.d0) then 
      dpsi=chi(1)      

      else

      if(t.eq.0.d0) then
      do 101 i=1,maxci
        if(l(i).eq.0) then
        vl(i)=1.d0            
        else
        vl(i)=0.d0
        end if
 101  continue
      else     
      do 102 i=1,maxci
        vl(i)=t**l(i)         
 102  continue
      end if

      if(u.eq.0.d0) then
      do 103 i=1,maxci
        if(m(i).eq.0) then
        vm(i)=1.d0            
        else
        vm(i)=0.d0
        end if
 103  continue
      else     
      do 104 i=1,maxci
        vm(i)=u**m(i)         
 104  continue
      end if

      if(s.eq.1.d0) then
      do 105 i=1,maxci
        if(j(i).eq.0) then
        vj(i)=1.d0            
        else
        vj(i)=0.d0
        end if
 105  continue
      else     
      do 106 i=1,maxci
        vj(i)=dlog(s)**j(i)         
 106  continue
      end if

      dpsi=0.d0
      do 100 i=1,maxci
        dpsi=dpsi+chi(i)*vl(i)*vm(i)*vj(i)*s**n(i)
100   continue 
      dpsi=dpsi*dexp(-0.5d0*s)                  
      end if                         
       
     psimorgan=dpsi*3.687695640501d0
!      psimorgan=dpsi
      return
      end

