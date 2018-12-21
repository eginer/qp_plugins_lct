!!!!!!!!!!first term!!!!!!!!

 double precision function d0delta(rs,xi)
 BEGIN_DOC
  ! eq 49 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi
 d0delta=(0.70605d0+0.12927d0*xi**2)*rs
 return
 end


 double precision function d_d0delta(xi)
 BEGIN_DOC
  ! Derivative of d0delta with respect to rs
 END_DOC
 implicit none
 double precision rs,xi
 d_d0delta=(0.70605d0+0.12927d0*xi**2)
 return
 end

 double precision function d_xi_d0delta(rs,xi)
 BEGIN_DOC
  ! Derivative of d0delta with respect with xi
 END_DOC
 implicit none
 double precision rs,xi
 d_xi_d0delta=(2d0*0.12927d0*xi)*rs
 return
 end


 double precision function denominator_delta(rs,xi,mu)
 BEGIN_DOC
  ! Denominator of DELTA_{LR-SR} (denominator of eq 42 of Panziani et al, PRB 73, 155111 (2006)) 
 END_DOC
 implicit none
 double precision rs,xi,mu,d0delta
 denominator_delta=(1d0+d0delta(rs,xi)**2*mu**2)**4d0
 return
 end 


 double precision function d_denominator_delta(rs,xi,mu)
 BEGIN_DOC
  ! Derivative of denominator_delta with respect to rs  
 END_DOC
 implicit none
 double precision rs,xi,mu,d0delta,d_d0delta
 d_denominator_delta=8d0*mu**2*d0delta(rs,xi)*d_d0delta(xi)*(1+d0delta(rs,xi)**2*mu**2)**3
 return
 end

 double precision function d_xi_denominator_delta(rs,xi,mu)
 BEGIN_DOC
  ! Derivative of denominator_delta with respect to xi  
 END_DOC
 implicit none
 double precision rs,xi,mu,d0delta,d_xi_d0delta
 d_xi_denominator_delta=8d0*mu**2*d0delta(rs,xi)*d_xi_d0delta(rs,xi)*(1+d0delta(rs,xi)**2*mu**2)**3
 return
 end

 double precision function delta_2(rs)
 BEGIN_DOC
  ! eq 48 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs
 delta_2=0.073867d0*rs**(3d0/2d0)
 return
 end 

 double precision function d_delta_2(rs)
 BEGIN_DOC
  ! Derivative of delta_2 with respecto to rs
 END_DOC
 implicit none
 double precision rs
 d_delta_2=0.1108005d0*rs**(1d0/2d0)
 return
 end 

 double precision function d_1st_deltaterm(rs,mu,xi)
 BEGIN_DOC
  ! Derivative of delta_2*mu**2/(1+d0**2*mu**2)**4 with respect to rs (1st part of equation 42 of Panziani et al, PRB 73, 155111 (2006))
 END_DOC
 implicit none
 double precision rs,mu,xi,d_delta_2,denominator_delta,d_denominator_delta,delta_2
 d_1st_deltaterm=d_delta_2(rs)*mu**2/denominator_delta(rs,xi,mu)-d_denominator_delta(rs,xi,mu)*delta_2(rs)*mu**2/denominator_delta(rs,xi,mu)**2
 return
 end  


 double precision function d_xi_1st_deltaterm(rs,mu,xi)
 BEGIN_DOC
  ! Derivative of delta_2*mu**2/(1+d0**2*mu**2)**4 with respect to xi (1st part of equation 42 of Panziani et al, PRB 73, 155111 (2006))
 END_DOC
 implicit none
 double precision rs,mu,xi,d_delta_2,denominator_delta,d_xi_denominator_delta,delta_2
 d_xi_1st_deltaterm=-d_xi_denominator_delta(rs,xi,mu)*delta_2(rs)*mu**2/denominator_delta(rs,xi,mu)**2
 return
 end

!!!!!!!!2nd terms!!!!!


 double precision function C3t_delta(rs,xi)
 BEGIN_DOC
  ! eq 47 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,g0f,pi
 pi=dacos(-1.d0)
 C3t_delta=-(1d0-xi**2)*g0f(rs)*(2d0*sqrt(2d0)-1d0)/(2d0*sqrt(pi)*rs**3)
 return
 end


 double precision function d_C3t_delta(rs,xi)
 BEGIN_DOC
  !Derivative of C3t_delta with respect to rs
 END_DOC
 implicit none        
 double precision rs,xi,g0f,pi,g0d
 pi=dacos(-1.d0)      
 d_C3t_delta=-(1d0-xi**2)*((2d0*sqrt(2d0)-1d0)/(2d0*sqrt(pi)))*(g0d(rs)/rs**3-g0f(rs)*3d0/(rs**4))
 return
 end

 double precision function d_xi_C3t_delta(rs,xi)
 BEGIN_DOC
  !Derivative of C3t_delta with respect to xi
 END_DOC
 implicit none
 double precision rs,xi,g0f,pi,g0d
 pi=dacos(-1.d0)
 d_xi_C3t_delta = 2d0*xi*g0f(rs)*(2d0*sqrt(2d0)-1d0)/(2d0*sqrt(pi)*rs**3)
 return
 end
 

 double precision function D3_delta(rs)
 BEGIN_DOC
  ! eq 34 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs
 D3_delta=exp(-0.31d0*rs)*(-4.95d0*rs+rs**2)/(rs**3)
 return
 end 

 double precision function d_D3_delta(rs)
 BEGIN_DOC
  ! Derivative of D3_delta with respect to rs 
 END_DOC
 implicit none
 double precision rs
 d_D3_delta=exp(-0.31d0*rs)*(-4.95d0+2d0*rs)/(rs**3)-3*exp(-0.31d0*rs)*(-4.95d0*rs+rs**2)/(rs**4)-0.31d0*exp(-0.31d0*rs)*(-4.95d0*rs+rs**2)/(rs**3)
 return
 end

 double precision function c5_delta(rs,xi)
 BEGIN_DOC
  ! eq 29 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,dpol,D3_delta
 c5_delta=((1d0+xi)/2d0)**2*dpol(rs*(2d0/(1d0+xi))**(1d0/3d0))+((1d0-xi)/2d0)**2*dpol(rs*(2d0/(1d0-xi))**(1d0/3d0))+(1d0-xi**2)*D3_delta(rs)
 return
 end 

 double precision function d_c5_delta(rs,xi)
 BEGIN_DOC
  ! Derivative c5_delta with respect with rs
 END_DOC
 implicit none
 double precision rs,xi,dpold,d_D3_delta
 d_c5_delta=((1d0+xi)/2d0)**(5d0/3d0)*dpold(rs*(2d0/(1d0+xi))**(1d0/3d0))+((1d0-xi)/2d0)**(5d0/3d0)*dpold(rs*(2d0/(1d0-xi))**(1d0/3d0))+(1d0-xi**2)*d_D3_delta(rs)
 return
 end 

 double precision function d_xi_c5_delta(rs,xi)
 BEGIN_DOC
  ! Derivative c5_delta with respect with xi
 END_DOC
 implicit none
 double precision rs,xi,dpold,D3_delta,dpol,part1,part2,part3
 part1 = ((1d0+xi)/2d0)*dpol(rs*(2d0/(1d0+xi))**(1d0/3d0))-rs*(1d0+xi)**(2d0/3d0)*dpold(rs*(2d0/(1d0+xi))**(1d0/3d0))/(3d0*2d0**(5d0/3d0)) 
 part2 = -((1d0-xi)/2d0)*dpol(rs*(2d0/(1d0-xi))**(1d0/3d0))+rs*(1d0-xi)**(2d0/3d0)*dpold(rs*(2d0/(1d0-xi))**(1d0/3d0))/(3d0*2d0**(5d0/3d0)) 
 part3 = -2d0*xi*D3_delta(rs)
 d_xi_c5_delta= part1 + part2 + part3
 return
 end 

 double precision function C5t_delta(rs,xi)
 BEGIN_DOC
  ! eq 47 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,pi,c5_delta
 pi=dacos(-1.d0)
 C5t_delta=-3d0*c5_delta(rs,xi)*(3d0-sqrt(2d0))/(20d0*sqrt(2d0*pi)*rs**3)
 return
 end

 double precision function d_C5t_delta(rs,xi)
 BEGIN_DOC
  ! Derivative of C5t_delta with respect to rs
 END_DOC
 implicit none        
 double precision rs,xi,pi,d_c5_delta,c5_delta
 pi=dacos(-1.d0)                                           
 d_C5t_delta=-3d0*(3d0-sqrt(2d0))*(d_c5_delta(rs,xi)/rs**3-c5_delta(rs,xi)*3d0*rs**2d0/(rs**6))/(20d0*sqrt(2d0*pi))
 return      
 end    

 double precision function d_xi_C5t_delta(rs,xi)
 BEGIN_DOC
  ! Derivative of C5t_delta with respect to xi
 END_DOC
 implicit none
 double precision rs,xi,pi,d_xi_c5_delta
 pi=dacos(-1.d0)
 d_xi_C5t_delta=-3d0*d_xi_c5_delta(rs,xi)*(3d0-sqrt(2d0))/(20d0*sqrt(2d0*pi)*rs**3)
 return
 end

 double precision function delta_3(rs,xi)
 BEGIN_DOC
  ! eq 43 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,C3t_delta,d0delta,C5t_delta
 delta_3=4d0*C3t_delta(rs,xi)*d0delta(rs,xi)**6+C5t_delta(rs,xi)*d0delta(rs,xi)**8
 return
 end

 double precision function d_delta_3(rs,xi)
 BEGIN_DOC
  ! Derivative of delta_3 with respect to rs 
 END_DOC
 implicit none
 double precision rs,xi,d_C3t_delta,d0delta,C3t_delta,d_d0delta,C5t_delta,d_C5t_delta
 d_delta_3=4d0*d_C3t_delta(rs,xi)*d0delta(rs,xi)**6+4d0*6d0*C3t_delta(rs,xi)*d_d0delta(xi)*d0delta(rs,xi)**5+d_C5t_delta(rs,xi)*d0delta(rs,xi)**8+8d0*C5t_delta(rs,xi)*d_d0delta(xi)*d0delta(rs,xi)**7
 return
 end

 double precision function d_xi_delta_3(rs,xi)
 BEGIN_DOC            
  ! Derivative of delta_3 with respect to xi
 END_DOC
 implicit none
 double precision rs,xi,d_xi_C3t_delta,d0delta,C3t_delta,d_xi_d0delta,C5t_delta,d_xi_C5t_delta
 d_xi_delta_3=4d0*d_xi_C3t_delta(rs,xi)*d0delta(rs,xi)**6+4d0*6d0*C3t_delta(rs,xi)*d_xi_d0delta(rs,xi)*d0delta(rs,xi)**5+d_xi_C5t_delta(rs,xi)*d0delta(rs,xi)**8+8d0*C5t_delta(rs,xi)*d_xi_d0delta(rs,xi)*d0delta(rs,xi)**7
 return   
 end 

 double precision function d_2nd_deltaterm(rs,mu,xi)
 BEGIN_DOC
  ! Derivative of delta_3*mu**3/(1+d0**2*mu**2)**4 with respect to rs (2nd part of equation 42 of Panziani et al, PRB 73, 155111 (2006))
 END_DOC
 implicit none
 double precision rs,mu,xi,d_delta_3,denominator_delta,d_denominator_delta,delta_3
 d_2nd_deltaterm=d_delta_3(rs,xi)*mu**3/denominator_delta(rs,xi,mu)-d_denominator_delta(rs,xi,mu)*delta_3(rs,xi)*mu**3/denominator_delta(rs,xi,mu)**2
 return
 end 


 double precision function d_xi_2nd_deltaterm(rs,mu,xi)
 BEGIN_DOC
  ! Derivative of delta_3*mu**3/(1+d0**2*mu**2)**4 with respect to xi (2nd part of equation 42 of Panziani et al, PRB 73, 155111 (2006))
 END_DOC
 implicit none
 double precision rs,mu,xi,d_xi_delta_3,denominator_delta,d_xi_denominator_delta,delta_3
 d_xi_2nd_deltaterm=d_xi_delta_3(rs,xi)*mu**3/denominator_delta(rs,xi,mu)-d_xi_denominator_delta(rs,xi,mu)*delta_3(rs,xi)*mu**3/denominator_delta(rs,xi,mu)**2
 return
 end  

!!!!!!!!third terms!!!!!


 double precision function C2_delta(rs,xi)
 BEGIN_DOC
  ! eq 30 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,pi,g0f
 pi=dacos(-1.d0)
 C2_delta=-3d0*(1d0-xi**2)*(g0f(rs)-0.5d0)/(8d0*rs**3)
 return
 end

 double precision function d_C2_delta(rs,xi)
 BEGIN_DOC
  ! Derivative of C2_delta with respect with rs
 END_DOC
 implicit none
 double precision rs,xi,pi,g0f,g0d
 pi=dacos(-1.d0)
 d_C2_delta=-3d0*(1d0-xi**2)*(g0d(rs)/rs**3-(g0f(rs)-0.5d0)*3d0*rs**2d0/(rs**6))/8d0
 return
 end

 double precision function d_xi_C2_delta(rs,xi)
 BEGIN_DOC  
  ! Derivative of C2_delta with respect with xi
 END_DOC
 implicit none
 double precision rs,xi,pi,g0f
 pi=dacos(-1.d0)
 d_xi_C2_delta=3d0*2d0*xi*(g0f(rs)-0.5d0)/(8d0*rs**3)
 return
 end

 double precision function D2_delta(rs)
 BEGIN_DOC
  ! eq 33 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs
 D2_delta=exp(-0.547d0*rs)*(-0.388d0*rs+0.676d0*rs**2)/(rs**2)
 return
 end

 double precision function d_D2_delta(rs)
 BEGIN_DOC
  ! Derivative of D2_delta with respect to rs 
 END_DOC
 implicit none
 double precision rs
 d_D2_delta=exp(-0.547d0*rs)*(-0.388+2d0*0.676d0*rs)/(rs**2)-2d0*exp(-0.547d0*rs)*(-0.388d0*rs+0.676d0*rs**2)/(rs**3)-0.547d0*exp(-0.547d0*rs)*(-0.388d0*rs+0.676d0*rs**2)/(rs**2)
 return
 end

 double precision function phi_8_delta(xi)
 BEGIN_DOC
  ! eq 14 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision xi
 phi_8_delta=0.5d0*((1d0+xi)**(8d0/3d0)+(1d0-xi)**(8d0/3d0))
 return
 end

 double precision function d_xi_phi_8_delta(xi)
 BEGIN_DOC
  ! derivative of phi_8_delta with respect to xi
 END_DOC
 implicit none
 double precision xi
 d_xi_phi_8_delta=0.5d0*((8d0/3d0)*(1d0+xi)**(5d0/3d0)-(8d0/3d0)*(1d0-xi)**(5d0/3d0))
 return
 end

 double precision function c4_delta(rs,xi)
 BEGIN_DOC
  ! eq 28 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,dpol,D2_delta,pi,alpha,phi_8_delta
 pi=dacos(-1.d0)
 alpha=(4d0/(9d0*pi))**(1d0/3d0)
 c4_delta=((1d0+xi)/2d0)**2*dpol(rs*(2d0/(1d0+xi))**(1d0/3d0))+((1d0-xi)/2d0)**2*dpol(rs*(2d0/(1d0-xi))**(1d0/3d0))+(1d0-xi**2)*D2_delta(rs)-phi_8_delta(xi)/(5d0*alpha**2*rs**2)
 return
 end

 double precision function d_c4_delta(rs,xi)
 BEGIN_DOC
  ! Derivative of c4_delta with respect to rs
 END_DOC
 implicit none
 double precision rs,xi,pi,dpold,d_D2_delta,alpha,phi_8_delta
 pi=dacos(-1.d0)
 alpha=(4d0/(9d0*pi))**(1d0/3d0)
 d_c4_delta=((1d0+xi)/2d0)**(5d0/3d0)*dpold(rs*(2d0/(1d0+xi))**(1d0/3d0))+((1d0-xi)/2d0)**(5d0/3d0)*dpold(rs*(2d0/(1d0-xi))**(1d0/3d0))+(1d0-xi**2)*d_D2_delta(rs)+2d0*phi_8_delta(xi)/(5d0*alpha**2*rs**3)
 return
 end

 double precision function d_xi_c4_delta(rs,xi)
 BEGIN_DOC
  ! Derivative of c4_delta with respect to xi
 END_DOC
 implicit none
 double precision rs,xi,pi,dpold,dpol,D2_delta,alpha,d_xi_phi_8_delta,part1,part2,part3
 pi=dacos(-1.d0)
 alpha=(4d0/(9d0*pi))**(1d0/3d0)
 part1 = ((1d0+xi)/2d0)*dpol(rs*(2d0/(1d0+xi))**(1d0/3d0))-rs*(1d0+xi)**(2d0/3d0)*dpold(rs*(2d0/(1d0+xi))**(1d0/3d0))/(3d0*2d0**(5d0/3d0))
 part2 = -((1d0-xi)/2d0)*dpol(rs*(2d0/(1d0-xi))**(1d0/3d0))+rs*(1d0-xi)**(2d0/3d0)*dpold(rs*(2d0/(1d0-xi))**(1d0/3d0))/(3d0*2d0**(5d0/3d0))
 part3 = -2d0*xi*D2_delta(rs)-d_xi_phi_8_delta(xi)/(5d0*alpha**2*rs**2)
 d_xi_c4_delta= part1 + part2 + part3
 return
 end

 double precision function capital_c4_delta(rs,xi)
 BEGIN_DOC
  ! eq 30 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,c4_delta
 capital_c4_delta=-9d0*c4_delta(rs,xi)/(64d0*rs**3)
 return
 end

 double precision function d_capital_c4_delta(rs,xi)
 BEGIN_DOC
  ! Derivative of capital_c4_delta with respect to rs
 END_DOC
 implicit none
 double precision rs,xi,pi,c4_delta,d_c4_delta
 d_capital_c4_delta=(-9d0/64d0)*(d_c4_delta(rs,xi)/rs**3-3*c4_delta(rs,xi)/rs**4)
 return
 end

 double precision function d_xi_capital_c4_delta(rs,xi)
 BEGIN_DOC  
  ! Derivative of capital_c4_delta with respect to xi
 END_DOC
 implicit none
 double precision rs,xi,pi,c4_delta,d_xi_c4_delta
 d_xi_capital_c4_delta =-9d0*d_xi_c4_delta(rs,xi)/(64d0*rs**3) 
 return
 end

double precision function delta_4(rs,xi)
 BEGIN_DOC
  ! eq 44 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,d0delta,C2_delta,capital_c4_delta
 delta_4=4d0*C2_delta(rs,xi)*d0delta(rs,xi)**6+capital_c4_delta(rs,xi)*d0delta(rs,xi)**8
 return
 end

 double precision function d_delta_4(rs,xi)
 BEGIN_DOC
  ! Derivative of delta_4 with respect to rs 
 END_DOC
 implicit none
 double precision rs,xi,d_C3t_delta,d0delta,C2_delta,d_C2_delta,d_d0delta,capital_c4_delta,d_capital_c4_delta
 d_delta_4=4d0*d_C2_delta(rs,xi)*d0delta(rs,xi)**6+4d0*6d0*C2_delta(rs,xi)*d_d0delta(xi)*d0delta(rs,xi)**5+d_capital_c4_delta(rs,xi)*d0delta(rs,xi)**8+8d0*capital_c4_delta(rs,xi)*d_d0delta(xi)*d0delta(rs,xi)**7
 return
 end

 double precision function d_xi_delta_4(rs,xi)
 BEGIN_DOC
  ! Derivative of delta_4 with respect to xi 
 END_DOC
 implicit none
 double precision rs,xi,d_xi_C3t_delta,d0delta,C2_delta,d_xi_C2_delta,d_xi_d0delta,capital_c4_delta,d_xi_capital_c4_delta
 d_xi_delta_4=4d0*d_xi_C2_delta(rs,xi)*d0delta(rs,xi)**6+4d0*6d0*C2_delta(rs,xi)*d_xi_d0delta(rs,xi)*d0delta(rs,xi)**5+d_xi_capital_c4_delta(rs,xi)*d0delta(rs,xi)**8+8d0*capital_c4_delta(rs,xi)*d_xi_d0delta(rs,xi)*d0delta(rs,xi)**7
 return
 end

 double precision function d_3rd_deltaterm(rs,mu,xi)
 BEGIN_DOC
  ! Derivative of delta_4*mu**4/(1+d0**2*mu**2)**4 with respect to rs (3rd part of equation 42 of Panziani et al, PRB 73, 155111 (2006))
 END_DOC
 implicit none
 double precision rs,mu,xi,d_delta_4,denominator_delta,d_denominator_delta,delta_4
 d_3rd_deltaterm=d_delta_4(rs,xi)*mu**4/denominator_delta(rs,xi,mu)-d_denominator_delta(rs,xi,mu)*delta_4(rs,xi)*mu**4/denominator_delta(rs,xi,mu)**2
 return
 end

 double precision function d_xi_3rd_deltaterm(rs,mu,xi)
 BEGIN_DOC
  ! Derivative of delta_4*mu**4/(1+d0**2*mu**2)**4 with respect to (3rd part of equation 42 of Panziani et al, PRB 73, 155111 (2006))
 END_DOC
 implicit none
 double precision rs,mu,xi,d_xi_delta_4,denominator_delta,d_xi_denominator_delta,delta_4
 d_xi_3rd_deltaterm=d_xi_delta_4(rs,xi)*mu**4/denominator_delta(rs,xi,mu)-d_xi_denominator_delta(rs,xi,mu)*delta_4(rs,xi)*mu**4/denominator_delta(rs,xi,mu)**2
 return
 end 

!!!!!!!!fourth term!!!!!

 double precision function delta_5(rs,xi)
 BEGIN_DOC
  ! eq 45 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,d0delta,C3t_delta
 delta_5=C3t_delta(rs,xi)*d0delta(rs,xi)**8
 return
 end

 double precision function d_delta_5(rs,xi)
 BEGIN_DOC
  ! Derivative of delta_5 with respect to rs 
 END_DOC
 implicit none
 double precision rs,xi,d_C3t_delta,d0delta,C3t_delta,d_d0delta
 d_delta_5=d_C3t_delta(rs,xi)*d0delta(rs,xi)**8+8d0*C3t_delta(rs,xi)*d_d0delta(xi)*d0delta(rs,xi)**7
 return
 end

 double precision function d_xi_delta_5(rs,xi)
 BEGIN_DOC
  ! Derivative of delta_5 with respect to xi
 END_DOC
 implicit none
 double precision rs,xi,d_xi_C3t_delta,d0delta,C3t_delta,d_xi_d0delta
 d_xi_delta_5=d_xi_C3t_delta(rs,xi)*d0delta(rs,xi)**8+8d0*C3t_delta(rs,xi)*d_xi_d0delta(rs,xi)*d0delta(rs,xi)**7
 return
 end

 double precision function d_4th_deltaterm(rs,mu,xi)
 BEGIN_DOC
  ! Derivative of delta_5*mu**5/(1+d0**2*mu**2)**4 with respect to rs (4th part of equation 42 of Panziani et al, PRB 73, 155111 (2006))
 END_DOC
 implicit none
 double precision rs,mu,xi,d_delta_5,denominator_delta,d_denominator_delta,delta_5
 d_4th_deltaterm=d_delta_5(rs,xi)*mu**5/denominator_delta(rs,xi,mu)-d_denominator_delta(rs,xi,mu)*delta_5(rs,xi)*mu**5/denominator_delta(rs,xi,mu)**2
 return
 end

 double precision function d_xi_4th_deltaterm(rs,mu,xi)
 BEGIN_DOC
  ! Derivative of delta_5*mu**5/(1+d0**2*mu**2)**4 with respect to xi (4th part of equation 42 of Panziani et al, PRB 73, 155111 (2006))
 END_DOC
 implicit none
 double precision rs,mu,xi,d_xi_delta_5,denominator_delta,d_xi_denominator_delta,delta_5
 d_xi_4th_deltaterm=d_xi_delta_5(rs,xi)*mu**5/denominator_delta(rs,xi,mu)-d_xi_denominator_delta(rs,xi,mu)*delta_5(rs,xi)*mu**5/denominator_delta(rs,xi,mu)**2
 return
 end


!!!!!!!!fifth term!!!!!

 double precision function delta_6(rs,xi)
 BEGIN_DOC
  ! eq 45 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,d0delta,C2_delta
 delta_6=C2_delta(rs,xi)*d0delta(rs,xi)**8
 return    
 end 
      
 double precision function d_delta_6(rs,xi)
 BEGIN_DOC
  ! Derivative of delta_6 with respect to rs 
 END_DOC
 implicit none
 double precision rs,xi,d_C2_delta,d0delta,C2_delta,d_d0delta
 d_delta_6=d_C2_delta(rs,xi)*d0delta(rs,xi)**8+8d0*C2_delta(rs,xi)*d_d0delta(xi)*d0delta(rs,xi)**7
 return
 end

 double precision function d_xi_delta_6(rs,xi)
 BEGIN_DOC
  ! Derivative of delta_6 with respect to xi 
 END_DOC
 implicit none
 double precision rs,xi,d_xi_C2_delta,d0delta,C2_delta,d_xi_d0delta
 d_xi_delta_6=d_xi_C2_delta(rs,xi)*d0delta(rs,xi)**8+8d0*C2_delta(rs,xi)*d_xi_d0delta(rs,xi)*d0delta(rs,xi)**7
 return
 end
    
 double precision function d_5th_deltaterm(rs,mu,xi)
 BEGIN_DOC
  ! Derivative of delta_6*mu**6/(1+d0**2*mu**2)**4 with respect to rs(5th part of equation 42 of Panziani et al, PRB 73, 155111 (2006))
 END_DOC
 implicit none
 double precision rs,mu,xi,d_delta_6,denominator_delta,d_denominator_delta,delta_6
 d_5th_deltaterm=d_delta_6(rs,xi)*mu**6/denominator_delta(rs,xi,mu)-d_denominator_delta(rs,xi,mu)*delta_6(rs,xi)*mu**6/denominator_delta(rs,xi,mu)**2
 return
 end

 double precision function d_xi_5th_deltaterm(rs,mu,xi)
 BEGIN_DOC
  ! Derivative of delta_6*mu**6/(1+d0**2*mu**2)**4 with respect to xi (5th part of equation 42 of Panziani et al, PRB 73, 155111 (2006))
 END_DOC
 implicit none
 double precision rs,mu,xi,d_xi_delta_6,denominator_delta,d_xi_denominator_delta,delta_6
 d_xi_5th_deltaterm=d_xi_delta_6(rs,xi)*mu**6/denominator_delta(rs,xi,mu)-d_xi_denominator_delta(rs,xi,mu)*delta_6(rs,xi)*mu**6/denominator_delta(rs,xi,mu)**2
 return
 end
  

!!!!!!!!DELTA TOTAL ET DERIVATIVE!!!!!!!



 double precision function delta_barth(rs,xi,mu)
 BEGIN_DOC
  ! eq 42 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,mu,delta_2,delta_3,delta_4,delta_5,delta_6,denominator_delta
 delta_barth=(delta_2(rs)*mu**2+delta_3(rs,xi)*mu**3+delta_4(rs,xi)*mu**4+delta_5(rs,xi)*mu**5+delta_6(rs,xi)*mu**6)/denominator_delta(rs,xi,mu)
 return    
 end 

 double precision function d_delta_barth(rs,xi,mu)
 BEGIN_DOC
  ! Derivative with respect of rs of eq 42 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none
 double precision rs,xi,mu,d_1st_deltaterm,d_2nd_deltaterm,d_3rd_deltaterm,d_4th_deltaterm,d_5th_deltaterm
 d_delta_barth=d_1st_deltaterm(rs,mu,xi)+d_2nd_deltaterm(rs,mu,xi)+d_3rd_deltaterm(rs,mu,xi)+d_4th_deltaterm(rs,mu,xi)+d_5th_deltaterm(rs,mu,xi)
 return
 end

 double precision function d_xi_delta_barth(rs,xi,mu)
 BEGIN_DOC
  ! Derivative with respect of xi of eq 42 of Panziani et al, PRB 73, 155111 (2006) 
 END_DOC
 implicit none        
 double precision rs,xi,mu,d_xi_1st_deltaterm,d_xi_2nd_deltaterm,d_xi_3rd_deltaterm,d_xi_4th_deltaterm,d_xi_5th_deltaterm
 d_xi_delta_barth=d_xi_1st_deltaterm(rs,mu,xi)+d_xi_2nd_deltaterm(rs,mu,xi)+d_xi_3rd_deltaterm(rs,mu,xi)+d_xi_4th_deltaterm(rs,mu,xi)+d_xi_5th_deltaterm(rs,mu,xi)
 return
 end 

 double precision function wignerseitz_radius(rho)
 BEGIN_DOC
  ! Wigner Seitz radius 
 END_DOC
 implicit none
 double precision pi,rho
 pi=dacos(-1.d0) 
 wignerseitz_radius = (4d0*pi*rho/3d0)**(-1d0/3d0)
 return
 end  

 double precision function d_wignerseitz_radius(rho)
 BEGIN_DOC
  ! derivative of Wigner Seitz radius 
 END_DOC 
 implicit none
 double precision pi,rho
 pi=dacos(-1.d0) 
 d_wignerseitz_radius = -(6d0**(2d0/3d0)*rho**(4d0/3d0)*pi**(1d0/3d0))**(-1) 
 return
 end 

 double precision function d_xi_rhoa(rhoa,rhob)
 BEGIN_DOC
  ! derivative of xi with respect to rhoa 
 END_DOC
 implicit none
 double precision rhoa,rhob
 d_xi_rhoa = 2d0*rhob/(rhoa+rhob)**2 
 return
 end

 double precision function d_xi_rhob(rhoa,rhob)
 BEGIN_DOC
  ! derivative of xi with respect to rhob 
 END_DOC
 implicit none
 double precision rhoa,rhob
 d_xi_rhob = 2d0*rhoa/(rhoa+rhob)**2
 return
 end

 double precision function d_total_deltarho_rhoa(rhoa,rhob,mu)
 BEGIN_DOC
  ! derivative of Delta_LR-SR*rho(r) with respect to rho_alpha
 END_DOC
 implicit none
 double precision d_wignerseitz_radius,wignerseitz_radius,d_xi_rhoa,xi,mu,d_delta_barth,d_xi_delta_barth,delta_barth,rhot,rs,drs,dxi_dna,rhoa,rhob
 rhot=rhoa+rhob
 rs = wignerseitz_radius(rhot)
 drs= d_wignerseitz_radius(rhot)
 xi= (rhoa-rhob)/(rhoa+rhob)
 dxi_dna=d_xi_rhoa(rhoa,rhob)
 d_total_deltarho_rhoa=(d_delta_barth(rs,xi,mu)*drs+d_xi_delta_barth(rs,xi,mu)*dxi_dna)*rhot+delta_barth(rs,xi,mu)  
 return
 end

 double precision function d_total_deltarho_rhob(rhoa,rhob,mu)
 BEGIN_DOC
  ! derivative of Delta_LR-SR*rho(r) with respect to rho_beta
 END_DOC
 implicit none
 double precision d_wignerseitz_radius,wignerseitz_radius,d_xi_rhob,xi,mu,d_delta_barth,d_xi_delta_barth,delta_barth,rhot,rs,drs,dxi_dnb,rhoa,rhob
 rhot=rhoa+rhob
 rs = wignerseitz_radius(rhot)
 drs= d_wignerseitz_radius(rhot)
 xi= (rhoa-rhob)/(rhoa+rhob)
 dxi_dnb=d_xi_rhob(rhoa,rhob)
 d_total_deltarho_rhob=(d_delta_barth(rs,xi,mu)*drs+d_xi_delta_barth(rs,xi,mu)*dxi_dnb)*rhot+delta_barth(rs,xi,mu)
 return
 end

!C****************************************************************************
      subroutine ESRC_MD_LDAERF_barth (mu,rho_a,rho_b,dospin,e)
!C*****************************************************************************
!C     Short-range spin-dependent LDA correlation functional with multideterminant reference
!C       for OEP calculations from Section V of 
!C       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
!C
!C     Input: rhot   : total density
!C            rhos   : spin density
!!C            mu     : Interation parameter
!C            dospin : use spin density
!C
!C     Ouput: e      : energy
!C
!     Created: 26-08-11, J. Toulouse
!C*****************************************************************************
      implicit none

      double precision, intent(in) :: rho_a,rho_b,mu
      logical, intent(in)          :: dospin
      double precision, intent(out):: e
 
      double precision             :: e1
      double precision             :: rhoa,rhob
      double precision             :: rhot, rhos
      double precision             :: xi,rs,pi,delta_barth
      pi=dacos(-1.d0)
      rhoa=max(rho_a,1.0d-15)
      rhob=max(rho_b,1.0d-15)
      rhot = rhoa + rhob
      rhos = rhoa - rhob
    
      call ec_only_lda_sr(mu,rho_a,rho_b,e1)
 
     !call DELTA_LRSR_LDAERF (rhot,rhos,mu,dospin,e)

      !rhoa=max((rhot+rhos)*.5d0,1.0d-15)
      !rhob=max((rhot-rhos)*.5d0,1.0d-15)
      xi=(rhoa-rhob)/(rhoa+rhob)      
      
      rs = (4d0*pi*rhot/3d0)**(-1d0/3d0)

      e = delta_barth(rs,xi,mu)*rhot


      e = e1 + e
     
      end 


