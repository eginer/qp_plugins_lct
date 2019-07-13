      subroutine DsplIntegr (X,

     .             Y,

     .             nbPoints,

     .             xlow,

     .             xup,

     .             Integral,

     .             error)

      implicit none

      integer error, Nbpoints

      double precision X(NbPoints), Y(NbPoints), Xlow, Xup

      double precision Integral



* local

      double precision Ypp(Nbpoints+4), delta3, delta

      double precision Yp1, YpN, DDERIVLAG4

      double precision Y1, YN, DINTERLAG4

      integer NbInterv

      integer i, k, point, Icur, Ifirst, Ilast

      integer Mxpoints, NbIntPnt

      double precision YL(Nbpoints+4), XL(Nbpoints+4), A,B,C,D,XX

      double precision YlL(Nbpoints+4), XLl(Nbpoints+4)

      double precision threshrange

      parameter(threshrange=1.d-15)



      logical debug

      debug=.true.

      debug=.false.

* interpolate Y by Cubic-Splines



      if(xup.lt.xlow) then

       print *,'DsplIntegr: integration from r=',xlow,' until r=',xup

       print *,' DsplIntegr: nbPoints =',nbPoints

       print *,' DsplIntegr: Xup =',Xup,' Xlow=',xlow

       print *,' DsplIntegr: ERROR Xup < Xlow '

       stop ' DsplIntegr: ERROR Xup < Xlow '

      end if

     

      if(NbPoints.lt.4) then

       print *,' DsplIntegr: ERROR Nbpoints =',NbPoints,' < 4 '

       stop ' DsplIntegr: ERROR nb of points < 4 '

      end if



      call DLOCATE (X, NbPoints, Xlow, Ifirst)

      if(debug) print *,' DsplIntegr: beginning Ifirst ',Ifirst

      Icur=Max(2,Ifirst-1) 

      Icur=Min(NbPoints-2,ICur)



      if(debug) then

       do i=1,Nbpoints

      print *,'  '

      print *,' DsplIntegr: X=',X(i),' Y=',Y(i)

       end do

      print *,' DsplIntegr: beginning Icur ',Icur

      print *,' DsplIntegr: beginning Xlow ',Xlow

      print *,' DsplIntegr: beginning interpolation from ',Icur-1,

     .        ' to ',Icur+2

      end if



        y1=DINTERLAG4 (y(Icur-1), y(Icur), y(Icur+1), y(Icur+2),

     .                 x(Icur-1), x(Icur), x(Icur+1), x(Icur+2),

     .                 xlow)



        yp1=DDERIVLAG4 (y(Icur-1), y(Icur), y(Icur+1), y(Icur+2),

     .                 x(Icur-1), x(Icur), x(Icur+1), x(Icur+2),

     .                 xlow)

      if(debug) then

      print *,' DsplIntegr: Dlocate ',xlow,' at ',Ifirst,' F=',Y1

      end if



* Last point:

      call DLOCATE (X, NbPoints, Xup, Ilast)

      Icur=Min(NbPoints-2,Ilast)

      Icur=Max(2,Icur) 



      if(debug) then

      print *,' DsplIntegr: end Ilast ',Ilast

      print *,' DsplIntegr: end Icur ',Icur

      print *,' DsplIntegr: end Xup  ',Xup 

      print *,' DsplIntegr: end interpolation from ',Icur-1,

     .        ' to ',Icur+2

      end if



      if(Ilast.eq.0) then

         Ilast=Ilast+1

      else if(Xup .gt. x(Ilast)) then

         Ilast=Ilast+1

      end if



        yN=DINTERLAG4 (y(Icur-1), y(Icur), y(Icur+1), y(Icur+2),

     .                 x(Icur-1), x(Icur), x(Icur+1), x(Icur+2),

     .                 xup)



        ypN=DDERIVLAG4 (y(Icur-1), y(Icur), y(Icur+1), y(Icur+2),

     .                 x(Icur-1), x(Icur), x(Icur+1), x(Icur+2),

     .                 xup)

      nbIntPnt=Ilast-Ifirst+1

      if(debug) then

      print *,' DsplIntegr: Dlocate ',Xup,' at ',Ilast,' F=',YN      

      print *,' derivative at 1=',yp1

      print *,' derivative at N=',ypN

      print *,' Nb of Integration points=',nbIntPnt

      end if



      k=1

      do i=Ifirst+1,Ilast-1

         k=k+1

         YL(k)=Y(i)

         XL(k)=X(i)

      end do

      NbInterv=nbIntPnt-1



      XL(1)=Xlow

      XL(nbIntPnt)=Xup



      YL(1)=Y1

      YL(nbIntPnt)=YN



*     if(debug) then

*     call DWRITEXGR ('y.spl', 'spline', 'yL', 'X', 'Y',

*    .                     1, nbIntPnt, XL, YL)

*     end if



      call DSPLINE(XL,

     .            YL,

     .            nbIntPnt,

     .            yp1, ypN,

     .            Ypp)



*     if(debug) then

*     call DWRITEXGR ('y2.spl', 'spline', 'y''''', 'X', 'Y',

*    .                     1, nbIntPnt, XL, Ypp)

*     end if



* A= (x(i+1)-x)/delta(i)

* B= (x-x(i))/delta(i)

* C= (A^3-A)*delta(i)^2/6

* D= (B^3-B)*delta(i)^2/6

* Y(x)= A*y(i) + B*y(i+1) + C*ypp(i) + D*ypp(i+1)



* dA = -dx/delta

* dB =  dx/delta

* Int of A on [0.d0.1d0]= delta/2

* Int of B on [0.d0.1d0]= delta/2

* Int of A^3 on [0.d0.1d0]= delta/4

* Int of B^3 on [0.d0.1d0]= delta/4

* Int of C on [0.d0.1d0]= -delta(i)^3/24

* Int of D on [0.d0.1d0]= -delta(i)^3/24



* Int of Y on [x(i)..x(i+1)]=

*  (y(i) + y(i+1))/2 + (delta(i)^2)(ypp(i) + ypp(i+1))/24



* integrate Cubic-Splines on each of the nbIntPnt-1 intervals

      integral=0.d0

      do i=1, NbInterv

        delta=xL(i+1)-xL(i)

        delta3=delta*delta*delta

        integral=integral+

     . (delta*(yL(i) + yL(i+1))/2.d0) - 

     . (delta3*(ypp(i) + ypp(i+1))/24.d0)

      end do



      if(debug) then

      print *,' integral=',integral



      do i=1, NbInterv

        delta=xL(i+1)-xL(i)

        XX= xL(i)+(delta/2.d0)

        XLL(i)=XX

        A= (xL(i+1)-xX)/delta

        B= (xX-xL(i))/delta

        C= (A**3-A)*(delta**2)/6.d0

        D= (B**3-B)*(delta**2)/6.d0

        YLL(i)= A*yL(i) + B*yL(i+1) + C*ypp(i) + D*ypp(i+1)

      end do



* YLL : spline interpolated function at mid-points:

*     call DWRITEXGR ('f.spl', 'spline', 'int', 'X', 'Y',

*    .                     1, nbinterv, XLL, YLL)



      end if ! debug



      return

      end



      SUBROUTINE dspline(x,y,n,yp1,ypn,y2)

      INTEGER n,NMAX

      double precision yp1,ypn,x(n),y(n),y2(n)

      PARAMETER (NMAX=500)

      INTEGER i,k

      double precision p,qn,sig,un,u(N)

      

      if(n.eq.0) then

        write(6,'(a)') ' error_dspline> n=0'

        stop ' error_dspline> n=0'

      end if

      if (yp1.gt..99d30) then

        y2(1)=0.d0

        u(1)=0.d0

      else

        y2(1)=-0.5d0

        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)

      endif

      do 11 i=2,n-1

        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))

        p=sig*y2(i-1)+2.d0

        y2(i)=(sig-1.d0)/p

        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+

     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*

     *u(i-1))/p

11    continue

      if (ypn.gt..99d30) then

        qn=0.d0

        un=0.d0

      else

        qn=0.5d0

        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))

      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

      do 12 k=n-1,1,-1

        y2(k)=y2(k)*y2(k+1)+u(k)

12    continue

      return

      END


* =======================================================================
* *

      double precision function DDERIVLAG4 (f1, f2, f3, f4, 

     .                                      x1, x2, x3, x4, 

     .                                      x)

* -----------------------------------------------------------------------
* *

      implicit none

      double precision f1, f2, f3, f4, f

      double precision x1, x2, x3, x4, x



       F = f1*(x-x3)*(x-x4)/((x1-x2)*(x1-x3)*(x1-x4))   ! x1

     .   + f1*(x-x2)*(x-x4)/((x1-x2)*(x1-x3)*(x1-x4))   ! x1

     .   + f1*(x-x2)*(x-x3)/((x1-x2)*(x1-x3)*(x1-x4))   ! x1

     .   + f2*(x-x3)*(x-x4)/((x2-x1)*(x2-x3)*(x2-x4))   ! x2

     .   + f2*(x-x1)*(x-x4)/((x2-x1)*(x2-x3)*(x2-x4))   ! x2

     .   + f2*(x-x1)*(x-x3)/((x2-x1)*(x2-x3)*(x2-x4))   ! x2

     .   + f3*(x-x2)*(x-x4)/((x3-x1)*(x3-x2)*(x3-x4))   ! x3

     .   + f3*(x-x1)*(x-x4)/((x3-x1)*(x3-x2)*(x3-x4))   ! x3

     .   + f3*(x-x1)*(x-x2)/((x3-x1)*(x3-x2)*(x3-x4))   ! x3

     .   + f4*(x-x2)*(x-x3)/((x4-x1)*(x4-x2)*(x4-x3))   ! x4

     .   + f4*(x-x1)*(x-x3)/((x4-x1)*(x4-x2)*(x4-x3))   ! x4

     .   + f4*(x-x1)*(x-x2)/((x4-x1)*(x4-x2)*(x4-x3))   ! x4



*       write(*,*) ' x=',x,' x1=',x1,' x2=',x2,' x3=',

*    .                       x3,' x4=',x4



*       write(*,*) ' f=',f,' f1=',f1,' f2=',f2,' f3=',

*    .                       f3,' f4=',f4



*      write (*,*) ' at X=',X,' F''=',F



      DDERIVLAG4=F

      return

      end



* =======================================================================
* *

      double precision function DINTERLAG4 (f1, f2, f3, f4, 

     .                                      x1, x2, x3, x4, 

     .                                      x)

* -----------------------------------------------------------------------
* *

      implicit none

      double precision f1, f2, f3, f4, f

      double precision x1, x2, x3, x4, x

      logical debug



      debug=.true.

      debug=.false.



       F = f1*(x-x2)*(x-x3)*(x-x4)/((x1-x2)*(x1-x3)*(x1-x4))   ! x1

     .   + f2*(x-x1)*(x-x3)*(x-x4)/((x2-x1)*(x2-x3)*(x2-x4))   ! x2

     .   + f3*(x-x1)*(x-x2)*(x-x4)/((x3-x1)*(x3-x2)*(x3-x4))   ! x3

     .   + f4*(x-x1)*(x-x2)*(x-x3)/((x4-x1)*(x4-x2)*(x4-x3))   ! x4



      DINTERLAG4=F

      if(debug) then

        print *,' DINTERLAG4: X1=',X1,' X2=',X2,' X3=',X3,' X4=',X4

        print *,' DINTERLAG4: F1=',F1,' F2=',F2,' F3=',F3,' F4=',F4

        print *,' DINTERLAG4: X=',X,' F=',F

      end if

      return

      end


      subroutine DLOCATE(xx,n,x,j)

      INTEGER j,n

      double precision x,xx(n)

      INTEGER jl,jm,ju

      jl=0

      ju=n+1

10    if(ju-jl.gt.1)then

        jm=(ju+jl)/2

        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then

          jl=jm

        else

          ju=jm

        endif

      goto 10

      endif

      if(x.eq.xx(1))then

        j=1

      else if(x.eq.xx(n))then

        j=n-1

      else

        j=jl

      endif

      return

      END
