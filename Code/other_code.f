cccccccccccccccccccccccccccccccccccccccccccc
      subroutine caleulang(mx,vv,ising)     
cccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer irdummy,ising,ising_1
      real(8) det
      real(8) vv(3)
      real(8) mx(3,3),rm(3,3),um(3,3)
      real(8) pi,r2g,g2r
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180
cccccccccccccccccccccccccccccccccccccccccccc
      ising_1=0
!      call pdecomposition(mx,um,rm,ising_1)
      call polar_decomp(mx,um,rm,ising_1)
      if(ising_1/=0)then
         write(6,*) 'pdecompsition of mx is failure'
         call pm(mx,3,3)
         ising=1
         return
      endif
      call icams_q2eang(transpose(rm),vv(1),vv(2),vv(3))
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_misori(v1,v2,ang)
ccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) v1(3),v2(3),ang,x1,x2
      real(8) qm1(3,3),qm2(3,3),dqm(3,3)
      real(8) pi,r2g,g2r
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180
      call icams_eang2q(v1(1),v1(2),v1(3),qm1)
      call icams_eang2q(v2(1),v2(2),v2(3),qm2)
      dqm=matmul(qm2,transpose(qm1))
      x1=dqm(1,1)+dqm(2,2)+dqm(3,3)
      x2=(x1-1.d0)/2
      if(dabs(x2)>1.d0)x2=1.d0*sign(1.d0,x2)
      ang=dabs(pi/2-dasin(x2)) *r2g
      return
      end

c****************************************************************
      subroutine pdecomposition(mx,umx,rmx,ising)
c****************************************************************
      implicit none
      integer ising
      real(8) mx(3,3),ce(3,3),rmx(3,3),umx(3,3),iumx(3,3)
      real(8) eb1(3,3),eb2(3,3),eb3(3,3)
      real(8) ev1,ev2,ev3,det
      ising=0
      ce=matmul(transpose(mx),mx)
      call spectral(ce,ev1,ev2,ev3,eb1,eb2,eb3,ising)
      if(ev1<=0 .or. ev2<=0 .or. ev3<=0 .or. ising/=0)then
         write(6,*) 'eigen value of ce <0'
         print*, ev1,ev2,ev3
         ising=1
         return
      endif
      umx=dsqrt(ev1)*eb1+dsqrt(ev2)*eb2+dsqrt(ev3)*eb3
      iumx=1/dsqrt(ev1)*eb1+1/dsqrt(ev2)*eb2+1/dsqrt(ev3)*eb3
      rmx=matmul(mx,iumx)
      return 
      end

      subroutine icams_determ(a,det)
c********************************************************************
c     this routine calculates the determinant of a
c********************************************************************
      implicit none
      real(8) a(3,3),v1,v2,v3,det
      v1= a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      v2= a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      v3= a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      det= v1-v2+v3
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_eang2q(p1,p,p2,qm)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     rotate from x[100],y[010],z[001] to v1,v2,v3
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) qm(3,3)
      real(8) p1,p,p2,xp1,xp,xp2
      real(8) c1,c,c2,s1,s,s2
      real(8) pi,r2g,g2r
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180
      xp1=p1*g2r
      xp =p *g2r
      xp2=p2*g2r
      c1=dcos(xp1)
      s1=dsin(xp1)
      c =dcos(xp )
      s =dsin(xp )
      s2=dsin(xp2)
      c2=dcos(xp2)
      qm(1,1)=+c1*c2-s1*s2*c
      qm(1,2)=+s1*c2+c1*s2*c
      qm(1,3)=+s2*s
      qm(2,1)=-c1*s2-s1*c2*c
      qm(2,2)=-s1*s2+c1*c2*c
      qm(2,3)=+c2*s
      qm(3,1)=+s1*s
      qm(3,2)=-c1*s
      qm(3,3)=+c
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_q2eang(qm,phi1,phi,phi2)
cccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) qm(3,3)
      real(8) phi1,phi,phi2
      real(8) sqhkl,squvw,sqhk,val
      real(8) pi,r2g,g2r,tol
      tol=1.d-15
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180
c---------------------------------------------
c
c             v1   v2    v3
c
c           | u   v2_1    h  |
c     qm =  | v   v2_2    k  |  with v2 = v3 x v1   
c           | w   v2_3    l  | 
c
c             100   010   001
c          x|  u   v2_1    h  |
c     qm = y|  v   v2_2    k  |    
c          z|  w   v2_3    l  | 
c
c---------------------------------------------


c---------------------------------------------
c
c    if the roation tensor is defined as following
c    qm_ij:= (dx/dx)_ij = dx_i/dx_j
c
c               x     y     x
c          v1|  u   v2_1    h  |
c     qm'= v2|  v   v2_2    k  |    
c          v3|  w   v2_3    l  | 
c    the ratation must be take 
c
c    qm=transpose(qm')
c
c---------------------------------------------
      squvw=dsqrt(qm(1,1)**2+qm(2,1)**2+qm(3,1)**2)
      sqhkl=dsqrt(qm(1,3)**2+qm(2,3)**2+qm(3,3)**2)
      sqhk =dsqrt(qm(1,3)**2+qm(2,3)**2           )

      val=qm(3,3)/sqhkl
      if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
      phi=dacos(val)

      if(phi < tol) then
         phi2=0.0
         val=qm(1,1)/squvw
         if(qm(2,1) <= 0.d0) then
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
            phi1=dacos(val)
         else
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
c           phi1=2*pi-dacos(val)
            phi1=-dacos(val)
         endif
      else
c        val=qm(2,3)/sqhk
         val=qm(2,3)/dsin(phi)
         if(qm(1,3) >= 0.d0) then
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
            phi2=dacos(val)
         else
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
c           phi2=2*pi-dacos(val)
            phi2=-dacos(val)
         endif
         val=-qm(3,2)/dsin(phi)
         if(qm(3,1) >= 0.d0) then
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
            phi1=dacos(val)
         else
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
c           phi1=2*pi-dacos(val)
            phi1=-dacos(val)
         endif
      endif
      phi1=phi1*r2g
      phi=phi*r2g
      phi2=phi2*r2g
      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_q2eang_old(qm,phi1,phi,phi2)
cccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) qm(3,3)
      real(8) phi1,phi,phi2,x1,x2,x3
      real(8) pi,r2g,g2r,tol
      tol=1.d-15
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180
c---------------------------------------------
c            v1   v2=v3xv1    v3
c
c           | u   v2_1    h  |
c     qm =  | v   v2_2    k  |    
c           | w   v2_3    l  | 
c
c---------------------------------------------

      x1=dsqrt(qm(1,1)**2+qm(2,1)**2+qm(3,1)**2)
      x2=dsqrt(qm(1,2)**2+qm(2,2)**2+qm(3,2)**2)
      x3=dsqrt(qm(1,3)**2+qm(2,3)**2+qm(3,3)**2)
      qm(:,1)=qm(:,1)/x1
      qm(:,2)=qm(:,2)/x1
      qm(:,3)=qm(:,3)/x1

      x1=qm(3,3)
      if(dabs(x1)>1.d0)x2=1.d0*sign(1.d0,x1)

      phi=dacos(x2)

      if(phi<tol)then
         phi2=0.d0

         x1=qm(1,1)
         if(dabs(x1)>1.d0)x2=1.d0*sign(1.d0,x1)

         if(qm(1,2)>=0.d0)then
            phi1=dacos(x2)
         else
c           phi1=2*pi-dacos(val)
            phi1=-dacos(x2)
         endif
      else

         x1=qm(2,3)/dsin(phi)
         if(dabs(x1)>1.d0)x2=1.d0*sign(1.d0,x1)

         if(qm(1,3)>=0.d0) then
            phi2=dacos(x2)
         else
c           phi2=2*pi-dacos(val)
            phi2=-dacos(x2)
         endif

         x1=-qm(3,2)/dsin(phi)
         if(dabs(x1)>1.d0)x2=1.d0*sign(1.d0,x1)

         if(qm(3,1) >= 0.d0) then
            phi1=dacos(x2)
         else
c           phi1=2*pi-dacos(val)
            phi1=-dacos(x2)
         endif
      endif
      phi1=phi1*r2g
      phi=phi*r2g
      phi2=phi2*r2g
      return
      end




cccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_angax2qm(ang,u,v,w,qm)
cccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) qm(3,3)
      real(8) s,c,u2,v2,w2,ang,u,v,w,x1
      real(8) pi,r2g,g2r
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180

      x1=dsqrt(u**2+v**2+w**2)
      u2=u/x1
      v2=v/x1
      w2=w/x1
      s=dsin(ang)
      c=dcos(ang)

      qm(1,1)=(1-u2**2)*c+u2**2
      qm(2,2)=(1-v2**2)*c+v2**2
      qm(3,3)=(1-w2**2)*c+w2**2

      qm(1,2)=u2*v2*(1-c)+w2*s
      qm(2,1)=u2*v2*(1-c)-w2*s

      qm(1,3)=u2*w2*(1-c)-v2*s
      qm(3,1)=u2*w2*(1-c)+v2*s

      qm(2,3)=v2*w2*(1-c)+u2*s
      qm(3,2)=v2*w2*(1-c)-u2*s


      return
      end

c
c
c********************************************************************** 
      subroutine hi(m,hi1m,hi2m,hi3m)
c**** hauptinvarianten hi1m, hi2m, hi3m der 3x3 matrix m
      implicit none
      real(8) m(3,3),hi1m,hi2m,hi3m 
      hi1m = m(1,1)+m(2,2)+m(3,3)
      hi2m =(m(1,1)+m(2,2)+m(3,3))**2/2.d0
     &      -m(1,1)**2/2.d0
     &      -m(2,2)**2/2.d0
     &      -m(3,3)**2/2.d0
     &      -m(1,2)*m(2,1)
     &      -m(1,3)*m(3,1)
     &      -m(2,3)*m(3,2)
      hi3m =+m(1,1)*m(2,2)*m(3,3)
     &      +m(2,1)*m(1,3)*m(3,2)
     &      +m(3,1)*m(1,2)*m(2,3)
     &      -m(1,1)*m(2,3)*m(3,2)
     &      -m(2,1)*m(1,2)*m(3,3)
     &      -m(3,1)*m(1,3)*m(2,2)
      return  
      end


c**********************************************************************
      subroutine norm(m,zm,sm,no)
c**********************************************************************
      implicit none
      integer i,j,zm,sm
      double precision m(zm,sm),no  
      no=0.d0
      do i=1,zm
      do j=1,sm 
         no=no+m(i,j)**2.d0 
      enddo 
      enddo
      no=dsqrt(no)
      return 
      end  


c**********************************************************************
      subroutine wm(m,ni,nj)
      implicit none
      integer ni,nj,i,j
      real(8) m(ni,nj)
      write(6,*) 
      do i=1,ni
         write(6,100) (m(i,j),j=1,nj)          
      enddo
      write(6,*) 
100   format(48e20.8)       
      return
      end


c**********************************************************************
      subroutine pm(m,ni,nj)
      implicit none
      integer ni,nj,i,j
      real(8) m(ni,nj)
      write(*,*) 
      do i=1,ni
         write(*,100) (m(i,j),j=1,nj)          
      enddo
      write(*,*) 
c100   format(48e20.8)        
100   format(48e12.4)       
c100   format(48f6.3)       
      return
      end

c**********************************************************************
      subroutine pv(m,ni)
      implicit none
      integer ni,i
      real(8) m(ni)
      write(*,100) (m(i),i=1,ni)          
c      write(6,*) 
c100   format(48e20.8)       
100   format(48e12.4)       
      return
      end


c******************************************** 
      subroutine gaussj(a,n,b,ising)
c******************************************** 
      implicit none
      integer   n,ising
      integer   i,icol,irow,j,k,l,ll,i1,i2
      integer   indxc(n),indxr(n),ipiv(n)
      real(8)   a(n,n),b(n,n),c(n,n),vx(n)
      real(8)   big,dum,pivinv,x1,x2,x3
c-------------------------------------------
c      write(6,*) 'coming to jd'
c      call pm(a,3,3)
c      call flush(6)

      ising=0
      c=a
      ipiv=0
c-----loop for the pivot procedures, from 1 to n
      do i=1,n
         big=0.d0
         do j=1,n
         do k=1,n
            if(ipiv(j).ne.1 .and. ipiv(k)==0)then
               if(dabs(a(j,k)).ge.big)then
                  big=dabs(a(j,k))
                  irow=j
                  icol=k
               endif
            endif
            if(ipiv(j).ne.1 .and. ipiv(k).gt.1)then
c              print*,'sigular matrix in gauss_jordan'
c              write(6,*)'sigular matrix in gauss_jordan'
c              call flush(6)
               ising=1
               return
            endif
         enddo
         enddo
c--------check whether the pivot element is zero or not
         if(a(irow,icol)==0.d0)then
c	      print*,'sigular matrix in gauss_jordan'
c	      print*,'indices is:', irow,icol
c           write(6,*)'sigular matrix in gauss_jordan'
c           write(6,*)'indices is:', irow,icol
c           call flush(6)
            ising=1
            return
         endif
c-------------------------------------------------------
c        if one component is selected as pivot element
c        the second indice is important, so it is marked
c        from 0 to 1 in ipiv(:) array
c        after convert, it is the row number
c-------------------------------------------------------
         ipiv(icol)=ipiv(icol)+1
c--------record the row and collum number for ith pivot element
         indxr(i)=irow
         indxc(i)=icol
c--------change pivot element to diagonal position
         if(irow.ne.icol)then
            vx=a(irow,:)
            a(irow,:)=a(icol,:)
            a(icol,:)=vx
         endif
c--------eliminate the elements besides a(icol,icol)
         pivinv=1.d0/a(icol,icol)
         a(icol,icol)=1.d0
         a(icol,:)=a(icol,:)*pivinv 
         do i2=1,n
            if(i2.ne.icol)then
               dum=a(i2,icol)
               a(i2,icol)=0.d0
               a(i2,:)=a(i2,:)-a(icol,:)*dum
            endif
         enddo
      enddo

c-----after maximum pivot strategy elimination
c-----rearrage the left matrix
      do l=n,1,-1
         if(indxr(l).ne.indxc(l))then
            vx=a(:,indxr(l))
            a(:,indxr(l))=a(:,indxc(l))
            a(:,indxc(l))=vx
         endif
      enddo
      b=a
      a=c

      return
      end


c***************************************************************c
c      subroutine spolar (f,u,r)				        c
c							  	        c
c *performs the polar decomposition of the	               c
c   tensor f=ru using a cayley-hamilton theorem.		        c
c *then, u=r^t f is found					        c
c***************************************************************c
      subroutine polar_decomp(f,u,r,ising)
      implicit none
      integer i,j,k,iflag1,iflag2,ising
      real(8) f(3,3),u(3,3),r(3,3),rt(3,3),c(3,3),cs(3,3),uinv(3,3)
      real(8) x1,x2,x3,f1,f2,c1,c2,c3,u1,u2,u3,b,b1,b2,b3,b4
c
      c=matmul(transpose(f),f)
      cs=matmul(c,c)

c-----1st, 2st, 3st invariant for tensor c
      c1=c(1,1)+c(2,2)+c(3,3)
      c2=(c1**2.d0-(cs(1,1)+cs(2,2)+cs(3,3)))/2
      c3=+c(1,1)*(c(2,2)*c(3,3)-c(2,3)*c(3,2))
     1	 -c(1,2)*(c(2,1)*c(3,3)-c(2,3)*c(3,1))
     2	 +c(1,3)*(c(2,1)*c(3,2)-c(2,2)*c(3,1))

c-----3st invariant for tensor u
      u3=dsqrt(c3)

      x1= 2.0**5.0 /27.0 * (2.0*c1**3.0-9.0*c1*c2+27.0*c3)
      x2= 2.**10.0 /27.0 * (4.0*c2**3.0 - c1**2.0*c2**2.0 +
     1	  4.0*c1**3.0*c3 - 18.0 * c1*c2*c3 + 27.0 * c3**2.0)

      if(x2<0)x2=0
      f1=x1+dsqrt(x2)
      iflag1=0
      if(f1<0)iflag1=1
      f2=x1-dsqrt(x2)
      iflag2=0
      if(f2<0)iflag2=1

      if(iflag1==1) f1=-f1
      if(iflag2==1) f2=-f2

      x3= -2.0/3.0*c1 + f1**(1.0/3.0) + f2**(1.0/3.0)
      if(iflag1==1) x3= -2.0/3.0*c1 - f1**(1.0/3.0) + f2**(1.0/3.0)
      if(iflag2==1) x3= -2.0/3.0*c1 + f1**(1.0/3.0) - f2**(1.0/3.0)

c-----1st, 2st invariant for tensor u
      b=-2.0*c1
      if(x3==b)then
         u1=dsqrt(c1+2.0*dsqrt(c2))
      else

         x1=dsqrt(2.0*c1+x3)
         if(x1==0)then
            ising=1
            return
         endif      
         u1= 0.5 * ( x1 + dsqrt(2.0*c1 - x3 + 16.0*dsqrt(c3)/x1) )
      endif
      u2=dsqrt(c2+2.0*u3*u1)

      b1= u3**2.0 * (u3+u1*c1) + u1**2.0 * (u1*c3+u3*c2)
      if(b1==0)then
         ising=1
         return
      endif      

      b2= u1*(u1*u2-u3) / b1
      b3=-(u1*u2-u3) * (u3+u1*c1) / b1
      b4= (u2*u3*(u3+u1*c1) + u1**2.0 * (u2*c2+c3))/b1

      uinv=b2*cs + b3*c
      uinv(1,1)=uinv(1,1)+b4
      uinv(2,2)=uinv(2,2)+b4
      uinv(3,3)=uinv(3,3)+b4

      r=matmul(f,uinv)
      u=matmul(transpose(r),f)

      return
      end


c**********************************************************************
      subroutine spectral(m,ew1,ew2,ew3,eb1,eb2,eb3,ising)
c**********************************************************************
      implicit none
      integer ising
      real(8) m(3,3),eb1(3,3),eb2(3,3),eb3(3,3),ew1,ew2,ew3,
     &        hi1m,hi2m,hi3m,tol,r,s,t,p,q,rho,phi,y1,y2,y3,d1,d2,d3,
     &        e(3,3),m1(3,3),m2(3,3),m3(3,3)
      real(8) x1,x2,x3,pi
      ising=0
      tol=1.d-15
      pi=dacos(-1.d0)
      call hi(m,hi1m,hi2m,hi3m)
      r=-hi1m
      s= hi2m
      t=-hi3m
      p=s-r**2.d0/3.d0
      q=2.d0/27.d0*r**3.d0-r*s/3.d0+t
      rho=dsqrt(-3.d0*p**3.d0)/9.d0
      if(dabs(rho)<=tol)then
         ising=1
         return         
      endif
      !======================================
      !  modification for overcome overflow
      !======================================
      if(dabs(rho)<=tol) rho=sign(1.d0,rho)*tol
      x1=-q/rho/2.d0
      if(dabs(x1)>1.d0) x1=1.d0*sign(1.d0,x1)
      phi=dacos(x1)
      y1=2*rho**(1.d0/3.d0)*dcos(phi/3.d0)
      y2=2*rho**(1.d0/3.d0)*dcos(phi/3.d0+2.d0/3.d0*pi)
      y3=2*rho**(1.d0/3.d0)*dcos(phi/3.d0+4.d0/3.d0*pi)
      ew1=y1-r/3.d0
      ew2=y2-r/3.d0
      ew3=y3-r/3.d0
      e(1,:)=[1,0,0]
      e(2,:)=[0,1,0]
      e(3,:)=[0,0,1]
      eb1=0.d0
      eb2=0.d0
      eb3=0.d0
      if(dabs(ew1-ew2)<tol .and. dabs(ew1-ew3)<tol)then
         eb1=e
      elseif(dabs(ew2-ew3)<tol)then
         eb1=matmul(m-ew2*e,m-ew3*e)/((ew1-ew2)*(ew1-ew3))
         eb2=e-eb1
      elseif(dabs(ew1-ew3)<tol)then
         eb2=matmul(m-ew1*e,m-ew3*e)/((ew2-ew1)*(ew2-ew3))
         eb1=e-eb2
      elseif(dabs(ew1-ew2)<tol)then
         eb3=matmul(m-ew1*e,m-ew2*e)/((ew3-ew1)*(ew3-ew2))
         eb1=e-eb3
      else
         eb1=matmul(m-ew2*e,m-ew3*e)/((ew1-ew2)*(ew1-ew3))
         eb2=matmul(m-ew1*e,m-ew3*e)/((ew2-ew1)*(ew2-ew3))
         eb3=matmul(m-ew1*e,m-ew2*e)/((ew3-ew1)*(ew3-ew2))
      endif

      return
      end


c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   fast fourier transformation code                            +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine fourn(datx,nn,ndim,isg)
      implicit none
      integer isg,ndim,nn(ndim)
      integer i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1
      integer k2,n,nprev,nrem,ntot,ifp2,ip1,ip2,ip3,k1
      real(8) datx(*)
      real(8) tempi,tempr
      real(8) theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then

            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=datx(i3)
                tempi=datx(i3+1)
                datx(i3)=datx(i3rev)
                datx(i3+1)=datx(i3rev+1)
                datx(i3rev)=tempr
                datx(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif

          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isg*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*dsin(0.5d0*theta)**2
          wpi=dsin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*datx(k2)-sngl(wi)*datx(k2+1)
                tempi=sngl(wr)*datx(k2+1)+sngl(wi)*datx(k2)
                datx(k2)=datx(k1)-tempr

                datx(k2+1)=datx(k1+1)-tempi
                datx(k1)=datx(k1)+tempr
                datx(k1+1)=datx(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      end

c**********************************************
      subroutine icams_conv33to6(am,ibx1,ibx2,av)
c**********************************************
      implicit none
      integer i
      integer ibx1(9),ibx2(9)
      real(8) am(3,3),av(6)
      do i=1,6
         av(i)=am(ibx1(i),ibx2(i))
      enddo
      return
      end


c***********************************************
      subroutine icams_conv33to9(am,ibx1,ibx2,av)
c***********************************************
      implicit none
      integer i
      integer ibx1(9),ibx2(9)
      real(8) am(3,3),av(9)
      do i=1,9
         av(i)=am(ibx1(i),ibx2(i))
      enddo
      return
      end


c**************************************************
      subroutine icams_conv6to33(av,ibx1,ibx2,am)
c**************************************************
      implicit none
      integer i
      integer ibx1(9),ibx2(9)
      real(8) am(3,3),av(6)
      do i=1,6
         am(ibx1(i),ibx2(i))=av(i)
      enddo
      do i=7,9
         am(ibx1(i),ibx2(i))=av(i-3)
      enddo
      return
      end


c====================================================================================================
c     computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of size n
c     by n, stored in a physical np by np array. on output, elements of a above the diagonal are
c     destroyed. d returns the eigenvalues of a in its first n elements. v is a matrix with the same
c     logical and physical dimensions as a, whose columns contain, on output, the normalized
c     eigenvectors of a. nrot returns the number of jacobi rotations that were required.
c====================================================================================================
      subroutine jacobi(a,n,np,d,v,nrot,ising)
      implicit none
      integer n,np,nrot,nmax
      real(8) a(np,np),d(np),v(np,np)
      parameter (nmax=500)
      integer i,ip,iq,j,ising
      real(8) c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)
      ising=0
      do ip=1,n  !initialize to the identity matrix.
         do iq=1,n
            v(ip,iq)=0.
         enddo
         v(ip,ip)=1.
      enddo 
      do ip=1,n
         b(ip)=a(ip,ip) !initialize b and d to the diagonal of a.
         d(ip)=b(ip)
         z(ip)=0.  !this vector will accumulate terms of the form tapq
      enddo        !as in equation (11.1.14).
      nrot=0
      do i=1,50
         sm=0.
         do ip=1,n-1  !sum off-diagonal elements.
         do iq=ip+1,n
            sm=sm+dabs(a(ip,iq))
         enddo 
         enddo 
         if(sm==0.)return   
         !===================!
         ! sucessful return  !
         !===================!
         if(i.lt.4)then       
            tresh=0.2*sm/n**2  !...on the first three sweeps.
         else
            tresh=0.           !...thereafter.
         endif
         do ip=1,n-1
         do iq=ip+1,n
            g=100.*dabs(a(ip,iq))
            !after four sweeps, skip the rotation if the off-diagonal element is small.
            if((i.gt.4).and.(dabs(d(ip))+g==dabs(d(ip)))
     *         .and.(dabs(d(iq))+g==dabs(d(iq))))then
               a(ip,iq)=0.
            else if(dabs(a(ip,iq)).gt.tresh)then
               h=d(iq)-d(ip)
               if(dabs(h)+g==dabs(h))then
                  t=a(ip,iq)/h          !t = 1/(2*theta)
               else
                  theta=0.5*h/a(ip,iq)  !equation (11.1.10).
                  t=1./(dabs(theta)+dsqrt(1.+theta**2))
                  if(theta.lt.0.)t=-t
               endif
               c=1./dsqrt(1+t**2)
               s=t*c
               tau=s/(1.+c)
               h=t*a(ip,iq)
               z(ip)=z(ip)-h
               z(iq)=z(iq)+h
               d(ip)=d(ip)-h
               d(iq)=d(iq)+h
               a(ip,iq)=0.
               do j=1,ip-1  !case of rotations 1 = j < p.
                  g=a(j,ip)
                  h=a(j,iq)
                  a(j,ip)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo 
               do j=ip+1,iq-1 !case of rotations p < j < q.
                  g=a(ip,j)
                  h=a(j,iq)
                  a(ip,j)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo 
               do j=iq+1,n !case of rotations q < j = n.
                  g=a(ip,j)
                  h=a(iq,j)
                  a(ip,j)=g-s*(h+g*tau)
                  a(iq,j)=h+s*(g-h*tau)
               enddo
               do j=1,n
                  g=v(j,ip)
                  h=v(j,iq)
                  v(j,ip)=g-s*(h+g*tau)
                  v(j,iq)=h+s*(g-h*tau)
               enddo 
               nrot=nrot+1
            endif
         enddo 
         enddo 
         do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip) !update d with the sum of tapq ,
            z(ip)=0.    !and reinitialize z.
         enddo
      enddo 

      write(*,*) 'too many iterations in jacobi, give up'
      ising=1

      return
      end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+                                                           +
c+     calculate stress due to dislocation line              +
c+                                                           +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine stress_dislocation(p1,p2,t,b,mu,nu,p,sm)
      implicit none
      integer i,j,k,l,m,n,ip
      real(8) p1(3),p2(3),p(3),t(3),b(3),mu,nu,sm(3,3),sv(6)
      real(8) r(2,3),rl(2,3),rp(2,3),ry(2,3)
      real(8) c_r(2),c_rl(2),c_rp(2),c_ry(2)
      real(8) stf(2,3,3,3)
      real(8) x1,x2,x3,mx1(3,3),mx2(3,3),vx(3)
      real(8) eij(3,3),eijk(3,3,3)
      real(8) c1,c2,c3,c4,c5
      real(8) z1,z2,z3,z4,z5,z6,z7

      eij(1,:)=[1,0,0]
      eij(2,:)=[0,1,0]
      eij(3,:)=[0,0,1]
      eijk=0
      eijk(1,2,3)= 1
      eijk(1,3,2)=-1
      eijk(2,1,3)=-1
      eijk(2,3,1)= 1
      eijk(3,1,2)= 1
      eijk(3,2,1)=-1

      r(1,:)=p-p1
      c_r(1)=dsqrt(sum(r(1,:)**2))
      c_rl(1)=dot_product(r(1,:),t)
      rp(1,:)=r(1,:)-c_rl(1)*t;          
      ry(1,:)=r(1,:)+c_r(1)*t
      c_ry(1)=dsqrt(sum(ry(1,:)**2))

      r(2,:)=p-p2
      c_r(2)=dsqrt(sum(r(2,:)**2))
      c_rl(2)=dot_product(r(2,:),t)
      rp(2,:)=r(2,:)-c_rl(2)*t;          
      ry(2,:)=r(2,:)+c_r(2)*t
      c_ry(2)=dsqrt(sum(ry(2,:)**2))

c      print*,c_r
c      print*,c_rl
c      print*,c_ry
c      read*

      stf=0
      do ip=1,2
         c1=mu/(3.14*c_ry(ip)**2)
         c2=1/(1-nu)
         vx=0
         do i=1,3
         do j=1,3
         do k=1,3
            vx(k)=vx(k)+eijk(k,i,j)*ry(ip,i)*t(j)
         enddo
         enddo
         enddo
         vx=vx/(2*(1-nu))
         c4=2/c_ry(ip)**2
         c5=c_rl(ip)/c_r(ip)
         do i=1,3
         do j=1,3
         do k=1,3
            z1=( dot_product(eijk(i,k,:),ry(ip,:))*t(j)
     &          +dot_product(eijk(j,k,:),ry(ip,:))*t(i))/2
            z2=( dot_product(eijk(i,k,:),t       )*ry(ip,j)
     &          +dot_product(eijk(j,k,:),t       )*ry(ip,i))/2
            z3=eij(i,j)+t(i)*t(j)
            z4=rp(ip,i)*ry(ip,j)+rp(ip,j)*ry(ip,i)
            z5=ry(ip,i)*ry(ip,j)
            stf(ip,i,j,k)=c1*(z1-c2*z2-vx(k)*(z3+c4*(z4+c5*z5)))         
         enddo
         enddo
         enddo
      enddo
      do i=1,3
      do j=1,3
         sm(i,j)=dot_product(stf(2,i,j,:)-stf(1,i,j,:),b)
      enddo
      enddo

      return
      end




c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   crystal plasticity fast fourier transformation code         +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine smooth(nx,ny,nz,l_size,rve_x)
      implicit none
      integer i,j,k,k1,ising
      integer nx,ny,nz,ix,iy,iz,mx,my,mz
      integer, parameter :: ax=0
      integer, parameter :: ay=0
      integer, parameter :: az=0
      real(8) fct_smooth,l_size
      real(8) wsfft((nx+ax*2)*(ny+ay*2)*(nz+az*2)*2)
      real(8) rve_v((nx+ax*2),(ny+ay*2),(nz+az*2))
      real(8) rve_r((nx+ax*2),(ny+ay*2),(nz+az*2))
      real(8) rve_i((nx+ax*2),(ny+ay*2),(nz+az*2))
      real(8) rve_x(nx,ny,nz)
      real(8) x1,x2,x3,fqc(3)

      mx=nx+ax*2
      my=ny+ay*2
      mz=nz+az*2

!      fct_smooth = 4*3.14**2*10.d0
!      fct_smooth = 4*3.14**2*20.d0
!      fct_smooth = 4*3.14**2*40.d0
!      fct_smooth = 4*3.14**2*100.d0
!      fct_smooth = 4*3.14**2*50.d0
!      fct_smooth = 4*3.14**2*30.d0

      fct_smooth = 4*3.14**2*60.d0
!      fct_smooth = 4*3.14**2*40.d0
!      fct_smooth = 4*3.14**2*40.d0*1.d-6

!      fct_smooth = 4*3.14**2*4.d19*l_size**2

c-----transfer stress to frequency space
      rve_v=0
      rve_v( ax+1:ax+nx, ay+1:ay+ny, az+1:az+nz ) = rve_x 

      wsfft=0               
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; wsfft(k1)=rve_v(ix,iy,iz)
         k1=k1+1; wsfft(k1)=0               
      enddo
      enddo
      enddo
      call fourn(wsfft,[mx,my,mz],3,1)
      wsfft=wsfft/(mx*my*mz)
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; rve_r(ix,iy,iz)=wsfft(k1) !--real part
         k1=k1+1; rve_i(ix,iy,iz)=wsfft(k1) !--image part
      enddo
      enddo
      enddo

c-----calculate equilibrium value in frequency space
      do ix=1,mx
      do iy=1,my
      do iz=1,mz
         if(ix<=mx/2+1)then
            fqc(1) = (ix-1.)/mx
         else
            fqc(1) = (ix-(mx+1.))/mx
         endif
         if(iy<=my/2+1)then
            fqc(2) = (iy-1.)/my
         else
            fqc(2) = (iy-(my+1.))/my
         endif
         if(iz<=mz/2+1)then
            fqc(3) = (iz-1.)/mz
         else
            fqc(3) = (iz-(mz+1.))/mz
         endif
         x1= 1 + fct_smooth*sum(fqc**2)
         rve_r(ix,iy,iz) = rve_r(ix,iy,iz)/x1
         rve_i(ix,iy,iz) = rve_i(ix,iy,iz)/x1
      enddo
      enddo
      enddo

c-----transfer new value from frequency space to physical space
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; wsfft(k1)=rve_r(ix,iy,iz)
         k1=k1+1; wsfft(k1)=rve_i(ix,iy,iz)
      enddo
      enddo
      enddo
      call fourn(wsfft,[mx,my,mz], 3, -1)
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; rve_v(ix,iy,iz)=wsfft(k1) 
         k1=k1+1
      enddo
      enddo
      enddo

      rve_x = rve_v( ax+1:ax+nx, ay+1:ay+ny, az+1:az+nz ) 
c
      return
      end

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   output date to ovito                                        +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine output(nx,ny,nz,sigr)
         implicit none
         integer ix,iy,iz,nx,ny,nz
         real(8) sigr(6,nx,ny,nz),v3(3)

         open(888,file='../stress_in.chkpt')

         write(888,'("#f a 1 1 1 3 1 1")')
         write(888,'("#c np ng nl x y z sig11 sig22 sig33 
     &   sig12 sig13 sig23")')
         write(888,'("#x 1 0 0")')
         write(888,'("#y 0 1 0")')
         write(888,'("#z 0 0 1")')
         write(888,'("#e")')

         do ix=1,nx
         do iy=1,ny
         do iz=1,nz
            v3=[ix,iy,iz]*2.d0  
            write(888,'(3i10,100e12.3)') 
     &         1,1,1,v3,sigr(1,ix,iy,iz),sigr(2,ix,iy,iz),
     &         sigr(3,ix,iy,iz),sigr(4,ix,iy,iz),sigr(5,ix,iy,iz),
     &         sigr(6,ix,iy,iz)
         enddo
         enddo
         enddo
         
         close(888)

         return
         end

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> @author 
   !> philipp engels, anxin ma
   ! description: 
   !> returns frobenius norm of matrix 
   !> @param[in]   matrix    matrix
   !> @param[in]   rows
   !> @param[in]   columns
   !> @param[out]  norm      frobenius norm    
c--------------------------------------------------------------------------
      subroutine icams_matrixnorm(matrix,rows,columns,norm)
c--------------------------------------------------------------------------
        implicit none
        integer :: i,j,rows, columns
        real(8), intent(in) :: matrix(rows,columns)
        real(8), intent(out) :: norm  
        norm = 0.d0
        do i=1,rows
          do j=1,columns 
             norm=norm+matrix(i,j)**2.d0 
          enddo 
        enddo
        norm=dsqrt(norm)
        return 
      end subroutine
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> @author 
   !> martin boeff
   ! description: 
   !> creates identity matrix of dimension n
   !> @param[in]   n                dimension of matrix
   !> @param[out]  mi               identity matrix
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine identitymatrix(n,mi)
c--------------------------------------------------------------------------
       implicit none
       integer,intent(in) :: n
       real(8),intent(out) :: mi(n,n)
       !local variables
       integer :: i            !loop integer

       !start assembly
       mi=0.0
       do i=1,n
         mi(i,i)=1.0                        
       enddo
       return
      end

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine thermal_ex(dther,ndi,ntens,dtemp,temp,ther)
       implicit none
       integer ndi
       integer :: i
       integer ntens
       real(8) dtemp,temp           ! increntment of temparature
       real(8) dther(ntens)      ! thermal deformation gradient
       real(8) ther,mf    
       
c	 calculate  thermal strainss
       if(ntens .eq. 3) then
         do i=1,ndi                  
            dther(i)=ther*(dtemp)
         enddo           
         dther(ntens)=0.d0
       else
         do i=1,ndi                 
            dther(i)=ther*(dtemp)
         enddo  
         do i=ndi+1,ntens
            dther(i)=0.d0
         enddo
       endif
       return
      end
ccccccccc calculate stress
      subroutine new_stress(dth,ds,ndi,ntens,st)
       implicit none
       integer ndi
       integer :: i
       integer :: j
       integer ntens
       real(8) dth(ntens)      ! thermal deformation gradient
       real(8) ds(ntens,ntens)  
       real(8) st(ntens)

       do  i=1,ntens
         do  j=1,ntens
            st(j)=st(j)+ds(j,i)*(dth(i))
         enddo
       enddo
       return
      end
c****************************** calculate work plastic deformation for CP model***************************************
c****************************** The following definition is used to couple the crystal plasticity and phase field models ***********************
C Wcpdot= sigma(total slip system)(shear stress * plastic slip rate)
c      subroutine cal_workpl(wpls)
c      implicit none
c      real wpls
c      integer :: i
c      integer :: j 
c      integer Nslp_mx
c      real(8) tau(Nslp_mx)
c     real(8) dgmdt(Nslp_mx)
c      integer :: N_slip = 12
c      real(8) IVB_bk(Nslp_mx)

c      do i=1,N_slip
c         ss = ss + dabs(tau(i)-IVB_bk(i))*dgmdt(i)
c      enddo