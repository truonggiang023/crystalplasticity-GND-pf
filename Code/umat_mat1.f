C******************************************************HETVAL*********************************************************
      subroutine umat_mat1(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c
      include 'aba_param.inc'
c
      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3),jstep(4)

      dimension eelas(ntens),eplas(ntens),flow(ntens),olds(ntens),oldpl(ntens), etherm(ntens)
      dimension dtherm(ntens), table(2,2000), stressl(ntens), stranl(ntens), ddsddel(ntens,ntens)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0)
      data newton,toler/10,1.d-6/

            e = props(1)
            xnu = props(2)
            b7 = props(3) ! σ_0 
            b8 = props(4)  ! σ_delta
            b9 = props(5) ! m
            b10 = props(6) ! alpha
            b11 = props(7) ! shear modulus
            b12 = props(8) ! b
            b13 = props(9) ! k
            b14 = props(10) !L
            k_t = props(11) ! transformation ratio
            alfa_aust = props(12) ! thermal coeficient austenite to martensite
            steel =   props(13) ! thermal coeficient
            ms = props(14) ! start tranformation from austenite to martensite
            mf = props(15) ! finsish tranformation martensite
            phase = props(16) ! define phase
c --------------------- define phase ferrite = 1, martensite = 2 --------------------------------
c --------------------- phase field damage model -------------------------------------------------
            xl=props(17) ! phase field length scale
            dd=(1.d0 - predef(1))**2 ! interface indicator functionss
            gc= dd*(props(20) - props(19)) + props(19)  ! fracture toughness at interface in phase                  
            kflags=int(props(18)) ! solution flag (0: monolithic, 1: staggered)
            phi= temp+dtemp + 1d-03
            psit=statev(170)
            g = (1.d0-phi)**2 ! degradation function

c-------------------------------------
c--------------------------------------
            if (phase.eq.1) then
             ther = steel
            else if (phase.eq.2) then
                  if ((predef(1)+dpred(1)).lt.mf) then
                        ther = alfa_aust
                  else 
                        ther = steel
                  endif
            endif
c -----------------------------------
            xmu= e/(2.d0*(1.d0+xnu))                 ! #shear modulus:   ��=e/[2*(1+��)]
            xlm= xnu*e/((1.d0+xnu)*(1.d0-2.d0*xnu))  ! #lames constant: ��=e��/[(1+��)(1-2��)]
            xval=e/(1.d0-xnu**2)    
            eg2= e/(1.d0+xnu)
            eg = eg2/2.d0
            ebulk3=e/(1.d0-2.d0*xnu)
c----------------
            if(ntens .eq. 3) then
             do i = 1, ntens
    		  do j = 1, ntens
    			 ddsdde(i,j)= zero
    		  enddo
    		 enddo
             do i=1, ndi
              do j=1, ndi
                  ddsdde(i,j)= xval*xnu
              enddo
              ddsdde(i,i)= xval*1.d0
             enddo
c     #   ddsdde(3,3)
             ddsdde(ntens,ntens)= 1.d0*xmu
              else
c     #   ddsdde(1:3, 1:3)
             do i=1, ndi
              do j=1, ndi
                  ddsdde(i,j)=xlm
              enddo
             ddsdde(i,i)=xlm+2.d0*xmu
             enddo
c     #   ddsdde(4:6, 4:6)
             do i=ndi+1, ntens
                  ddsdde(i,i)=xmu
             enddo
            endif
c -------------------------------------------------------------------
c	 calculate  thermal strain
            if ((predef(1)+dpred(1)).lt.mf) then
                  if(ntens .eq. 3) then
                        do i=1,ndi                  
                              dtherm(i)=ther*(dpred(1))
                        enddo           
                        dtherm(ntens)=zero
                  else
                        do i=1,ndi                 
                              dtherm(i)=ther*(dpred(1))
                        enddo
                        do i=ndi+1,ntens
                              dtherm(i)=zero
                        enddo
                  endif
            else
                  do i=1,ntens
                        dtherm(i)=zero
                  enddo
            endif      
c ------------------------------------------------------------------
c calculate stress
            do  i=1,ntens
             do  j=1,ntens
                  stress(j)=stress(j)+ddsdde(j,i)*(dstran(i) + dtherm(i))
             enddo
            enddo
ccccccc begining a new step tensile test ccccccccc            
c     get  linear elastic stiffness tensor
            stranl = stran+dstran
c     the undamaged cauchy stress
            stressl = matmul(ddsdde,stranl)
c     recover elastic and plastic strains
            do  i=1,ntens
                  eelas(i)=statev(i+100)+dstran(i)
                  etherm(i)=statev(i+ntens+100)+dtherm(i)
                  eplas(i)=statev(i+2*ntens+100)
                  statev(i+5*ntens+100)=eplas(i)
            enddo
            eqplas=statev(1+3*ntens+100)
            statev(150) = eqplas
c calculate equivalent von mises stress
            if (nprops.gt.2.and.props(3).gt.0.0) then
             if(ntens .eq. 3) then
                  smises = ((stress(1)-stress(2)))**2 +(stress(2))**2+((stress(1)))**2
                  smises = smises + 6.d0*stress(ntens)**2
                  smises = sqrt(smises/2.d0)
             else
                  smises = ((stress(1)-stress(2)))**2+((stress(2)-stress(3)))**2+((stress(1)-stress(3)))**2
                  do i=ndi+1,ntens
                        smises = smises + 6.d0*stress(i)*stress(i)
                  enddo
                  smises = sqrt(smises/2.d0)
             endif
c             print*,smises
c             smisesther = statev(150)
c             if (time(2) .eq. 1.d0) then
c                        smisesther = smises
c             endif     
C             statev(150) = smisesther
c get yeild stress from the specified hardening curve
             nvalue=2000
             if (kstep.eq.1) then
                  do i= 0,nvalue-1
                  eqplas_hard = i*0.005
                  if ((phase.eq.1).and.((predef(1)+dpred(1)).ge.mf)) then
                        flow_curve=b7+b8+b10*b9*b11*sqrt(b12)*sqrt((1-exp(-b9*b13*eqplas_hard))/(b14*b13))
                        table(1,i+1)=flow_curve
                        table(2,i+1)=eqplas_hard
                  else if ((phase.eq.1).and.((predef(1)+dpred(1)).lt.mf)) then
                        flow_curve1f=b7+b8+b10*b9*b11*sqrt(b12)*sqrt((1-exp(-b9*b13*eqplas_hard))/(b14*b13))
                        flow_curve= flow_curve1f*(0.00131*(predef(1)+dpred(1))+1.026)
                        table(1,i+1)=flow_curve 
                        table(2,i+1)=eqplas_hard                                    
                  else if ((phase.eq.2).and.((predef(1)+dpred(1)).ge.mf)) then
                        flow_curve=b7+b8+b10*b9*b11*sqrt(b12)*sqrt((1-exp(-b9*b13*eqplas_hard))/(b14*b13))
                        table(1,i+1)=flow_curve 
                        table(2,i+1)=eqplas_hard
                  else if ((phase.eq.2).and.((predef(1)+dpred(1)).lt.mf)) then
c                  flow_curve=600 + 1200*((eqplas_hard)**(0.4))                                      
                        flow_curve1a=650 + 1200*((eqplas_hard)**(0.4))
                        flow_curve2a=flow_curve1a*(0.00131*(predef(1)+dpred(1))+1.026)
                        flow_curve1m=b7+b8+b10*b9*b11*sqrt(b12)*sqrt((1-exp(-b9*b13*eqplas_hard))/(b14*b13))
                        flow_curve2m= flow_curve1m*(0.00028*(predef(1)+dpred(1))+1.026)               
c                 flow curve martensite differenf vs ferrite 
                        f_m = 1.0 - exp(-k_t*(ms+(predef(1)+dpred(1))))
                        flow_curve = ((1.0-f_m)*flow_curve2a + f_m*flow_curve2m)
                        table(1,i+1)=flow_curve
                        table(2,i+1)=eqplas_hard
                  endif
                  enddo
             else
                  do i= 0,nvalue-1
                  eqplas_hard = i*0.005
                  if ((phase.eq.1).and.((predef(1)+dpred(1)).ge.mf)) then
                        flow_curve=b7+b8+b10*b9*b11*sqrt(b12)*sqrt((1-exp(-b9*b13*eqplas_hard))/(b14*b13))
                        table(1,i+1)=flow_curve
                        table(2,i+1)=eqplas_hard
                  else if ((phase.eq.1).and.((predef(1)+dpred(1)).lt.mf)) then
                        flow_curve1f=b7+b8+b10*b9*b11*sqrt(b12)*sqrt((1-exp(-b9*b13*eqplas_hard))/(b14*b13))
                        flow_curve= flow_curve1f
                        table(1,i+1)=flow_curve 
                        table(2,i+1)=eqplas_hard                                    
                  else if ((phase.eq.2).and.((predef(1)+dpred(1)).ge.mf)) then
                        flow_curve=b7+b8+b10*b9*b11*sqrt(b12)*sqrt((1-exp(-b9*b13*eqplas_hard))/(b14*b13))
                        table(1,i   +1)=flow_curve 
                        table(2,i+1)=eqplas_hard
                  else if ((phase.eq.2).and.((predef(1)+dpred(1)).lt.mf)) then
c                  flow_curve=600 + 1200*((eqplas_hard)**(0.4))                                      
                        flow_curve1a=650 + 1200*((eqplas_hard)**(0.4))
                        flow_curve2a=flow_curve1a
                        flow_curve1m=b7+b8+b10*b9*b11*sqrt(b12)*sqrt((1-exp(-b9*b13*eqplas_hard))/(b14*b13))
                        flow_curve2m= flow_curve1m            
c                 flow curve martensite differenf vs ferrite 
                        f_m = 1.0 - exp(-k_t*(ms+(predef(1)+dpred(1))))
                        flow_curve = ((1.0-f_m)*flow_curve2a + f_m*flow_curve2m)
                        table(1,i+1)=flow_curve
                        table(2,i+1)=eqplas_hard
                  endif
                  enddo             
             endif      
             nvalue=2000
C             print*, table
             call uhard(syiel0,hard,eqplas,table,nvalue,predef(1)+dpred(1))
c             write(6) 'the value of i is: ', hard
c get yeild stress from the specified hardening curve
c sf = sy*(1.d0+e*eqplas/sy)**xn
             if ((smises.gt.(1.d0+toler)*syiel0)) then !
c calculate the flow dicrections
              if(ntens .eq. 3) then
                  sh = (stress(1)+stress(2))/3.d0
                  onesy=one/smises
                  do i=1,ndi
                        flow(i)=onesy*(stress(i)-sh)
                  enddo
                  flow(ntens)=stress(ntens)*onesy
c ----------------------------------------------
                  syield=syiel0
                  deqpl = 0.d0
                  do kewton = 1,newton
                        rhs = smises-(3.d0*eg)*deqpl-syield
                        deqpl = deqpl + rhs/((3.d0*eg)+hard)
                        call uhard(syield,hard,eqplas+deqpl,table,nvalue,predef(1)+dpred(1))
                        if (abs(rhs).lt.toler*syiel0) goto 90
                  enddo
                  if (kewton.eq.newton) write(7,*)'warning:plasticity loop failed'
90                continue
                  effhrd=(3.d0*eg*hard)/(3.d0*eg+hard)
              else
                  sh = (stress(1)+stress(2)+stress(3))/3.d0
                  onesy=one/smises
                  do i=1,ndi
                        flow(i)=onesy*(stress(i)-sh)
                  enddo
                  do i=ndi+1,ntens
                        flow(i)=stress(i)*onesy
                  enddo
                  syield=syiel0
                  deqpl = 0.d0
                  do kewton = 1,newton
                        rhs = smises-(3.d0*eg)*deqpl-syield
                        deqpl = deqpl + rhs/((3.d0*eg)+hard)
                        call uhard(syield,hard,eqplas+deqpl,table,nvalue,predef(1)+dpred(1))
                        if (abs(rhs).lt.toler*syiel0) goto 100
                  enddo
                  if (kewton.eq.newton) write(7,*)'warning:plasticity loop failed'
100               continue
                  effhrd=3.d0*eg*hard/(3.d0*eg+hard)
              endif
c ---------------------------------------
c update stresses and strains
              if(ntens .eq. 3) then
                  do i=1,ndi
                        stress(i)=flow(i)*syield+sh
                        eplas(i)=eplas(i)+3.0d0*flow(i)*deqpl/2.0d0
                        eelas(i)=eelas(i)-3.0d0*flow(i)*deqpl/2.0d0 
                  enddo
                        stress(ntens)=flow(ntens)*syield
                        eplas(ntens)=eplas(ntens)+3.0d0*flow(i)*deqpl 
                        eelas(ntens)=eelas(ntens)-3.0d0*flow(i)*deqpl 
                  eqplas=eqplas+deqpl
                  spd=deqpl*(syiel0+syield)/2.0d0
              else
                  do i=1,ndi
                        stress(i)=flow(i)*syield+sh
                        eplas(i)=eplas(i)+3.0d0*flow(i)*deqpl/2.0d0
                        eelas(i)=eelas(i)-3.0d0*flow(i)*deqpl/2.0d0 
                  enddo
                  do i=ndi+1,ntens
                        stress(i)=flow(i)*syield
                        eplas(i)=eplas(i)+3.0d0*flow(i)*deqpl
                        eelas(i)=eelas(i)-3.0d0*flow(i)*deqpl 
                  enddo
                  eqplas=eqplas+deqpl
                  spd=deqpl*(syiel0+syield)/2.0d0
              endif
c ---------jacobian
              if(ntens .eq. 3) then
                  effg=eg*syield/smises
                  effg2=2.0d0*effg
                  effg3=3.0d0*effg2/2.0d0
                  efflam=(ebulk3-effg2)/3.0d0
                  do i=1,ndi
                   do  j=1,ndi
                              ddsdde(j,i)=efflam
                   enddo
                              ddsdde(i,i)=effg2+efflam
                  enddo
                   ddsdde(ntens,ntens)=effg
                  do  i=1,ntens
                   do  j=1,ntens
                              ddsdde(j,i)=ddsdde(j,i)+flow(j)*flow(i)*(effhrd-effg3)
                   enddo
                  enddo            
              else
                  effg=eg*syield/smises
                  effg2=2.0d0*effg
                  effg3=3.0d0*effg2/2.0d0
                  efflam=(ebulk3-effg2)/3.0d0
                  do i=1,ndi
                   do  j=1,ndi
                              ddsdde(j,i)=efflam
                   enddo
                              ddsdde(i,i)=effg2+efflam
                  enddo
                  do  i=ndi+1,ntens
                              ddsdde(i,i)=effg
                  enddo
                  do  i=1,ntens
                   do  j=1,ntens
                              ddsdde(j,i)=ddsdde(j,i)+flow(j)*flow(i)*(effhrd-effg3)
                   enddo
                  enddo
              endif
c     store strains in state variable array
             do i=1,ntens
                  statev(i+100)= eelas(i)
                  statev(i+2*ntens+100)=eplas(i)
                  statev(i+ntens+100)=etherm(i)
             enddo
                  statev(1+3*ntens+100)=eqplas
                  statev(100)= eqplas
                  statev(106) = eqplas
             endif
            endif
c   calculate stress equivalance plastic
            if (ntens .eq. 3) then
                  smises1 = ((stress(1)-stress(2)))**2 +(stress(2))**2+((stress(1)))**2
                  smises1 = smises1 + 6.d0*stress(ntens)**2
                  smises1 = sqrt(smises1/2.d0)
            else
                  smises1 = ((stress(1)-stress(2)))**2+((stress(2)-stress(3)))**2+((stress(1)-stress(3)))**2
                  do i=ndi+1,ntens
                        smises1 = smises1 + 6.d0*stress(i)*stress(i)
                  enddo
                  smises1 = sqrt(smises1/2.d0)
            endif
            statev(130) = smises1
ccccccc after going down to 0.0 cccccccccccccccccccccccccccccccccc
c     compute the strain energy density
c            if (kstep.eq.3.d0) then
            if (statev(1+3*ntens+100).le.0) then
                        psi=0.d0
                        do i=1,ntens
                              psi=psi+abs(stressl(i)*(stranl(i))*0.5d0)
                        enddo
                        wp=0
            else
                        psi = statev(160)
                        wp = statev(131)
                        wp1 = 0
                        do i=1,ntens
                              wp1= wp1 + abs(stress(i)*abs(statev(i+2*ntens+100)- statev(i+5*ntens+100)))
                        enddo
                        wp = wp + wp1
                        psi=psi+wp
            endif
            h=abs(max(psit,psi))
            stress=stress*g
            ddsdde=ddsdde*g
            statev(170)=h
            statev(124)=h
            statev(125)=xl
            statev(126)=gc
            statev(127)=dd
            if (kflags.eq.1) statev(124)=psit ! staggered scheme
            statev(131) = wp
            statev(160) = psi
c            endif
            statev(169)= phase
      return
      end
c --------------------------------------
c --------------------------------------
      subroutine uhard(syield,hard,eqplas,table,nvalue,temp)
c
      include 'aba_param.inc'
      dimension table(2,nvalue)
c     
c   if more than one entry, search table
c

c    set yield stress to last value of table, hardening to zero
      syield=table(1,nvalue)
      hard=0.0      
      if(nvalue.gt.1) then
        do i=1,nvalue-1
           eqpl1=table(2,i+1)
           if(eqplas.lt.eqpl1) then
             eqpl0=table(2,i)
             if(eqpl1.le.eqpl0) then
                write(6,1)
 1              format(//,30x,'***error - plastic strain must be ',
     1                 'entered in ascending order')
                call xit
              endif
c
c           current yield stress and hardening
c
            deqpl=eqpl1-eqpl0
c            if (temp.ge.(528.d0)) then
c                  syiel0=table(1,i)*0.7
c                  syiel1=table(1,i+1)*0.7
c            else if (temp.le.(298.d0)) then
c                 syiel0=(table(1,i))
c                 syiel1=(table(1,i+1))
c            else if ((temp.gt.(298.d0)).and.(temp.lt.(528.d0))) then
c                  syiel0=(table(1,i))+((temp-298.d0)*(((table(1,i)*0.7) - ((table(1,i))))/(230.d0)))
c                  syiel1=(table(1,i+1))+((temp-298.d0)*(((table(1,i+1)*0.7) - ((table(1,i+1))))/(230.d0)))
c            endif
            syiel0=table(1,i)
            syiel1=table(1,i+1)
            dsyiel=syiel1-syiel0
            hard=dsyiel/deqpl
            syield=syiel0+(eqplas-eqpl0)*hard
            goto 20
            endif
        enddo
 20         continue
      endif
      return
      end