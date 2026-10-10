        include 'umat_mat0.f'
        include 'umat_mat1.f'
C*******************************************************************************************
Ccccccccccccccccccccccccccccccccc
        subroutine hetval(cmname,temp,time,dtime,statev,flux,predef,dpred,kstep)

        include 'aba_param.inc'
            
        character*80 cmname
            
        dimension temp(2),statev(*),predef(*),time(2),flux(2),dpred(*)

c        if (predefined.ge.0.d0) then    
        phi=temp(1)
        h=statev(124)
        xl=statev(125)
        gc=statev(126)
        wp = statev(131)
        a = wp
        if (a.gt.0) then
                  flux(1)=-((phi/(xl**2))-((2.d0*(1.d0-phi)*(h+a))/(gc*xl)))
                  flux(2)=-((1.d0/(xl**2))+2.d0*(h+a)/(gc*xl))
        else
                  flux(1)=-((phi/(xl**2))-((2.d0*(1.d0-phi)*(h))/(gc*xl)))
                  flux(2)=-((1.d0/(xl**2))+2.d0*(h)/(gc*xl))
        endif
c       endif
        return
        end
c********************* umat subroutine **********************
        subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,stran,dstran,
     2 time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,celent,
     4 dfgrd0,dfgrd1,noel,npt,kslay,kspt,kstep,kinc);
c
        include 'aba_param.inc'

        character*80 cmname;

        dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
c
        if (cmname(1:13).ne."MARTENSITE_J2") then
                 call umat_mat0(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,stran,dstran,
     2 time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,celent,
     4 dfgrd0,dfgrd1,noel,npt,kslay,kspt,kstep,kinc);
        else if(cmname(1:13).eq."MARTENSITE_J2") then
                 call umat_mat1(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,stran,dstran,
     2 time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,celent,
     4 dfgrd0,dfgrd1,noel,npt,kslay,kspt,kstep,kinc);
        endif
        return
        end

