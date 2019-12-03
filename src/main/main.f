ccc   calculates CEP cross section
      function cs(rarr,wgt)
      implicit double precision (a-z)
      complex*16 wt(10)
      double precision rarr(10),ran(5)
      integer i,p,icut,outl

      include 'gencuts.f'
      include 'pi.f'
      include 'unweighted.f'
      include 'bin.f'
      include 'x.f'
      include 'vars.f'
      include 'mom.f'
      include 'mres.f'
      include 'mp.f'
      include 'eff.f'
      include 'wmax.f'
      include 'alpvars.f'
      include 'exp.f'
      include 'gax.f'
      include 'mn.f'
      include 'adecay.f'
      
      wtt=0d0

      mx=mres    

      mpp=mp
     
      r2=rarr(2)
      r3=rarr(3)
      r4=rarr(4)
      r5=rarr(5)

      r1=rann2()
         
      phi1=2d0*pi*r1
      phi2=2d0*pi*r2+phi1
  
      ptmax=dsqrt(3d0)           
      ptmin=0d0
      
      xgmin=(mx/rts)**2
      
      ypmax=dlog(xgmin**2*mpp**2+ptmax**2)
      ypmin=dlog(xgmin**2*mpp**2+ptmin**2)

      yppmax=dlog(xgmin**2*mn**2+ptmax**2)
      yppmin=dlog(xgmin**2*mn**2+ptmin**2)
      
      yp=(ypmax-ypmin)*r3+ypmin
      ypp=(yppmax-yppmin)*r4+yppmin
      
      pt1sq=dexp(yp)-xgmin**2*mpp**2
      pt2sq=dexp(ypp)-xgmin**2*mn**2
      
      if(pt1sq.lt.0d0)then
         pt1sq=0d0
      endif
      
      if(pt2sq.lt.0d0)then
         pt2sq=0d0
      endif

      pt2x=dsqrt(pt2sq)*dcos(phi2)
      pt2y=dsqrt(pt2sq)*dsin(phi2)
      pt1x=dsqrt(pt1sq)*dcos(phi1)
      pt1y=dsqrt(pt1sq)*dsin(phi1)
         
      ptxx=(pt1x+pt2x)**2+(pt1y+pt2y)**2
      rmx=dsqrt(ptxx+mx**2)
      
      ymax1=ymax
      ymin1=ymin
      ycut=dlog(rts/rmx)

      ymin=-ycut
      ymax=ycut
      
      if(ycut.lt.0d0)goto 777
         
      if(ymin.gt.0d0.and.ycut.lt.ymin)goto 777
      if(ymax.lt.0d0.and.ycut.gt.ymax)goto 777
      
      if(ymax.gt.ycut) ymax=ycut
      if(-ymin.gt.ycut) ymin=-ycut

      ry=rarr(1)
      yx=ymin+(ymax-ymin)*ry
      
      wty=ymax-ymin
      
      ymin=ymin1
      ymax=ymax1
      
      x1=rmx*dexp(yx)/rts       ! photon 1 mom. fraction 
      x2=rmx*dexp(-yx)/rts      ! photon 2 mom. fraction

      x1t=mx*dexp(yx)/rts       ! photon 1 mom. fraction 
      x2t=mx*dexp(-yx)/rts      ! photon 2 mom. fraction
     
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      aa1=(1d0-x1)*rts/dsqrt(2d0)
      aa2=(1d0-x2)*rts/dsqrt(2d0)
      cc1=0.5d0*(pt2sq+mn**2)
      cc2=0.5d0*(pt1sq+mp**2)

c     impose massive on-shell condition by solving
c                   p1+ + cc1/p2- = aa1
c                   p2- + cc2/p1+ = aa2 
      
      root1sq=(cc1-cc2-aa1*aa2)**2-4d0*cc2*aa1*aa2
      root2sq=(cc2-cc1-aa1*aa2)**2-4d0*cc1*aa1*aa2
      if(root1sq.le.0d0.or.root2sq.le.0d0)then
         wtt=0d0
         goto 777
      endif
      p1p=(cc2-cc1+aa1*aa2+dsqrt(root1sq))/(2d0*aa2)
      p2m=(cc1-cc2+aa1*aa2+dsqrt(root2sq))/(2d0*aa1)
      
      p1m=(pt1sq+mp**2)/(2d0*p1p)
      p2p=(pt2sq+mn**2)/(2d0*p2m)
      
      if(p1p.lt.0d0.or.p2m.lt.0d0)then
         wtt=0d0
         goto 777
      endif
      
      q(1,3)=pt1x
      q(2,3)=pt1y
      q(3,3)=(p1p-p1m)/dsqrt(2d0)
      q(4,3)=(p1p+p1m)/dsqrt(2d0)

      q(1,4)=pt2x
      q(2,4)=pt2y
      q(3,4)=(p2p-p2m)/dsqrt(2d0)
      q(4,4)=(p2p+p2m)/dsqrt(2d0)

      do i=1,4
         q(i,5)=q(i,1)+q(i,2)-q(i,3)-q(i,4)
      enddo

      if((q(4,5)**2-q(3,5)**2-q(2,5)**2-q(1,5)**2).lt.0d0)then
         wtt=0d0
         goto 777
      endif

cccccccccccccccccc
cccc decays
ccccccccccccccccccc

         wt2=1d0

         call axionps
         r6=rarr(6)
         r7=rarr(7)
         call twobodyn(1,10,11,12,r6,r7,wt2)

         if(adecay)then
         
         xdmax=lsh+dvol
         xdmin=lsh
                  
         xd=xdmin+(xdmax-xdmin)*r5
         
         call gamdecay
         
         if(xd/ald.gt.50d0)then
            wtdecay=0d0
         else
            wtdecay=1d0/ald*dexp(-xd/ald)
         endif
         
         endif
         
ccccccccccccccccccc  cuts ccccccccccccccccc

         neff0=neff0+1

         if(gencuts)then
            call cut(icut)
            if(icut.eq.0)goto 777
         endif
         
         neff=neff+1

ccccccccc
         
         call cscalc(pt1x,pt1y,pt2x,pt2y,wtt)
         
         wtt=wtt*wt2
         wtt=wtt*(ypmax-ypmin)*(yppmax-yppmin)
         wtt=wtt*(xgmin**2*mpp**2+pt1sq)*(xgmin**2*mn**2+pt2sq)
         wtt=wtt*wty

         beta=dsqrt(1d0-2d0*(mp**2+mn**2)/s+(mp**2-mn**2)**2/s**2)
         wtt=wtt/beta
         betat=dabs(q(3,3)/q(4,3)-q(3,4)/q(4,4))/2d0
         betat=betat*4d0*q(4,3)*q(4,4)/s
         wtt=wtt/betat
         wtt=wtt*2d0/mx
         wtt=wtt*pi/2d0/mx**3

         if(adecay)wtt=wtt*wtdecay*(xdmax-xdmin)

ccccccccccc         
         
 888     val=wtt*wgt
         if(bin)then
         call binit(val)
         endif
     
         if(calcmax)then
            if(wmax.lt.wtt*wgt*ren)then
               wmax=wtt*wgt*ren
               iw=iw+1
            endif
         endif

         if(unw)then
            runw=rann2()
            call unweight(wtt*wgt*ren,runw)
         endif

 777     cs=wtt
 
      return
      end

