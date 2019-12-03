ccccc EPA form factors
      subroutine formfacgam(t1,t2,out)
      implicit double precision(a-y)
      integer i1,i2

      include 'mn.f'
      include 'mp.f'
      include 'pi.f'
      include 'x.f'
      include 'vars.f'
      include 'mom.f'
      include 'ions.f'
      include 'btype.f'
      
ccccccc kinematics - get photon Q^2 values
      
      x1t=(q(4,1)-q(4,3))/q(4,1)
      x2t=(q(4,2)-q(4,4))/q(4,2)

      qsq1tt=(x1t**2*mp**2+t1)/(1d0-x1t)
      qsq2tt=(x2t**2*mn**2+t2)/(1d0-x2t)    

      qsq1=(q(4,1)-q(4,3))**2-(q(3,1)-q(3,3))**2-(q(2,1)-q(2,3))**2
     &        -(q(1,1)-q(1,3))**2
      qsqp1=(q(4,2)-q(4,4))**2-(q(3,2)-q(3,4))**2-(q(2,2)-q(2,4))**2
     &     -(q(1,2)-q(1,4))**2

      qsq=-qsq1
      qsqp=-qsqp1
     
      if(qsq.lt.0d0)qsq=qsq1tt
      if(qsqp.lt.0d0)qsqp=qsq2tt

ccccccccc proton form factor

      ge=1d0/(1d0+qsq/0.71d0)**4
      gm=ge*7.78d0
      fe=(4d0*mp**2*ge+qsq*gm)/(4d0*mp**2+qsq)
      fm=gm

      if(btype.eq.'prot')then
         fe1=fe
         fm1=fm
      elseif(btype.eq.'elec')then
         fe1=1d0
         fm1=1d0
      endif
      
ccccccc heavy ion form factor - 1512.03069 and refs. therein

      fe=tpint(1,dsqrt(qsqp))**2
      fm=0d0
      
      fm2=fm
      fe2=fe

ccccccc proton(ion) - photon vertices - evaluated using FORM
      
         ma=mx
         
         p1_q1=q(4,1)*(q(4,1)-q(4,3))-q(3,1)*(q(3,1)-q(3,3))
     &        -q(2,1)*(q(2,1)-q(2,3))-q(1,1)*(q(1,1)-q(1,3))

         p2_q2=q(4,2)*(q(4,2)-q(4,4))-q(3,2)*(q(3,2)-q(3,4))
     &        -q(2,2)*(q(2,2)-q(2,4))-q(1,2)*(q(1,2)-q(1,4))

         p1_q2=q(4,1)*(q(4,2)-q(4,4))-q(3,1)*(q(3,2)-q(3,4))
     &        -q(2,1)*(q(2,2)-q(2,4))-q(1,1)*(q(1,2)-q(1,4))

         p2_q1=q(4,2)*(q(4,1)-q(4,3))-q(3,2)*(q(3,1)-q(3,3))
     &        -q(2,2)*(q(2,1)-q(2,3))-q(1,2)*(q(1,1)-q(1,3))
         
         q1pcq2p=-(q(2,1)-q(2,3))*(q(1,2)-q(1,4))
     &        +(q(1,1)-q(1,3))*(q(2,2)-q(2,4))

cccccccc
         
      epacc =
     &  - 1./2.*qsq*qsqp*ma**4 - qsq*qsqp**2*ma**2 - 1./2.*qsq*qsqp**3
     &  - qsq**2*qsqp*ma**2 + qsq**2*qsqp**2 - 1./2.*qsq**3*qsqp

      epacd =
     & qsq*ma**4*mn**2 + 2*qsq*qsqp*ma**2*mn**2 + qsq*qsqp**2*mn**2 + 2
     & *qsq**2*ma**2*mn**2 - 2*qsq**2*qsqp*mn**2 + qsq**3*mn**2 - 4*
     & p2_q1*p2_q2*qsq*ma**2 - 4*p2_q1*p2_q2*qsq*qsqp - 4*p2_q1*p2_q2*
     & qsq**2 - 4*p2_q1**2*qsq*qsqp - 4*p2_q2**2*qsq**2

      epadc =
     & qsqp*ma**4*mp**2 + 2*qsqp**2*ma**2*mp**2 + qsqp**3*mp**2 + 2*qsq
     & *qsqp*ma**2*mp**2 - 2*qsq*qsqp**2*mp**2 + qsq**2*qsqp*mp**2 - 4*
     & p1_q1*p1_q2*qsqp*ma**2 - 4*p1_q1*p1_q2*qsqp**2 - 4*p1_q1*p1_q2*
     & qsq*qsqp - 4*p1_q1**2*qsqp**2 - 4*p1_q2**2*qsq*qsqp
             
      betap=dsqrt(1d0-2d0*(mp**2+mn**2)/s+(mp**2-mn**2)**2/s**2)
      epadd1=-4d0*s**2*(q1pcq2p**2)*betap**2
      
      outdd=-epadd1
      outcd=-epacd
      outdc=-epadc
      outcc=-epacc
 
      if(outdd.lt.0d0)outdd=0d0
      if(outcd.lt.0d0)outcd=0d0
      if(outdc.lt.0d0)outdc=0d0
      if(outcc.lt.0d0)outcc=0d0

      outtot=outdd*fe1*fe2+outcd*fm1*fe2+outdc*fe1*fm2+outdd*fm1*fm2
      outtot=outtot/2d0/s**2
      outtot=outtot/pi**2/137d0**2
      outtot=outtot/qsq**2/qsqp**2

      if(qsq.lt.0d0)outtot=0d0
      if(qsqp.lt.0d0)outtot=0d0

      out=dsqrt(outtot)
      
      return
      end
