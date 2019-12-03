ccc   generates two-body decay to particles of mass m1,m2
      subroutine axionps
c      (e,theta)
      implicit double precision(a-y)
      double precision pcm(4),pboo(4),px(4),pboo1(4),pcm1(4),
     &     pboo2(4)
      integer i
      common/phic/phic
      
      include 'mom.f'
      include 'pi.f'
      include 'mres.f'
      include 'mp.f'
      include 'vars.f'
      include 'x.f'
      include 'xt.f'
      include 'tcheck.f'
      include 'alpvars.f'
      include 'gax.f'
      
      do i=1,4
         pcm(i)=q(i,5)
         px(i)=-q(i,2)
      enddo
      px(4)=q(4,2)
      ein=dsqrt(px(4)**2-px(3)**2-px(2)**2-px(1)**2)
      
      call boost(ein,px,pcm,pboo)
      
      cost=pboo(3)/dsqrt(pboo(4)**2-mres**2)

      if(atheta.lt.1d-8)atheta=dasin(sint)
      if(dabs(cost).gt.1d0)then
         atheta=1d-10
      else
         atheta=dacos(cost)
      endif
      
ccccccc      
      
      ae=pboo(4)
      do i=1,4
         q(i,10)=pboo(i)
      enddo

cccccc

      taualp=64d0*pi/gax**2/mres**3
      agamma=ae/dsqrt(dabs(pboo(4)**2-pboo(3)**2-pboo(2)**2-pboo(1)**2)) ! gamma
      abeta=dsqrt(1d0-1d0/agamma**2) ! beta      
      ald=abeta*agamma*taualp   ! decay length
      
      convd=1.974d-16    ! convert to m
      ald=ald*convd

      aldmax=ald/ae*ebeam/abeta

      Return
      end
