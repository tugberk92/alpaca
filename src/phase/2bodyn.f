ccc   generates two-body decay to particles of mass m1,m2
      subroutine twobodyn(j,in,i1,i2,rcost,rphi,wt)
      implicit double precision(a-y)
      integer in,i1,i2
      integer i,j
      double precision pcm(4),pboo(4),px(4)
  
      include 'mom.f'
      include 'wtinit.f'
      include 'partonmom4.f'
      include 'pi.f'

      ein=dsqrt(q(4,in)**2-q(3,in)**2-q(2,in)**2-q(1,in)**2)

      cost=-1d0+2d0*rcost
      phi=2d0*pi*rphi

      pcm(4)=ein/2d0
      pcm(3)=pcm(4)*cost
      pcm(2)=pcm(4)*dsqrt(1d0-cost**2)*dsin(phi)
      pcm(1)=pcm(4)*dsqrt(1d0-cost**2)*dcos(phi)

      do i=1,4
         px(i)=q(i,in)
      enddo 

      call boost(ein,px,pcm,pboo)
      
      do i=1,4
         q(i,i1)=pboo(i)
         q(i,i2)=px(i)-q(i,i1)
      enddo
      
      wt=1d0

      return
      end
