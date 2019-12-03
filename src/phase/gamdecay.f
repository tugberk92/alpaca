ccc   generates two-body decay to particles of mass m1,m2
      subroutine gamdecay
      implicit double precision(a-y)
      double precision pcm(4),pboo(4),px(4),pboo1(4),pcm1(4),
     &     pboo2(4)
      integer i
      
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
      include 'decvec.f'
      include 'exp.f'
      include 'hepevt.f'

      do i=1,4
         pboo(i)=q(i,10)
      enddo

      rlalp=xd
      do i=1,3
         ralp(i)=pboo(i)/dsqrt(pboo(4)**2-mres**2)*rlalp
      enddo

      do i=1,3
         vhep(i,6)=ralp(i)*1d3
         vhep(i,7)=vhep(i,6)
      enddo

      vhep(4,6)=dsqrt(ralp(1)**2+ralp(2)**2+ralp(3)**2)/abeta*1d3
      vhep(4,7)=vhep(4,6)
      
cccccc

      do i=1,4
         pboo(i)=q(i,11)
      enddo

      egam1=pboo(4)
      
      rlgam1=(lsh+dvol-ralp(3))*pboo(4)/pboo(3)
      do i=1,3
         rgam1(i)=pboo(i)/pboo(4)*rlgam1
      enddo

cccccc

      do i=1,4
         pboo(i)=q(i,12)
      enddo

      egam2=pboo(4)
      
      rlgam2=(lsh+dvol-ralp(3))*pboo(4)/pboo(3)
      do i=1,3
         rgam2(i)=pboo(i)/pboo(4)*rlgam2
      enddo
      
      Return
      end
