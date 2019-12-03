      subroutine rhoxycalc
      implicit double precision(a-y)       
      implicit complex*16(z)
      integer i
      
      include 'pi.f'
      include 'ions.f'
      include 'rhoxypars.f'
      include 'rho0.f'
      
ccccccccccc

      sum=0d0
      
      irho=900

      rzmax=3d0*rzg
      rzmin=0d0
      
      do i=1,irho+1

         rz=rzmin+(rzmax-rzmin)*dble(i-1)/dble(irho)
         
         rhoxyarr(1,i)=rz
         rhoxyarr(2,i)=rhozxy(dabs(rz))
   
      enddo


      return
      end
