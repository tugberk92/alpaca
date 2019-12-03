      subroutine rhonorm
      implicit double precision(a-y)
      integer i,itot

      include 'rho0.f'
      include 'ions.f'
      include 'pi.f'

      rho0z=1d0
      
      rmax=3d0*rzg
      itot=500
      hr=rmax/dble(itot)

      sumz=0d0

      do i=1,itot

         r=(dble(i)-0.5d0)*hr

         wtz=rhoz(r)
         wtz=wtz*4d0*pi
         wtz=wtz*r**2*hr

         sumz=sumz+wtz
         
      enddo

      rho0z=az/sumz

      return
      end
