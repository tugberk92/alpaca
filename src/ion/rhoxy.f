      function rhozxy(rxy)
      implicit double precision(a-y)
      integer i,itot

      include 'ions.f'
      
      itot=10000

      rzmax=10d0*rzg
      hz=rzmax/dble(itot)
      
      sum=0d0
      
      do i=1,itot

         rrz=(dble(i)-0.5d0)*hz

         r=dsqrt(rrz**2+rxy**2)
         wt=rhoz(r)*hz
         
         sum=sum+wt

      enddo

      rhozxy=sum*2d0

      return
      end
