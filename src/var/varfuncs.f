      function dphi(i,j)
      implicit double precision(a-y)
      double precision p1(4),p2(4)
      integer i,j,k

      include 'pi.f'
      include 'mom.f'

      do k=1,4
         p1(k)=q(k,i)
         p2(k)=q(k,j)
      enddo

      dot=(p1(1)*p2(1)+p1(2)*p2(2))
     &/dsqrt((p1(1)**2+p1(2)**2)*(p2(1)**2+p2(2)**2))

      if(dot.gt.1d0)then
         dphi=0d0
      elseif(dot.lt.-1d0)then
         dphi=pi
      else
         dphi=dacos(dot)
      endif

      return
      end
      

      function phi(i)
      implicit double precision(a-y)
      double precision p1(4)
      integer i,j

      include 'mom.f'
      include 'pi.f'

      do j=1,4
         p1(j)=q(j,i)
      enddo

      if(p1(1).gt.0d0)then
         phi=datan(p1(2)/p1(1))
      elseif(p1(2).gt.0d0)then
         phi=datan(p1(2)/p1(1))+pi
      else
         phi=datan(p1(2)/p1(1))-pi
      endif

      return
      end

      function phip(p)
      implicit double precision(a-y)
      double precision p(4)
  
      et1=dsqrt(p(1)**2+p(2)**2)
      sphi=p(1)/et1
      phip=dasin(sphi)

      return
      end


      function sdot(pa,pb)
      implicit double precision(a-y)
      double precision pa(4),pb(4)   

      sdot=(pa(4)*pb(4)-pa(3)*pb(3)-pa(2)*pb(2)-pa(1)*pb(1))

      return
      end

      subroutine cdot(pa,pb,zout)
      complex*16 pa(4),zout
      double precision pb(4)

      zout=pa(4)*pb(4)-pa(3)*pb(3)-pa(2)*pb(2)-pa(1)*pb(1)
      
      return
      end

      subroutine ccdot(pa,pb,zout)
      complex*16 pa(4),pb(4),zout

      zout=pa(4)*pb(4)-pa(3)*pb(3)-pa(2)*pb(2)-pa(1)*pb(1)
      
      return
      end
