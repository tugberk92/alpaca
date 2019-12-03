      subroutine cut(icut)
      implicit double precision(a-y)
      double precision p1(4),p2(4),p3(4)
      integer icut,jflag,i
      logical accut

      include 'gencuts.f'
      include 'vars.f'
      include 'mom.f'
      include 'pi.f'
      include 'alpvars.f'
      include 'exp.f'
      include 'decvec.f'
      include 'adecay.f'
      
      icut=0
      
      rg1=dsqrt((ralp(1)+rgam1(1))**2+(ralp(2)+rgam1(2))**2)
      rg2=dsqrt((ralp(1)+rgam2(1))**2+(ralp(2)+rgam2(2))**2)
      dgg=dsqrt((rgam1(1)-rgam2(1))**2+(rgam1(2)-rgam2(2))**2)
 
      if(adecay)then
      
         if(egam1.lt.emin/2d0)return
         if(egam2.lt.emin/2d0)return
         if(rg1.gt.rmax)return
         if(rg2.gt.rmax)return
         if(rg1.lt.rmin)return
         if(rg2.lt.rmin)return
         if(dgg.lt.dmin)return

      else

         if(atheta.gt.athetamax)return

      endif
      
      icut=1

      return
      end
 
