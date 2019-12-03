      subroutine cscalc(p1x,p1y,p2x,p2y,wtt)
      implicit double precision(a-y)
      implicit complex*16(z)
      complex*16 pp,mm,pm,mp

      include 'vars.f'

cccccccccccccc      
 
      call axion(mx,pp,mm,pm,mp)
      pinc=cdabs(pp)**2+cdabs(mm)**2           
     &     +cdabs(pm)**2+cdabs(mp)**2
      
cccccccccccccc      

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2
  
      call formfacgam(t11,t22,x00p)
      
      dbl=dsqrt(pinc)*x00p
      wtt=dbl**2

      return
      end
