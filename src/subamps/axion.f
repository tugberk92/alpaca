ccc   gamgam --> ALP subprocess amplitude
      subroutine axion(mx,pp,mm,pm,mp)
      implicit double precision (a-z)
      complex*16 pp,mm,pm,mp

      include 'norm.f'
      include 'gax.f'
      
      norma=gax/2d0*mx**2
      norma=norma*dsqrt(conv)
      norma=norma/2d0   ! to get correct M^2/4
      
      pp=norma
      mm=-norma
      pm=0d0
      mp=0d0
      
      return
      end







