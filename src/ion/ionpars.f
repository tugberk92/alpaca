      subroutine ionpars
      implicit double precision(a-y)
      integer fit
      
      include 'ions.f'
      include 'pi.f'
    
      ss0=0.9d0
      ss=0.9d0/0.1973d0
      
      rz=dsqrt((1.23d0*aa**(1d0/3d0)-0.6d0)**2+7d0*pi**2/3d0*0.52d0**2-
     &     5d0*ss0**2)  
       
      dz=0.55d0

      rzg=rz/0.1973d0
      dzg=dz/0.1973d0
 
      return
      end
