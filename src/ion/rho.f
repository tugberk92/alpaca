      function rhoz(r)
      implicit double precision(a-y)

      include 'rho0.f'
      include 'ions.f'
      
      rhoz=rho0z/(1d0+dexp((r-Rzg)/dzg))

      return
      end

