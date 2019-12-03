      subroutine supinit
      implicit double precision (a-z)

      include 'mres.f'
      include 'pdg.f'
      include 'hepevt.f'
      include 'leshouches.f'
      include 'vegas.f'
      include 'gax.f'
            
      ndim=7
      
      pdgid(5)=90
      pdgid(6)=22
      pdgid(7)=22
      istup(5)=2  
      istup(6)=1  
      istup(7)=1
      mothup(1,6)=3
      mothup(2,6)=0
      mothup(1,7)=3
      mothup(2,7)=0
      icolup(1,6)=501
      icolup(2,6)=0
      icolup(1,7)=0
      icolup(2,7)=501   
      isthep(5)=2
      isthep(6)=1
      isthep(7)=1
      jdahep(1,5)=6
      jdahep(2,5)=7

      return
      end
