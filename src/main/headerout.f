ccc   prints out header information
      subroutine headerout
      implicit double precision(a-y)
      integer outl
    
      include 'unweighted.f'
      include 'vegas.f'
      include 'vegaspars.f'
      include 'genunw.f'
      include 'vars.f'  
      include 'record.f'
      include 'output.f'
      include 'prec.f'
      include 'ions.f'
      include 'gax.f'
      include 'exp.f'
      include 'gencuts.f'
      include 'adecay.f'

      call length(outtag,outl)

      print*,'*************************************************'
      print*,'***************  Alpaca v1.1  *******************'
      print*,'*************************************************'
      print*,'*   v1.01              DATE  20/05/19           *'
      print*,'*                                               *'
      print*,'*  Author: Lucian Harland-Lang                  *'
      print*,'*  (lucian.harland-lang@physics.ox.ac.uk        *'
      print*,'*                                               *'
      print*,'*  For details see                              *'
      print*,'*                                               *'
      print*,'*  "A fresh look at ALP searches in fixed       *'
      print*,'*   target experiments"                         *'
      print*,'*  L.A. Harland-Lang, J. Jaeckel, M. Spannowsky *'
      print*,'*  arXiv 1902.04878                             *'
      print*,'*                                               *'
      print*,'*  Available at :                               *'
      print*,'*  https://alpaca.hepforge.org                  *'
      print*,'*                                               *'
      print*,'*************************************************'
      print*,''
      print*,'********************  Input parameters *******************
     &**************'
      write(*,99)' *',ebeam,' :  Beam kinetic energy (GeV)'
      write(*,96)' *',outtag(1:outl)//'.dat',' :  Output file'
      write(*,99)' *',aa,' : Target mass number'
      write(*,99)' *',az,' : Target atomic number'
      write(*,99)' *',lsh,' : Shielding length (m)'
      write(*,99)' *',dvol,' : Decay volume (m)'
      write(*,99)' *',rmin,' : Inner radius (m)'
      write(*,99)' *',rmax,' : Outer radius (m)'
      write(*,99)' *',dmin,' : Minimum photon separation (m)'
      write(*,99)' *',emin,' : Minimum photon energy (GeV)'
      write(*,99)' *',athetamax,' : Maximum ALP theta'
      print*,'**********************************************************
     &**************'
      print*,''
      print*,'********************* Integration parameters  ******************
     &**************'
      write(*,97)' *',ncall,' :  Preconditioning calls'
      write(*,97)' *',itmx,' :  Preconditioning iterations'
      write(*,99)' *',prec,' :  Percentage accuracy'
      write(*,97)' *',ncall1,' :  Calls in first main iteration'
      write(*,97)' *',inccall,' :  Increase calls per iteration'
      write(*,97)' *',itend,' :  Maximum number of iterations'
      write(*,97)' *',iseed,' :  Random number seed'
      print*,'**********************************************************
     &**************'
      print*,''
      print*,'********************* Unweighted Events  *****************
     &**************'
      write(*,98)' *',genunw,' :  Generate unweighted events'
      if(genunw)then
         call length(erec,outl)
         write(*,97)' *',nev,' :  Number of events'
         write(*,96)' *',erec(1:outl),' : Record format'
      endif
      print*,'**********************************************************
     &**************'
      print*,''
      print*,'*****************************  Cuts  *********************
     &**************'
      write(*,98)' *',gencuts,' :  Generate cuts'
      write(*,98)' *',adecay,' : Include ALP decay'
      print*,'**********************************************************
     &**************'
      print*,''
      

 99   format(a,f24.4,8x,a)
 97   format(a,i24,8x,a)
 96   format(a,a24,8x,a)
 98   format(a,l24,8x,a)

      return
      end
