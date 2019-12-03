ccc   prints out header information
      subroutine headerlhe
      implicit double precision(a-y)
      integer outl
    
      include 'gencuts.f'
      include 'unweighted.f'
      include 'vegas.f'
      include 'vegaspars.f'
      include 'genunw.f'
      include 'vars.f'
      include 'pdg.f'   
      include 'record.f'
      include 'output.f'
      include 'ions.f'
      include 'gax.f'
      include 'exp.f'
      include 'adecay.f'
      
      call length(procn,outl)


      write(45,*)'*************************************************'
      write(45,*)'***************  Alpaca v1.00  ******************'
      write(45,*)'*************************************************'
      write(45,*)'*   v1.00              DATE  14/02/19           *'
      write(45,*)'*                                               *'
      write(45,*)'*  Author: Lucian Harland-Lang                  *'
      write(45,*)'*  (lucian.harland-lang@physics.ox.ac.uk        *'
      write(45,*)'*                                               *'
      write(45,*)'*  For details see                              *'
      write(45,*)'*                                               *'
      write(45,*)'*  "A fresh look at ALP searches in fixed       *'
      write(45,*)'*   target experiments"                         *'
      write(45,*)'*  L.A. Harland-Lang, J. Jaeckel, M. Spannowsky *'
      write(45,*)'*  arXiv 1902.04878                             *'
      write(45,*)'*                                               *'
      write(45,*)'*  Available at :                               *'
      write(45,*)'*  https://alpaca.hepforge.org                  *'
      write(45,*)'*                                               *'
      write(45,*)'*************************************************'
      write(45,*)''
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      write(45,*)'********************  Input parameters ***************
     &******************'
      write(45,99)' *',rts,' :  CMS collision energy (GeV)'
      call length(outtag,outl)
      write(45,96)' *',outtag(1:outl)//'.dat',' :  Output file'
      write(45,99)' *',aa,' : Target mass number'
      write(45,99)' *',az,' : Target atomic number'
      write(45,99)' *',lsh,' : Shielding length (m)'
      write(45,99)' *',dvol,' : Decay volume (m)'
      write(45,99)' *',rmin,' : Inner radius (m)'
      write(45,99)' *',rmax,' : Outer radius (m)'
      write(45,99)' *',dmin,' : Minimum photon separation (m)'
      write(45,99)' *',emin,' : Minimum photon energy (GeV)'
      write(45,99)' *',athetamax,' : Maximum ALP theta'
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      write(45,*)'********************* Integration parameters  ********
     &******************'
      write(45,97)' *',ncall,' :  Preconditioning calls'
      write(45,97)' *',itmx,' :  Preconditioning iterations'
      write(45,99)' *',prec,' :  Percentage accuracy'
      write(45,97)' *',ncall1,' :  Calls in first main iteration'
      write(45,97)' *',inccall,' :  Increase calls per iteration'
      write(45,97)' *',itend,' :  Maximum number of iterations'
      write(45,97)' *',iseed,' :  Random number seed'
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      write(45,*)'********************* Unweighted Events  *************
     &******************'
      write(45,98)' *',genunw,' :  Generate unweighted events'
      if(genunw)then
         call length(erec,outl)
         write(45,97)' *',nev,' :  Number of events'
         write(45,96)' *',erec(1:outl),' : Record format'
      endif
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      write(45,*)'************************ Cuts ************************
     &******************'
      write(45,98)' *',gencuts,' :  Generate cuts'
      write(45,98)' *',adecay,' : Include ALP decay'
      write(45,*)'******************************************************
     &******************'
      write(45,*)''

 99   format(a,f24.4,8x,a)
 97   format(a,i24,8x,a)
 96   format(a,a24,8x,a)
 98   format(a,l24,8x,a)

      return
      end
