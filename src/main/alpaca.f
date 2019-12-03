ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                         c 
c     Alpaca MC for ALP production in fixed               c
c     target experiments                                  c
c                                                         c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     **************************************************
c     *   v1.1              20/05/19                   *
c     *                                                *
c     *  Author: Lucian Harland-Lang                   *
c     *  (lucian.harland-lang@physics.ox.ac.uk)        *
c     *                                                *
c     *  For details see :                             *
c     *                                                *
c     *  "A fresh look at ALP searches in fixed        *
c     *   target experiments"                          *
c     *  L.A. Harland-Lang, J. Jaeckel, M. Spannowsky  *
c     *  arXiv:1902.04878                              *
c     *                                                *
c     *  Available at :                                *
c     *  https:://alpaca.hepforge.org                  *                       
c     *                                                *
c     **************************************************
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program alpaca
      implicit double precision (a-z)
      double precision px(4),pcm(4),pboo(4)
      integer i,j,k
      integer nhistmax
      integer outl
      integer iinc,ncallu
      logical histol
      character*100 dum
      integer idum
      COMMON /ranno/ idum
      
      include 'genunw.f'
      include 'pi.f'
      include 'pdg.f'
      include 'unweighted.f'
      include 'bin.f'
      include 'vars.f'
      include 'mom.f'
      include 'norm.f'
      include 'mres.f'
      include 'record.f'
      include 'mp.f'
      include 'output.f'
      include 'nhist.f'
      include 'vegas.f'
      include 'vegaspars.f'
      include 'hepevt.f'
      include 'leshouches.f'
      include 'prec.f'
      include 'wmax.f'
      include 'wtmax.f'
      include 'gencuts.f'
      include 'ions.f'
      include 'gax.f'
      include 'exp.f'
      include 'mn.f'
      include 'adecay.f'
      include 'btype.f'
 
ccccccc

      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)ebeam
      read(*,*)btype
      read(*,*)aa
      read(*,*)az
      read(*,*)lsh
      read(*,*)dvol
      read(*,*)rmin
      read(*,*)rmax
      read(*,*)dmin
      read(*,*)emin
      read(*,*)athetamax
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)mres
      read(*,*)gax
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)outtag
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)ncall
      read(*,*)itmx
      read(*,*)prec
      read(*,*)ncall1
      read(*,*)inccall
      read(*,*)itend
      read(*,*)iseed
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)genunw
      read(*,*)nev
      read(*,*)erec
      read(*,*)dum
      read(*,*)dum
      read(*,*)gencuts
      read(*,*)adecay

cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      call length(outtag,outl)

      open(45,file='evrecs/evrec'//outtag(1:outl)//'.dat')
      wmax=0d0
      evnum=0    

      iw=0

      if(btype.eq.'prot')then
         mp=0.938272046d0
      elseif(btype.eq.'elec')then
         mp=0.511d-3
      endif
      pi=dacos(-1d0)
      conv=389379d3
      rts=dsqrt(4d0*mp**2+2d0*mp*ebeam)
  
      mn=0.938272046d0*aa
      rts=dsqrt(mp**2+mn**2+2d0*mn*(ebeam+mp))

      do i=1,20
         jdahep(1,i)=0
         jdahep(2,i)=0
      enddo

cccccccccccccccccccccccccc
      
      call supinit
      call headerout

ccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      print*,'Initialising ion form factor...' 

      call ionpars
      call rhonorm
      call rhoxycalc
      call tpcalc

      print*,'...done!'

ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      s=rts**2

      beta=dsqrt(1d0-4d0*mp**2/s)

      q(1,1)=0d0
      q(2,1)=0d0
      q(4,1)=(s+mp**2-mn**2)/2d0/rts
      q(3,1)=dsqrt(q(4,1)**2-mp**2)
      
      q(1,2)=0d0
      q(2,2)=0d0
      q(4,2)=(s+mn**2-mp**2)/2d0/rts
      q(3,2)=-dsqrt(q(4,2)**2-mn**2)

      if(btype.eq.'prot')then
         pdgid(1)=2212
      elseif(btype.eq.'elec')then
         pdgid(1)=11
      endif

      pdgid(2)=1000000000
      pdgid(2)=pdgid(2)+nint(az)*10000
      pdgid(2)=pdgid(2)+nint(aa)*10
      
      pdgid(3)=pdgid(1)
      pdgid(4)=pdgid(2)
   
ccccccccccccccccccccccccccccccccccccccccccccccc
cccc     HEPEVT
ccccccccccccccccccccccccccccccccccccccccccccccc

      do i=1,4
         px(i)=-q(i,2)
      enddo
      px(4)=q(4,2)
      ein=dsqrt(px(4)**2-px(3)**2-px(2)**2-px(1)**2)
      
      do k=1,2
         do i=1,4
            pcm(i)=q(i,k)
         enddo
         call boost(ein,px,pcm,pboo)         
         do j=1,4
            phep(j,k)=pboo(j)
         enddo
         phep(5,k)=dsqrt(q(4,k)**2-q(3,k)**2-q(2,k)**2-q(1,k)**2)
      enddo
      do k=1,20
         do j=1,4
            vhep(j,k)=0d0
         enddo
      enddo
      isthep(1)=2
      isthep(2)=2
      isthep(3)=1
      isthep(4)=1
      jmohep(2,5)=2
      jdahep(1,1)=0
      jdahep(2,1)=0
      jdahep(1,2)=0
      jdahep(2,2)=0
      jdahep(1,3)=0
      jdahep(2,3)=0
      jdahep(1,4)=0
      jdahep(2,4)=0

ccccccccccccccccccccccccccccccccccccccccccccccc
ccccc Les Houches
ccccccccccccccccccccccccccccccccccccccccccccccc

      nprup=1
      idwtup=3
      pdfsup(1)=0
      pdfsup(2)=0
      pdfgup(1)=0
      pdfgup(2)=0
      idprup=0
      xwgtup=1d0
      aqedup=alpha
      
      do k=1,2
         do i=1,4
            pcm(i)=q(i,k)
         enddo
         call boost(ein,px,pcm,pboo) 
         do j=1,4
            pup(j,k)=pboo(j)
         enddo
         pup(5,k)=dsqrt(q(4,k)**2-q(3,k)**2-q(2,k)**2-q(1,k)**2)
      enddo
      
      istup(1)=-1
      istup(2)=-1
      istup(3)=1
      istup(4)=1   
      mothup(1,1)=0
      mothup(2,1)=0
      mothup(1,2)=0
      mothup(2,2)=0
      mothup(1,3)=0
      mothup(2,3)=0
      mothup(1,4)=0
      mothup(2,4)=0
      mothup(1,5)=0
      mothup(2,5)=0
      icolup(1,1)=0
      icolup(2,1)=0
      icolup(1,2)=0
      icolup(2,2)=0
      icolup(1,3)=0
      icolup(2,3)=0
      icolup(1,4)=0
      icolup(2,4)=0
      icolup(1,5)=0
      icolup(2,5)=0
      do i=1,20
         vtimup(i)=0
         spinup(i)=9
      enddo

      do i=1,2
         do j=1,5
            jmohep(i,j)=mothup(i,j)
         enddo
      enddo

ccccccccc

      if(adecay)then
         nup=7
      else
         nup=5
      endif
      nhep=nup
      nup=nup-2

cccccccccccc
     
      nhist=0
      nhistmax=20

ccccccccccccccc

      histol=.true.

ccccccc    initialise histograms

      if(histol)call inithist(nhistmax)

cccccccccccccccc

      do i=1,10        
         xu(i)=1d0
         xl(i)=0d0
      enddo

      ACC=-1D0
      NPRN=1

      idum=-abs(iseed)
      randum=rann2()
     
      ITMX1=1
      
      bin=.false.
      unw=.false.
      calcmax=.false.
      iinc=1

      print*,''
      print*,'**************************************************************
     &**************'
      print*,'                Vegas: initialisation run                '
      print*,'**************************************************************
     &**************'

      CALL VEGAS(cs,AVGI,SD,CHI2A)

      print*,''
      print*,'**************************************************************
     &**************'
      print*,'                Vegas : main run                '
      print*,'**************************************************************
     &**************'

      ITMX=ITMX1
      NCALL=NCALL1
      avgi1=avgi
      sd1=sd

      ncall=ncall*iinc
      inccall=inccall*iinc

 779  bin=.true.

      ncallu=ncall

      unw=.false.
      calcmax=.true.
      ren=1d0

      it=1
      itmx=1

      ren=dble(ncall)

      CALL VEGAS1(cs,AVGI,SD,CHI2A)

      prec=prec*1d-2
     
 777  if(dabs(sd/avgi).gt.prec)then

         it=it+1    
         ncall=ncall+inccall     
         ren=dble(ncall)

         CALL VEGAS2(cs,AVGI,SD,CHI2A)
      

         if(it.gt.itend)goto 778

         goto 777

      endif

 778  unw=.true.
      calcmax=.true.

  10  FORMAT(' cross section = ',G16.7,' +/-',G16.7,' ( ',F9.4,' )')

cccccccccccccc

 999  if(genunw)then

         avgio=avgi
         sdo=sd
      
      print*,''
      print*,'**************************************************************
     &**************'
      print*,'Generating unweighted events'
      print*,'**************************************************************
     &**************'
      print*,''

      ncall=ncallu
      if(ncall.lt.1000)ncall=1000
 566  ren=dble(ncall)
      itmx=1
   
      CALL VEGAS2(cs,AVGI,SD,CHI2A)   

      if(evnum.lt.nev)then

      print*,''
      print*,'**************************************************************
     &**************'
      write(6,100)evnum,nev
      print*,'**************************************************************
     &**************'
      
      endif

 100  format('  generated events so far = ',i7,'    total = ',i7)

      ncall=ncall+inccall
     
      if(evnum.lt.nev)goto 566

      call unwprint

      avgi=avgio
      sd=sdo

      endif

      close(55)
      close(45)

cccccccccccccccc
       
      call length(outtag,outl)
      open(56,file='outputs/output'//outtag(1:outl)//'.dat')

      write(56,*)'**********************************************************
     &**************'
      write(56,*)
      write(56,301)avgi,sd
      write(56,*)
      write(56,*)'********************  Input parameters *******************
     &**************'
      write(56,99)' *',ebeam,' :  Beam kinetic energy (GeV)'
      write(56,96)' *',outtag(1:outl)//'.dat',' :  Output file'
      write(56,99)' *',aa,' : Target mass number'
      write(56,99)' *',az,' : Target atomic number'
      write(56,99)' *',lsh,' : Shielding length (m)'
      write(56,99)' *',dvol,' : Decay volume (m)'
      write(56,99)' *',rmin,' : Inner radius (m)'
      write(56,99)' *',rmax,' : Outer radius (m)'
      write(56,99)' *',dmin,' : Minimum photon separation (m)'
      write(56,99)' *',emin,' : Minimum photon energy (GeV)'
      write(56,99)' *',athetamax,' : Maximum ALP theta'
      write(56,*)'**********************************************************
     &**************'
      write(56,*)''
      write(56,*)'****************** Integration parameters  ***************
     &**************'
      write(56,97)' *',ncall,' :  Preconditioning calls'
      write(56,97)' *',itmx,' :  Preconditioning iterations'
      write(56,99)' *',prec*100d0,' :  Percentage accuracy'
      write(56,97)' *',ncall1,' :  Calls in first main iteration'
      write(56,97)' *',inccall,' :  Increase calls per iteration'
      write(56,97)' *',itend,' :  Maximum number of iterations'
      write(56,97)' *',iseed,' :  Random number seed'
      write(56,*)''
      write(56,*)'********************* Unweighted Events  *****************
     &**************'
      write(56,98)' *',genunw,' :  Generate unweighted events'
      write(56,99)' *',wmax,' :  Maximum weight'
      if(genunw)then
         call length(erec,outl)
         write(56,97)' *',nev,' :  Number of events'
         write(56,96)' *',erec(1:outl),' :  Record format'
      endif
      write(56,*)'**********************************************************
     &**************'
      write(56,*)''
      write(56,*)'*****************************  Cuts  *********************
     &**************'
      write(56,98)' *',gencuts,' :  Generate cuts'
      write(56,98)' *',adecay,' : Include ALP decay'
      write(56,*)'**********************************************************
     &**************'
      write(56,*)''

 99   format(a,f24.4,8x,a)
 97   format(a,i24,8x,a)
 96   format(a,a24,8x,a)
 98   format(a,l24,8x,a)

      close(56)

      if(histol)then

      do j=1,nhist
           call histo2(j,0)
      enddo

      endif

      print*,''
      write(6,301)avgi,sd
      print*,''

 301  format(' Cross section = ',G16.7,' +/-',G16.7,' pb')

      call cpu_time(t2)
      print*,'time elapsed = ', t2, ' s'
      
      stop
      end
