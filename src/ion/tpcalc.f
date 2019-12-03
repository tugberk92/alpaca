      subroutine tpcalc
      implicit double precision(a-y)      
      integer i,j,jmax
      
      include 'tppars.f'
      include 'ions.f'
      
c      itp=1900
      itp=400

      qtmin=1d-3
      qtmax=2d0

      lqtmin=dlog(qtmin)
      lqtmax=dlog(qtmax)

      sum=0d0
      
      do i=1,itp+1

         lqt=lqtmin+(lqtmax-lqtmin)*dble(i-1)/dble(itp)
         qt=dexp(lqt)

         tparr(1,i)=lqt
         tparr(2,i)=tpz(qt)
         
      enddo

      return
      end
