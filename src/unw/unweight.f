ccc   writes event information to array for unweighted generation
      subroutine unweight(wt,r)
      implicit double precision(a-y)
      double precision px(4),pcm(4),pboo(4)
      integer i,j

      include 'record.f'
      include 'unweighted.f'
      include 'mom.f'
      include 'hepevt.f'
      include 'leshouches.f'
      include 'wmax.f'
      include 'alpvars.f'

      
      if(wt/wmax.gt.r)then
         
         evnum=evnum+1

         do i=1,4
            px(i)=-q(i,2)
         enddo
         px(4)=q(4,2)
         ein=dsqrt(px(4)**2-px(3)**2-px(2)**2-px(1)**2)
         
         do j=3,4

            do i=1,4
               pcm(i)=q(i,j)
            enddo
            
            call boost(ein,px,pcm,pboo)
            
            do i=1,4
               evrec(evnum,j,i)=pboo(i)
            enddo
            
         enddo

         do j=5,nup+2
            
            do i=1,4
               evrec(evnum,j,i)=q(i,j+5)
            enddo
            
         enddo

         evrec(evnum,9,1)=taualp*6.58d-25
         do i=1,4
            evrec(evnum,10,i)=vhep(i,6)
         enddo
         
      endif
            
      return
      end
