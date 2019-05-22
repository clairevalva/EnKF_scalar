module m_es
use mod_xyqgrid
use m_cyyreg
use m_cov
use m_getcaseid
use m_func
use m_normal
use m_tecmargpdf
use m_tecpdf
implicit none
contains
subroutine es(samples,xsampini,qsampini,dpert,nrsamp,esamp,              &
              beta,funcmode,nmda,alphageo,cdd,d,lcyyreg,sigw,sigq)
   integer, intent(in)  :: nrsamp
   integer, intent(in)  :: esamp
   integer, intent(in)  :: nmda
   real,    intent(out) :: samples(nrsamp,2)
   real,    intent(in)  :: xsampini(nrsamp) 
   real,    intent(in)  :: qsampini(nrsamp) 
   real,    intent(in)  :: dpert(nrsamp) 

   real, intent(in) :: alphageo
   real, intent(in) :: beta
   real, intent(in) :: d
   real, intent(in) :: cdd
   real, intent(in) :: sigw,sigq

   logical, intent(in) :: lcyyreg
   integer, intent(in) :: funcmode

   real, allocatable :: xsamp(:)
   real, allocatable :: qsamp(:)
   real, allocatable :: ysamp(:)

   integer n,i
   real Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,alphasum
   integer :: gradient=0
   real pert
   character(len=40) caseid

   allocate(xsamp(nrsamp))
   allocate(qsamp(nrsamp))
   allocate(ysamp(nrsamp))

   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')'ES analysis...'
   do i=1,nrsamp
      xsamp(i)=xsampini(i)
      qsamp(i)=qsampini(i)
      ysamp(i)=func(xsamp(i),beta,funcmode)+qsamp(i)
   enddo

   call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)

!   write(*,'(a,f10.4)')'Pyy from samples :',Cyy
   if (lcyyreg) cyy=cyyreg(Cxx,Cqq,Cyx,Cqy,Cqx)

   do i=1,nrsamp
      xsamp(i)=xsamp(i) + (cyx/(cyy+cdd))*(dpert(i)-ysamp(i))
      qsamp(i)=qsamp(i) + (cqy/(cyy+cdd))*(dpert(i)-ysamp(i))

!      if (updatemode==0) then
!         ysamp(i)=ysamp(i)+ (cyy/(cyy+cdd))*(dpert(i)-ysamp(i))
!      else
         ysamp(i)=func(xsamp(i),beta,funcmode) + qsamp(i)
!      endif

   enddo
   write(*,'(a,10g11.3)')'Es',xsamp(1:10)

   if (sigw < sigq) then
      do i=1,nrsamp
         ysamp(i)=ysamp(i)+sigq*normal()
      enddo
   endif
   call getcaseid(caseid,'ES',alphageo,nmda,esamp,gradient,beta,sigw,0)
   call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
   call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
   call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
   call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
   write(*,'(a)')'ES analysis completed'

   samples(:,1)=xsamp(:)
   samples(:,2)=qsamp(:)
end subroutine
end module
