module m_iniens
use mod_xyqgrid
use m_getcaseid
use m_normal
use m_tecpdf
use m_func
use m_tecmargpdf
use m_moments
implicit none
contains
subroutine iniens(xsampini,qsampini,ysampini,dpert,nrsamp,              &
                  siga,sigw,x0,beta,funcmode,alphageo,nmda,esamp,d,cdd)
   integer, intent(in)   :: nrsamp,funcmode,nmda,esamp
   real,    intent(out)  :: xsampini(nrsamp),qsampini(nrsamp),ysampini(nrsamp),dpert(nrsamp)
   real,    intent(in)   :: beta,x0,alphageo,siga,sigw,d,cdd

   real avex,adevx,sdevx,varx,skewx,kurtx
   integer :: gradient=0
   character(len=40) caseid
   integer i

   do i=1,nrsamp
      xsampini(i)=siga*normal()
   enddo
   call moments(xsampini,nrsamp,avex,adevx,sdevx,varx,skewx,kurtx)
   do i=1,nrsamp
      xsampini(i)=xsampini(i)-avex
   enddo

   do i=1,nrsamp
      qsampini(i)=sigw*normal()
   enddo
   call moments(qsampini,nrsamp,avex,adevx,sdevx,varx,skewx,kurtx)
   do i=1,nrsamp
      qsampini(i)=qsampini(i)-avex
   enddo

   do i=1,nrsamp
      xsampini(i)=x0+xsampini(i)
      ysampini(i)=func(xsampini(i),beta,funcmode) + qsampini(i)
   enddo

   do i=1,nrsamp
      dpert(i)=d+sqrt(cdd)*normal()
   enddo
   print *,'Sampling done'

   call getcaseid(caseid,'INI',alphageo,nmda,esamp,gradient,beta,sigw,0)
   call tecpdf(x,y,nx,ny,xsampini,ysampini,nrsamp,xa,ya,dx,dy,caseid)
   call tecmargpdf('x',xsampini,nrsamp,caseid,xa,xb,nx)
   call tecmargpdf('y',ysampini,nrsamp,caseid,ya,yb,ny)
   call tecmargpdf('q',qsampini,nrsamp,caseid,qa,qb,nx)

end subroutine
end module
