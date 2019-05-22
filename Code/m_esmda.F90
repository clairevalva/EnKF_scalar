module m_esmda
use m_cyyreg
use m_cov
use m_getcaseid
use m_func
use m_getalpha
use m_normal
use m_tecmargpdf
use m_tecpdf
use mod_xyqgrid
implicit none
contains
subroutine esmda(samples,xsampini,qsampini,nrsamp,esamp,&
                 beta,funcmode,nmda,alphageo,cdd,d,lcyyreg,sigw,sigq)
   integer, intent(in)  :: nrsamp
   integer, intent(in)  :: esamp
   real,    intent(out) :: samples(nrsamp,2)
   real,    intent(in)  :: xsampini(nrsamp) 
   real,    intent(in)  :: qsampini(nrsamp) 

   real, intent(in) :: beta
   real, intent(in) :: alphageo
   real, intent(in) :: d
   real, intent(in) :: cdd
   real, intent(in) :: sigw,sigq

   logical, intent(in) :: lcyyreg
   integer, intent(in) :: funcmode
   integer, intent(in) :: nmda

   real, allocatable :: xsamp(:)
   real, allocatable :: qsamp(:)
   real, allocatable :: ysamp(:)

   real, allocatable :: alpha(:)
   integer n,i
   real Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,alphasum
   integer :: gradient=0
   real pert
   character(len=40) caseid

   allocate(xsamp(nrsamp))
   allocate(qsamp(nrsamp))
   allocate(ysamp(nrsamp))

   allocate(alpha(nmda))

   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')'MDA analysis...'
   do i=1,nrsamp
      xsamp(i)=xsampini(i)
      qsamp(i)=qsampini(i)
      ysamp(i)=func(xsamp(i),beta,funcmode)+qsamp(i)
   enddo

      alphasum=0.0
      do n=1,nmda
         alpha(n)=getalpha(n,nmda,alphageo)
         alphasum=alphasum+1.0/alpha(n)
         print *,'alpha            :',alpha(n)
         call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
         write(*,'(a,f10.4)')'Cyy from samples :',Cyy
         if (lcyyreg) cyy=cyyreg(Cxx,Cqq,Cyx,Cqy,Cqx)
         write(*,'(i3,f10.2,a,2f13.5,e13.5)')n,alpha(n),', cxx= ',cxx,cyx,cqy

         do i=1,nrsamp
            pert=sqrt(alpha(n))*sqrt(cdd)*normal()
            xsamp(i)=xsamp(i) + (Cyx/(Cyy+alpha(n)*Cdd))*(d+pert-ysamp(i))
            qsamp(i)=qsamp(i) + (Cqy/(Cyy+alpha(n)*Cdd))*(d+pert-ysamp(i))
            ysamp(i)=func(xsamp(i),beta,funcmode)+qsamp(i)
         enddo
         call getcaseid(caseid,'MDA',alphageo,nmda,esamp,gradient,beta,sigw,n)
         call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
      enddo

      samples(:,1)=xsamp(:)
      samples(:,2)=qsamp(:)

      if (sigw < sigq) then
         do i=1,nrsamp
            ysamp(i)=ysamp(i)+sigq*normal()
         enddo
      endif
      call getcaseid(caseid,'MDA',alphageo,nmda,esamp,gradient,beta,sigw,0)
      call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
      call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
      call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
      call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
      write(*,'(a,f12.4)')'Alphasum=',alphasum
      write(*,'(a)')'ES-MDA analysis completed'
end subroutine
end module
