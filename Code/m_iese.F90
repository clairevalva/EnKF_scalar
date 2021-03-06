module m_iese
use mod_inistat
use mod_xyqgrid
use m_cyyreg
use m_cov
use m_getcaseid
use m_func
use m_iescostf
use m_normal
use m_tecmargpdf
use m_tecsampini
implicit none
logical, save :: liese                ! Run IESE or not
integer, save :: maxieseit            ! maximumn number of iterations
real,    save :: gamma_iese           ! step length used in IESE
real,    save :: iese_eps=0.0000001

contains 
   subroutine iese(samples,xf,qf,nrsamp,datum)
   implicit none
   integer, intent(in)  :: nrsamp            ! Number of samples
   real,    intent(out) :: samples(nrsamp,2) ! Returns posterior samples
   real,    intent(in)  :: xf(nrsamp)        ! Prior samples of x
   real,    intent(in)  :: qf(nrsamp)        ! Prior samples of q
   real,    intent(in)  :: datum             ! The measurement

   real, allocatable :: xsamp(:)
   real, allocatable :: qsamp(:)
   real, allocatable :: ysamp(:)
   integer, allocatable :: iconv(:)

   integer n,i,sumconv

   real Cxx,Cyy,Cqq,Cyx,Cqy,Cqx
   real Pxx,Pyy,Pqq,Pyx,Pqy,Pqx

   real Czz(2,2),CIzz(2,2)
   real Pzz(2,2),PIzz(2,2),Pzy(2,1),Pyz(1,2)

   real grad2(2)
   real H(2,2),HI(2,2)

   real zi(2),zf(2)

   real G(1,2),GT(2)
   real Ie(2,1)

   if (Cqq == 0.0) then
      print '(a)','SIESE only works for model error > zero'
      stop
   endif


   allocate(xsamp(nrsamp))
   allocate(qsamp(nrsamp))
   allocate(ysamp(nrsamp))
   allocate(iconv(nrsamp)) 

   Ie=0.0
   Ie(2,1)=1.0

   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)',advance='yes')'IESE analysis'

   do n=1,nrsamp
      iconv(n)=0
      xsamp(n)=xf(n)
      qsamp(n)=qf(n)
      ysamp(n)=func(xsamp(n),qsamp(n))
   enddo
   call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
   Czz(1,1)=Cxx;  Czz(2,2)=Cqq;  Czz(1,2)=Cqx;   Czz(2,1)=Cqx
   CIzz(1,1)=Cqq; CIzz(2,2)=Cxx; CIzz(1,2)=-Cqx; CIzz(2,1)=-Cqx
   if (Cqq > 0.0) CIzz=CIzz/(Cxx*Cqq-Cqx*Cqx)

   sumconv=0
   do i=1,maxieseit
      if (mod(i,1) == 0) then
         write(*,'(a,i4,a)',advance='no')'i=',i,'...'
         print *,'converged realizations=', sumconv,nrsamp,real(100*sumconv)/real(nrsamp)
      endif
      call cov(Pxx,Pyy,Pqq,Pyx,Pqy,Pqx,xsamp,ysamp,qsamp,nrsamp)
      Pzz(1,1)=Pxx;  Pzz(2,2)=Pqq;  Pzz(1,2)=Pqx;   Pzz(2,1)=Pqx
      PIzz(1,1)=Pqq; PIzz(2,2)=Pxx; PIzz(1,2)=-Pqx; PIzz(2,1)=-Pqx
      if (Pqq > 0.0) PIzz=PIzz/(Pxx*Pqq-Pqx*Pqx)

      Pzy(1,1)=Pyx; Pzy(2,1)=Pqy
      Pyz(1,1)=Pyx; Pyz(1,2)=Pqy

!      if (lcyyreg) Pyy=cyyreg(Pxx,Pqq,Pyx,Pqy,Pqx)

      G=matmul(Pyz,PIzz) - transpose(Ie)
      GT(1)=G(1,1); GT(2)=G(1,2)

      do n=1,nrsamp
         if (iconv(n) > 0) cycle ! do nothing for converged realizations

         zf(1)=xf(n)
         zf(2)=qf(n)

         zi(1)=xsamp(n)
         zi(2)=qsamp(n)

         H=CIzz + matmul( transpose(G) , G )/Cdd
         HI(1,1)=H(2,2); HI(2,2)=H(1,1); HI(1,2)=-H(2,1); HI(2,1)=-H(1,2)
         if (H(2,2) > 0.0) HI=HI/(H(1,1)*H(2,2)-H(1,2)*H(2,1))

         grad2 = matmul(CIzz,zi-zf) +GT*( ysamp(n) - datum - qsamp(n)) / Cdd

         zi(:) = zi(:) - gamma_iese*matmul(HI,grad2)

         xsamp(n)=zi(1)
         qsamp(n)=zi(2)
         ysamp(n)=func(xsamp(n),qsamp(n))

         if ((abs(grad2(1)) < iese_eps) .and. (abs(grad2(2)) < iese_eps)) iconv(n)=1
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sumconv=sum(iconv(1:nrsamp))
      if (sumconv == nrsamp) then
         write(*,'(a,i4,a)',advance='no')'i=',i,'...'
         print *,'converged realizations=', sumconv,nrsamp,real(100*sumconv)/real(nrsamp)
         print *,'Exiting IESE iterations'
         exit
      endif

   enddo
   samples(:,1)=xsamp(:)
   samples(:,2)=qsamp(:)


   deallocate(xsamp,ysamp,qsamp,iconv)
   write(*,'(a)')'IESE analysis completed'
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')

end subroutine
end module
