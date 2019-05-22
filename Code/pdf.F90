program iterative_smoothers
   use mod_inistat
   use mod_xyqgrid
   use m_marginalpdf
   use m_set_random_seed2
   use m_moments
   use m_avevar
   use m_func
   use m_steinkern
   use m_cyyreg
   use m_aaprojection
   use m_iniens
   use m_es
   use m_ies
   use m_sies
   use m_esmda
   use m_costf
   use m_tecpdf
   use m_omegafact
   use m_tecmargpdf
   use m_pseudoinv
   use m_pseudoinv2
   use m_printcostf
   use m_tecjointpdf
   use m_tecfunc
   use m_teccostens
   use m_tecsampini
   use m_normal
   use m_cov
   use m_integrals
   implicit none

   integer nrsamp,esamp



   logical :: lenstein=.true.  ! Run EnStein or not

   character(len=7) :: method(0:4)=['ES     ','IES    ','SIES   ', 'ESMDA  ', 'EnSTEIN']
   character(len=1) :: variable(1:2)=['x','q']

   real, allocatable :: kern(:,:)

   real ave,var
   real sump,fac,tmp
   real grad1,grad2(2)
   real cxx,cyy,cqq,cyx,cqy,cqx
   real pxx,Pyy,pqq,pyx,pqy,pqx,Pyymat(1,1)
   real cdd,G
   real WoodA,WoodB
   real hessian, cvec(2) , psca
   real Czz(2,2),Pzz(2,2),PIzz(2,2),CIzz(2,2),zi(2),zf(2),Pzy(2,1),Pyz(1,2)
   real n1


   real, allocatable :: sf(:),si(:,:),grad(:,:),newgrad(:,:)
   real :: xlength=1.0


   logical lmoderr, testpseudo


   real, allocatable :: xsamp(:),xsampini(:) 
   real, allocatable :: ysamp(:),ysampini(:)
   real, allocatable :: qsamp(:),qsampini(:)
   real, allocatable :: dpert(:)
   integer, allocatable :: iconv(:)
   real pert

   character(len=3) tag3
   character(len=2) tag2
   character(len=1) tag1
   character(len=40) caseid

   integer i,j,k,m,n,loc(1),ival,jval
   integer sumconv

   real jx,djx,djq,ddjx
   real dgx,dgq
   real aveA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IEnKS variables

   real, allocatable :: samples(:,:,:) ! for comparing ES, IES and SIES
   integer innovation
   logical :: lpseudoinv=.false.
   integer info
   logical :: rank1svd=.false.
   real start,finish
   integer :: ndim
   integer, parameter :: nrobs=1
   character(len=2) ca
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call set_random_seed2()
   print '(a,g13.5)','Normal number', normal()

   open(10,file='infile.in')
      read(10,*)esamp       ; print '(a,i4,i12)',   'number  of samples 10^x    :',esamp,10**esamp
      read(10,'(a)')ca      
      if (ca /= '#1') then
         print *,'#1: error in infile.in'
         stop
      endif
      read(10,*)xa          ; print '(a,f10.3)',    'xa                         :',xa
      read(10,*)xb          ; print '(a,f10.3)',    'xb                         :',xb
      read(10,*)qa          ; print '(a,f10.3)',    'qa                         :',qa
      read(10,*)qb          ; print '(a,f10.3)',    'qb                         :',qb
      read(10,*)ya          ; print '(a,f10.3)',    'ya                         :',ya
      read(10,*)yb          ; print '(a,f10.3)',    'yb                         :',yb
      read(10,'(a)')ca
      if (ca /= '#2') then
         print *,'#2: error in infile.in'
         stop
      endif
      read(10,*)x0          ; print '(a,f10.3)',    'prior estimate of x        :',x0
      read(10,*)siga        ; print '(a,f10.3)',    'err std dev of prior       :',siga
      read(10,*)d           ; print '(a,f10.3)',    'measurement of y           :',d
      read(10,*)sigo        ; print '(a,f10.3)',    'measurement std dev        :',sigo
      read(10,*)sigw        ; print '(a,f10.3)',    'model std dev              :',sigw
      read(10,'(a)')ca
      if (ca /= '#3') then
         print *,'#3: error in infile.in'
         stop
      endif
      read(10,*)beta        ; print '(a,f10.3)',    'model parameter beta       :',beta
      read(10,*)funcmode    ; print '(a,tr7,i3)',   'function to use            :',funcmode
      read(10,'(a)')ca
      if (ca /= '#4') then
         print *,'#4: error in infile.in'
         stop
      endif
      read(10,*)lcyyreg     ; print '(a,tr10,l1)',  'Regression for Cyy         :',lcyyreg
      read(10,'(a)')ca
      if (ca /= '#5') then
         print *,'#5: error in infile.in'
         stop
      endif
      read(10,*)lmda        ; print '(a,tr10,l1)',  'Run MDA                    :',lmda
      read(10,*)nmda        ; print '(a,tr7,i3)',   'number of mda iterations   :',nmda
      read(10,*)alphageo    ; print '(a,f10.3)',    'geometrical alpha value    :',alphageo
      read(10,'(a)')ca
      if (ca /= '#6') then
         print *,'#6: error in infile.in'
         stop
      endif
      read(10,*)lies        ; print '(a,tr10,l1)',  'Run IES                    :',lies
      read(10,*)maxiesit    ; print '(a,i5)',       'maxiesit                   :',maxiesit
      read(10,*)gamma_ies   ; print '(a,f10.3)',    'gamma_ies                  :',gamma_ies
      read(10,*)IESv        ; print '(a,i1)',       'IESv                       :',IESv
      read(10,'(a)')ca
      if (ca /= '#7') then
         print *,'#7: error in infile.in'
         stop
      endif
      read(10,*)lsies      ; print '(a,tr10,l1)',  'Run SIES                   :',lsies
      read(10,*)maxsiesit   ; print '(a,i5)',       'maxiesit                   :',maxsiesit
      read(10,*)gamma_sies  ; print '(a,f10.3)',    'gamma_sies                 :',gamma_sies
      read(10,*)ndim        ; print '(a,i5)',       'ndim                       :',ndim
   close(10)

   if (xa==xb) then
      xa=-5.0*siga
      xb= 5.0*siga
      print '(a,2f10.4)','xa and xb set to :',xa,xb
   endif

   if (qa==qb) then
      qa=-5.0*max(sigq,sigw)
      qb= 5.0*max(sigq,sigw)
      print '(a,2f10.4)','qa and qb set to :',qa,qb
   endif

!   if (funcmode == 0 .and. (nx /= 450 .or. ny /= 450)) stop 'check mod_dimensions'
!   if (funcmode == 2 .and. (nx /= 450 .or. ny /= 250)) stop 'check mod_dimensions'
         

!Include model errors in inverse calculation (Tarantola approach)
   cdd=sigo**2
   nrsamp=10**esamp
   nrsamp=1*nrsamp
   allocate(samples(nrsamp,2,0:4))
   allocate(dpert(nrsamp))
   allocate(qsamp(nrsamp),qsampini(nrsamp))
   allocate(xsamp(nrsamp),xsampini(nrsamp)) 
   allocate(ysamp(nrsamp),ysampini(nrsamp))
   allocate(iconv(nrsamp)) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Theroretical values for statistical moements (trustat.dat)
   call integrals(maxiesit)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid domain for plotting pdfs in x, q, and y
   call xyqgrid(x,y,q)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of the strong constraint cost function
   do i=1,nx
      prior(i)=exp(-0.5*(x(i)-x0)**2/siga**2)
      cost(i)=(x(i)-x0)**2/siga**2 + (func(x(i))-d)**2/cdd    
   enddo
   sump=sum(prior(:))*dx
   prior=prior/sump
   call tecfunc('costf',cost,x,nx,'x','Cost')
   call tecfunc('prior',prior,x,nx,'x','Prior')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Datum pdf
   do j=1,ny
      datum(j)=exp(-0.5*(y(j)-d)**2/cdd)
   enddo
   sump=sum(datum(:))*dy
   datum=datum/sump
   call tecfunc('datum',datum,y,ny,'y','Datum')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of the analytical joint unconditional pdf
   do j=1,ny
   do i=1,nx
      pdf(i,j)=exp( -0.5*(x(i)-x0)**2/siga**2            &
                    -0.5*(y(j)-func(x(i)))**2/max(sigw,sigq)**2 )
   enddo
   enddo
   sump=sum(pdf(:,:))*dx*dy
   pdf=pdf/sump
   call getcaseid(caseid,'PDFJ',-1.0,-1,esamp,sigw,0)
   call tecjointpdf(pdf,x,y,nx,ny,caseid)
   call marginalpdf(pdf,margx,margy,nx,ny,x,y,dx,dy)
   call tecfunc('margx',margx,x,nx,'x',caseid)
   call tecfunc('margy',margy,y,ny,'y',caseid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of the analytical joint conditional pdf
   do j=1,ny
   do i=1,nx
      pdf(i,j)=exp(-0.5*(x(i)-x0)**2/siga**2            &
                   -0.5*(y(j)-func(x(i)))**2/max(sigw,sigq)**2  &
                   -0.5*(y(j)-d)**2/cdd)    
   enddo
   enddo
   sump=sum(pdf(:,:))*dx*dy
   pdf=pdf/sump
   call getcaseid(caseid,'PDFC',-1.0,1,esamp,sigw,0)
   call tecjointpdf(pdf,x,y,nx,ny,caseid)
   call marginalpdf(pdf,margx,margy,nx,ny,x,y,dx,dy)

   margx(:)=0.0
   do i=1,nx
   do j=1,ny
      if (sigw==0.0) then
         margx(i)=margx(i)+exp(-0.5*(x(i)-x0 )**2/siga**2                   &
                               -0.5*(d-func(x(i)))**2/cdd**2)
      else 
         margx(i)=margx(i)+exp(-0.5*(x(i)-x0 )**2/siga**2                   &
                               -0.5*(q(j)-0.0)**2/max(sigw,0.0001)**2  & 
                               -0.5*(d-func(x(i))-q(j))**2/cdd**2)
      endif
   enddo
   enddo
   margx=margx/(sum(margx(:))*dx)
   call tecfunc('margx',margx,x,nx,'x',caseid)
   call tecfunc('margy',margy,y,ny,'y',caseid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ensemble initialization
   call iniens(xsampini,qsampini,ysampini,dpert,nrsamp,alphageo,nmda,esamp,cdd)

! ES 
   call es(samples(1:nrsamp,1:2,0),xsampini,qsampini,dpert,nrsamp,esamp,cdd)

   call printcostf(xsampini,xsamp,dpert,nrsamp,cdd)

! IES
   if (lies) then
      call ies(samples(1:nrsamp,1:2,1),xsampini,qsampini,nrsamp,esamp,cdd,dpert)
   endif

! Subspace IES
   if (lsies) then
      call sies(samples(1:nrsamp,1:2,2),xsampini,qsampini,nrsamp,esamp,cdd,dpert,nrobs,ndim)
   endif

! ESMDA
   if (lmda) then
      call esmda(samples(1:nrsamp,1:2,3),xsampini,qsampini,nrsamp,esamp,cdd)
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! E n S T E I N !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (lenstein) then
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
!    EnStein update
      write(*,'(a)',advance='yes')'EnStein analysis...'
      allocate(sf(2),si(2,nrsamp),grad(2,nrsamp),newgrad(2,nrsamp))
!      allocate(kern(nrsamp,nrsamp))

         xlength = 2.0/3.14159265 ; print '(a,f13.5)','xlength=',xlength
!        xlength = 1.0            ; print '(a,f13.5)','xlength=',xlength

      do n=1,nrsamp
!        Initial guess
         iconv(n)=0
         xsamp(n)=xsampini(n)
         qsamp(n)=qsampini(n)
         ysamp(n)=func(xsamp(n))+qsamp(n)
!        Prior in cost function 
         sf(1)=x0   
         sf(2)=0.0
      enddo
      call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
      Czz(1,1)=Cxx; Czz(2,2)=Cqq; Czz(1,2)=Cqx; Czz(2,1)=Cqx
      CIzz(1,1)=Cqq; CIzz(2,2)=Cxx; CIzz(1,2)=-Cqx; CIzz(2,1)=-Cqx
      CIzz=CIzz/(Cxx*Cqq-Cqx*Cqx)

      call getcaseid(caseid,'EnSTEIN',alphageo,nmda,esamp,sigw,0)
      sumconv=0
      do i=1,50 !maxiesit
         if (mod(i,10) == 0) then
            write(*,'(a,i4,a)')'i=',i,'...'
         endif
         call cov(Pxx,Pyy,Pqq,Pyx,Pqy,Pqx,xsamp,ysamp,qsamp,nrsamp)
         Pzz(1,1)=Pxx; Pzz(2,2)=Pqq; Pzz(1,2)=Pqx; Pzz(2,1)=Pqx
         PIzz(1,1)=Pqq; PIzz(2,2)=Pxx; PIzz(1,2)=-Pqx; PIzz(2,1)=-Pqx
         PIzz=PIzz/(Pxx*Pqq-Pqx*Pqx)
         Pzy(1,1)=Pyx; Pzy(2,1)=Pqy
         Pyz(1,1)=Pyx; Pyz(1,2)=Pqy

         do n=1,nrsamp
            si(1,n)=xsamp(n)
            si(2,n)=qsamp(n)
         enddo

         do n=1,nrsamp
! Prior term
            call dgemm('N','N',ndim,1,ndim,1.0,CIzz,ndim,si(:,n)-sf(:),ndim,0.0,grad(:,n),ndim)
! C_d^{-1} (y_j - d) with unperturbed data
            fac=(ysamp(n)-d)/sigo**2   
! Data term in EnStein  G=Pxx^{-1} Pxy
            call dgemm('N','N',ndim,1,ndim,fac,PIzz,ndim,Pzy,ndim,1.0,grad(:,n),ndim)
! Data term in Stein with analytic G
!            grad(1,n)=grad(1,n)+dfunc(si(1,n))*fac
!            grad(2,n)=grad(2,n)+fac
         enddo


         do n=1,nrsamp
            newgrad(:,n)=0.0
            do k=1,100
               call random_number(tmp)
               j=nint(tmp*real(nrsamp-1)+1.0)
               if (k==1) j=n
               newgrad(:,n)=newgrad(:,n)+&
                       steinkern(si(:,n),si(:,j),xlength,2)*  &
                       (grad(:,j) + (2.0*(si(:,j) - si(:,n))/xlength))
            enddo
            newgrad(:,n)=newgrad(:,n)/real(100.0)
         enddo

         do n=1,nrsamp
            si(:,n) = si(:,n) - 0.25*newgrad(:,n)
            xsamp(n)=si(1,n)
            qsamp(n)=si(2,n)
            ysamp(n)=func(xsamp(n)) + qsamp(n)
         enddo

         if (mod(i,i)==0) call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx,i)

      enddo
      samples(:,1,4)=xsamp(:)
      samples(:,2,4)=qsamp(:)

!    Recomputing ysamp with some noise for nicer plotting
      if (sigw < sigq) then
         do n=1,nrsamp
            ysamp(n)=ysamp(n)+sigq*normal()
         enddo
      endif
      call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
      call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
      call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
!      call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
      write(*,'(a)')'Stein analysis completed'
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! POST PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do n=1,min(ndim,2)
   do j=0,4
      call avevar(samples(1:nrsamp,n,j),nrsamp,ave,var)
      if (var > 0.0) write(*,'(4a,2g16.8,a,8g16.8)')method(j),' ',variable(n),' :',ave,sqrt(var),'  Samples:',samples(1:8,n,j)
   enddo
   enddo
end program 
