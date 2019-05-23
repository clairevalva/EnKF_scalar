program iterative_smoothers
   use mod_inistat
   use mod_xyqgrid
   use m_marginalpdf
   use m_set_random_seed2
   use m_moments
   use m_avevar
   use m_func
   use m_cyyreg
   use m_aaprojection
   use m_iniens
   use m_es
   use m_ies
   use m_sies
   use m_esmda
   use m_enstein
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

   character(len=7) :: method(0:5)=['INI    ','ES     ','IES    ','SIES   ', 'ESMDA  ', 'EnSTEIN']
   character(len=1) :: variable(1:2)=['x','q']

   real, allocatable :: xsampini(:) 
   real, allocatable :: ysampini(:)
   real, allocatable :: qsampini(:)
   real, allocatable :: dpert(:)

   real, allocatable :: samples(:,:,:) ! for comparing ES, IES and SIES

   character(len=40) caseid

   integer i,j,k,m,n,loc(1),ival,jval
   real jx,djx,djq,ddjx
   real dgx,dgq

   real ave,var
   real sump,fac,tmp
   real G
   real hessian, cvec(2) , psca


   character(len=2) ca
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call set_random_seed2()
   print *,'Normal number', normal()

   open(10,file='infile.in')
      read(10,*)esamp       ; print '(a,i4,i12)',   'number  of samples 10^x    :',esamp,10**esamp
      read(10,'(a)')ca      
      if (ca /= '#1') then
         print *,'#1: error in infile.in'
         stop
      endif
      read(10,*)xa            ; print '(a,f10.3)',    'xa                         :',xa
      read(10,*)xb            ; print '(a,f10.3)',    'xb                         :',xb
      read(10,*)qa            ; print '(a,f10.3)',    'qa                         :',qa
      read(10,*)qb            ; print '(a,f10.3)',    'qb                         :',qb
      read(10,*)ya            ; print '(a,f10.3)',    'ya                         :',ya
      read(10,*)yb            ; print '(a,f10.3)',    'yb                         :',yb
      read(10,'(a)')ca
      if (ca /= '#2') then
         print *,'#2: error in infile.in'
         stop
      endif
      read(10,*)x0            ; print '(a,f10.3)',    'prior estimate of x        :',x0
      read(10,*)siga          ; print '(a,f10.3)',    'err std dev of prior       :',siga
      read(10,*)d             ; print '(a,f10.3)',    'measurement of y           :',d
      read(10,*)sigo          ; print '(a,f10.3)',    'measurement std dev        :',sigo ; cdd=sigo**2
      read(10,*)sigw          ; print '(a,f10.3)',    'model std dev              :',sigw
      read(10,'(a)')ca
      if (ca /= '#3') then
         print *,'#3: error in infile.in'
         stop
      endif
      read(10,*)beta          ; print '(a,f10.3)',    'model parameter beta       :',beta
      read(10,*)funcmode      ; print '(a,tr7,i3)',   'function to use            :',funcmode
      read(10,'(a)')ca
      if (ca /= '#4') then
         print *,'#4: error in infile.in'
         stop
      endif
      read(10,*)lcyyreg       ; print '(a,tr10,l1)',  'Regression for Cyy         :',lcyyreg
      read(10,*)lactivepdf    ; print '(a,tr10,l1)',  'Activate pdf printing      :',lactivepdf
      read(10,'(a)')ca
      if (ca /= '#5') then
         print *,'#5: error in infile.in'
         stop
      endif
      read(10,*)lmda          ; print '(a,tr10,l1)',  'Run MDA                    :',lmda
      read(10,*)nmda          ; print '(a,tr7,i3)',   'number of mda iterations   :',nmda
      read(10,*)alphageo      ; print '(a,f10.3)',    'geometrical alpha value    :',alphageo
      read(10,'(a)')ca
      if (ca /= '#6') then
         print *,'#6: error in infile.in'
         stop
      endif
      read(10,*)lies          ; print '(a,tr10,l1)',  'Run IES                    :',lies
      read(10,*)maxiesit      ; print '(a,i5)',       'maxiesit                   :',maxiesit
      read(10,*)gamma_ies     ; print '(a,f10.3)',    'gamma_ies                  :',gamma_ies
      read(10,*)IESv          ; print '(a,i1)',       'IESv                       :',IESv
      read(10,'(a)')ca
      if (ca /= '#7') then
         print *,'#7: error in infile.in'
         stop
      endif
      read(10,*)lsies         ; print '(a,tr10,l1)',  'Run SIES                   :',lsies
      read(10,*)maxsiesit     ; print '(a,i5)',       'maxiesit                   :',maxsiesit
      read(10,*)gamma_sies    ; print '(a,f10.3)',    'gamma_sies                 :',gamma_sies
      read(10,'(a)')ca
      if (ca /= '#8') then
         print *,'#8: error in infile.in'
         stop
      endif
      read(10,*)lenstein      ; print '(a,tr10,l1)',  'Run EnSTEIN                :',lenstein
      read(10,*)maxensteinit  ; print '(a,i5)',       'maxiesit                   :',maxensteinit
      read(10,*)gamma_enstein ; print '(a,f10.3)',    'gamma_sies                 :',gamma_enstein
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
         

   nrsamp=10**esamp
   nrsamp=1*nrsamp

   allocate(qsampini(nrsamp), xsampini(nrsamp) , ysampini(nrsamp), dpert(nrsamp))
   allocate(samples(nrsamp,2,0:5)); samples=0.0


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
      cost(i)=(x(i)-x0)**2/siga**2 + (func(x(i))-d)**2/sigo**2    
   enddo
   sump=sum(prior(:))*dx
   prior=prior/sump
   call tecfunc('costf',cost,x,nx,'x','Cost')
   call tecfunc('prior',prior,x,nx,'x','Prior')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Datum pdf
   do j=1,ny
      datum(j)=exp(-0.5*(y(j)-d)**2/sigo**2)
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
                   -0.5*(y(j)-d)**2/sigo**2)    
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
                               -0.5*(d-func(x(i)))**2/sigo**2)
      else 
         margx(i)=margx(i)+exp(-0.5*(x(i)-x0 )**2/siga**2                   &
                               -0.5*(q(j)-0.0)**2/max(sigw,0.0001)**2  & 
                               -0.5*(d-func(x(i))-q(j))**2/sigo**2)
      endif
   enddo
   enddo
   margx=margx/(sum(margx(:))*dx)
   call tecfunc('margx',margx,x,nx,'x',caseid)
   call tecfunc('margy',margy,y,ny,'y',caseid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ensemble initialization
   call iniens(samples(1:nrsamp,1:2,0),xsampini,qsampini,ysampini,dpert,nrsamp,esamp)

! ES 
   call es(samples(1:nrsamp,1:2,1),xsampini,qsampini,dpert,nrsamp,esamp)

! IES
   if (lies) then
      call  ies(samples(1:nrsamp,1:2,2),xsampini,qsampini,nrsamp,esamp,dpert)
   endif

! Subspace IES
   if (lsies) then
      call sies(samples(1:nrsamp,1:2,3),xsampini,qsampini,nrsamp,esamp,dpert)
   endif

! ESMDA
   if (lmda) then
      call esmda(samples(1:nrsamp,1:2,4),xsampini,qsampini,nrsamp,esamp)
   endif

! EnSTEIN
   if (lenstein) then
      call enstein(samples(1:nrsamp,1:2,5),xsampini,qsampini,nrsamp,esamp)
   endif

! POST PROCESSING
!   call avevar(dpert,nrsamp,ave,var)
!   if (var > 0.0) write(*,'(4a,2g16.8,a,8g16.8)')'MEAS   ',' ',variable(1),' :',ave,sqrt(var),'  Samples:',dpert(1:8)
   do n=1,2
   do j=1,5
      call avevar(samples(1:nrsamp,n,j),nrsamp,ave,var)
      if (var > 0.0) write(*,'(4a,2g16.8,a,8g16.8)')method(j),' ',variable(n),' :',ave,sqrt(var),'  Samples:',samples(1:8,n,j)
   enddo
   enddo

! Printing some cost functions and the ES solution before and after (see also dead code in ies
!   call printcostf(samples(1:nrsamp,1,0),samples(1:nrsamp,1,1),dpert,nrsamp)

end program 
