PROGRAM GGsimul 


implicit none


! Project: Gender Gaps 
! Author: Clement Joubert
! Date: December 14 2015

! This code simulates salary distributions at the world bank for the purpose of decomposing pay gaps as des-
! cribed in :
! 	"gender-pay-gaps.pdf"

! Tasks:
! 	(1) Set program specifications
!	(2) Import salary distributions
!	(3) Loop over years and types simulating salary distributions
!	(4) Output simulations

! Specification variables 
character*200,parameter:: path=""
character*200:: pathinput
character*200:: pathoutput
integer,parameter:: tenure      = 0   				! toggles the Career Gap decomposition
integer,parameter:: yp          = 0				! treats yp as a separate category
integer,parameter:: year1       = 1987
integer,parameter:: yearT       = 2027 
integer,parameter:: lastyeardata= 2015 
integer,parameter:: nyrs        = yearT-year1+1 
integer,parameter:: nyrsdata    = lastyeardata-year1+1 
integer,parameter:: ntypes      = 24				! number of employee groups
integer,parameter:: nbins       = 500				! number of salary bins
real(8),parameter:: binsize     = 1.0
integer,parameter:: nclones     = 10
real(8),parameter:: seed        = 145677.0d-0
integer,parameter:: regdum      = 1

! Import_salary variables
integer:: nend(nyrs,ntypes), nentry(nyrs,ntypes), nexit(nyrs,ntypes),ninter(nyrs,ntypes),tempint,gap,gaptype(ntypes)
real(8):: muend(nyrs,ntypes),     sigend(nyrs,ntypes),&
        & muentry(nyrs,ntypes), sigentry(nyrs,ntypes),&
        & muexit(nyrs,ntypes),   sigexit(nyrs,ntypes), &
        & muint(nyrs,ntypes), sigint(nyrs,ntypes),tempr
real(8):: alphagrowth(nyrs,ntypes), &
         & betagrowth(nyrs,ntypes), &
         &  siggrowth(nyrs,ntypes), &
         & noincrease(nyrs,ntypes), &
         & entryrate(nyrs,ntypes), &
         & exitrate(nyrs,ntypes), &
         & alphaexit(nyrs,ntypes), betaexit(nyrs,ntypes), beta2exit(nyrs,ntypes)

! Output statistics
integer:: nend_sim(nyrs,ntypes),nint_sim(nyrs,ntypes),nintdif(nyrs,ntypes)
real(8):: muend_sim(nyrs,ntypes),muenddif(nyrs,ntypes),&
        & sigend_sim(nyrs,ntypes), sigenddif(nyrs,ntypes), &
        & muint_sim(nyrs,ntypes),muintdif(nyrs,ntypes),&
        & sigint_sim(nyrs,ntypes), sigintdif(nyrs,ntypes),nenddif(nyrs,ntypes)

! internal variables
integer:: ss,typecor(ntypes),yy,tt,bin,idum,newbin,i,yr,tp,nolegacy,noentrygap,noexitgap,nogrowthgap, nocompgap 

integer::   sal_sim(nyrs,ntypes,nbins),&
          &               sal_sim_int(nyrs,ntypes,nbins), sal_data(nyrs,ntypes,nbins), &
          &                 entry_sim(nyrs,ntypes,nbins),  entry_data(nyrs,ntypes,nbins),&
          &                  exit_sim(nyrs,ntypes,nbins), exit_data(nyrs,ntypes,nbins)
real(8)::median_sal(nbins),sal,raise,newsal,sig,mu,draw(1), temp(nbins),bub,blb
real(8),allocatable:: raisedraw(:),saldraw(:),noincreasedraw(:)
real(8)::probexit, res,test, test2,test3,test4,test5
character*6::specification

! misc.

    median_sal=(/(i,i=1,nbins)/)*binsize+binsize/2
    if (tenure.eq.1) then
        pathinput=adjustl(trim(path))//"./Estimation/"//"CareerGap/"
        pathoutput=adjustl(trim(path))//"Simulations/"//"CareerGap/"
    elseif (tenure.eq.0) then
        pathinput=adjustl(trim(path))//"./Estimation/"//"AggregateGap/"
        pathoutput=adjustl(trim(path))//"Simulations/"//"AggregateGap/"
    endif
    

! TASK 1 - Set program specifications
    do gap=1,2

    ! type correspondence for counterfactual (i.e. simulate type x with parameters from type y)
    
    if (gap.eq.1) then ! gender gap
            typecor=(/(0,i=1,(ntypes/2)),(i,i=1,(ntypes/2))/)
    else if (gap.eq.2) then ! country part gap
            typecor=(/(0,i=1,(ntypes/4)),(i,i=1,(ntypes/4)),(0,i=1,(ntypes/4)),(i+ntypes/2,i=1,(ntypes/4))/)
    endif
    
    do noexitgap=0,1
    do nogrowthgap=0,1
    do noentrygap=0,1
    do nocompgap=0,1
    do nolegacy=0,1
    
    idum = dint(-1*abs(seed))
    
    if (noentrygap+noexitgap+nogrowthgap+nolegacy+nocompgap.gt.0) then 
       ! cycle !skip counterfactuals
    endif

    if (noexitgap.lt.nogrowthgap.or.nogrowthgap.lt.noentrygap.or.noentrygap.lt.nocompgap.or.nocompgap.lt.nolegacy) then 
        cycle !do only preferred order
    endif

    print*, noexitgap,nogrowthgap, noentrygap, nocompgap, nolegacy
    write(fmt='(6I1)',unit=specification) gap,nolegacy, nocompgap, noentrygap, nogrowthgap, noexitgap


! TASK 2 - Import data salary distributions
!   note: eventually/alternatively may import whole data set of salaries
!
!   nend, muend, sigend:                   distrib. of salaries at end of period
!   nentry, muentry, sigentry:             distrib. of salaries among entrants
!   nexit, muexit, sigexit:                distrib. of salaries among attriters
!   alphagrowth, betagrowth, siggrowth:    distrib. of logsalary growth
!   alphaexit, betaexit, beta2exit:        coefficients in probability of exit model
   
    call import_salary  (nbins,year1,nyrs,nyrsdata, ntypes, &
                    &    ninter, muint, sigint, &         !   end = end of period
                    &    nend, muend, sigend, &         !   end = end of period
                    &    nentry, muentry, sigentry, &
                    &    nexit, muexit, sigexit, &
                    &    alphagrowth, betagrowth, siggrowth, &
                    &    alphaexit, betaexit, beta2exit, &
                    &    noincrease,pathinput,entry_data,exit_data,sal_data)
    close(1)
    
  
    ! Simulate "nclones" copies of each employee
    sal_data=sal_data*nclones
    entry_data=entry_data*nclones
    exit_data=exit_data*nclones
    nend=nend*nclones
    ninter=ninter*nclones
    nentry=nentry*nclones
    nexit=nexit*nclones
    ! nentry=sum(entry_data,3)
    ! nexit=sum(exit_data,3)
    
    ! Change units to 000s
!    muentry=muentry/1000
!    sigentry=sigentry/1000
!    muexit=muexit/1000
!    sigexit=sigexit/1000
!    muend=muend/1000
!    sigend=sigend/1000
!    muint=muint/1000
!    sigint=sigint/1000
    
!    do yy=1,nyrs
!    do tt=1,ntypes
!        muentry(yy,tt)=dot_product(entry_data(yy,tt,:),median_sal(:))
!        muentry(yy,tt)=muentry(yy,tt)/nentry(yy,tt)
!        muexit(yy,tt)=dot_product(exit_data(yy,tt,:),median_sal(:))
!        muexit(yy,tt)=muexit(yy,tt)/nexit(yy,tt)
!    enddo
!    enddo

    ! Compute entry and exit rates
    entryrate=nentry*1.0d-0/ninter
    exitrate=nexit*1.0d-0/ninter


! TASK 3 - Simulate salary distributions one-year-ahead for each type of employee
!
! Simulations are stored in sal_sim(ntypes,nyrs,nbins). 
! Entries correspond to the number of individuals of a given type, 
! with salaries within the range of a given bin in a given year


    ! 0. Initialize year 1 and entries using data
    sal_sim=0
    sal_sim_int=0
    entry_sim(1,:,:)=entry_data(1,:,:)

    ! 0bis. Use counterfactual distributions if applicable
    	
    do tt=1,ntypes
    if (typecor(tt).gt.0) then
        if (nocompgap==1) then
            nentry       (:,tt)=nentry      (:,typecor(tt))            
        endif
        if (noentrygap==1) then 
            muentry      (:,tt)=muentry     (:,typecor(tt))
            sigentry     (:,tt)=sigentry    (:,typecor(tt))
        endif
        if (noexitgap==1) then 
            alphaexit    (:,tt)=alphaexit   (:,typecor(tt))
            betaexit     (:,tt)=betaexit    (:,typecor(tt))
            beta2exit    (:,tt)=beta2exit   (:,typecor(tt))
        endif
        if (nogrowthgap==1) then 
            alphagrowth (:,tt)=alphagrowth  (:,typecor(tt))
            betagrowth  (:,tt)=betagrowth   (:,typecor(tt))
            siggrowth   (:,tt)=siggrowth    (:,typecor(tt))
            noincrease  (:,tt)=noincrease   (:,typecor(tt))
        endif
        if (nolegacy==1) then 
            entry_sim (1,tt,:)=entry_sim  (1,typecor(tt),:)
        endif
    endif
    enddo

    ! 1. Simulate salaries one year ahead
    ! 
    ! Use data for number of entries and exits as well as distrib. of salaries for
    ! entrants, exiters and distrib. of salary growth

    do yy=1,nyrs       
    yr=yy+year1-1       !yr denotes calendar year, while yy is a counter
    
    ! If yy is larger than last year of data, use parameters from last year of data

    if (yr.gt.lastyeardata) then
	nentry(yy,:)=nentry(yy-1,:)
	muentry(yy,:)=muentry(yy-1,:)
	sigentry(yy,:)=sigentry(yy-1,:)
	alphaexit(yy,:)=alphaexit(yy-1,:)
	betaexit(yy,:)=betaexit(yy-1,:)
	beta2exit(yy,:)=beta2exit(yy-1,:)
	alphagrowth(yy,:)=alphagrowth(yy-1,:)
	betagrowth(yy,:)=betagrowth(yy-1,:)
	siggrowth(yy,:)=siggrowth(yy-1,:)
	noincrease(yy,:)=noincrease(yy-1,:)
    endif
    

    do tt=1,ntypes

    ! a. Salary growth
    
        if (yy>1)    call salgrowth(nyrs,ntypes,binsize,nbins,sal_sim,sal_sim_int,idum,noincrease,alphagrowth,betagrowth,siggrowth,tt,yy)
       
    ! b. Entries and exits
 
        res=0.0d-0            
        do bin=1,nbins

            ! Add new entrants
            if (yy>1.or.tenure==1) then   !For 1987, use exact salary bins, for later years, simulate based on estimated average and sd 
                    blb= 0.0d-0 + (bin-1)*binsize
                    if (bin.eq.1) blb=-999999999999.0d-0
                    bub= bin*binsize
                if (sigentry(yy,tt).gt.0.0d-0)  then 
                    if (bin.eq.nbins) bub=999999999999.0d-0
                    call hist(entry_sim(yy,tt,bin),nentry(yy,tt),muentry(yy,tt),sigentry(yy,tt),blb,bub,res)
                elseif (sigentry(yy,tt).eq.0.0d-0.and.muentry(yy,tt).gt.blb.and.muentry(yy,tt).le.bub) then
                    entry_sim(yy,tt,bin)=nentry(yy,tt) !in case the distribution of salary entrants is degenerate
                else
                    entry_sim(yy,tt,bin)=0
                endif
            endif
            
            sal_sim(yy,tt,bin)=sal_sim_int(yy,tt,bin)+entry_sim(yy,tt,bin)
            
            ! For each salary bin, simulate how many employees leave
            tempr=probexit(median_sal(bin),alphaexit(yy,tt),betaexit(yy,tt),beta2exit(yy,tt),regdum)
            tempr=tempr*sal_sim(yy,tt,bin)
            tempint=idnint(tempr)
            exit_sim(yy,tt,bin)=tempint
            
            
            ! Combine salaries of entrants exiters and employee stock 
            sal_sim(yy,tt,bin)=sal_sim(yy,tt,bin)-exit_sim(yy,tt,bin)
        
        enddo
    enddo
    enddo
    
! TASK 4 - Compare actual and simulated end-of-year distributions
    
    call fit(pathoutput,sal_sim, sal_data, sal_sim_int,entry_sim,entry_data,exit_sim,exit_data,nend,ninter,nexit,nentry,binsize,nbins,nyrs, ntypes, muend,sigend,muint,sigint,muentry,muexit,specification,tenure)
    
enddo
enddo
enddo
enddo
enddo
enddo


print*,"the end"

END PROGRAM GGsimul
	
subroutine hist(freq,n,mu,sigma,blb,bub,res)

external alnorm

!*****************************************************************************80
!  Modified:
!
!    17 December 2015
!
!  Author:
!
!    Clement Joubert
!
!  Parameters:
!   Inputs
    integer:: n !is the size of the population
    real(8):: mu, sigma !are the parameters of the normal distribution
    real(8):: res !residual frequency from rounding the frequency of the previous bin
!
!    Output
    integer:: freq !frequency of the bin
!****************************************************************************80
    real(8):: binproba,blb,bub,sblb, sbub, alnorm, lcdf, ucdf

    
    sblb=(blb-mu)/sigma   !standardized bin lower bound
    sbub=(bub-mu)/sigma   !standardized bin upper bound
    
    lcdf=alnorm(sblb,.false.)   !cdf at lower bound
    ucdf=alnorm(sbub,.false.)   !cdf at upper bound
    
    binproba= ucdf-lcdf 
    
    if (blb.gt.18) then
        freq = nint(n*binproba + res)  
    else
        freq = 0
    endif
    res = n*binproba + res - freq
    
endsubroutine


function probexit(sal,alpha,beta,beta2,regdum)

!*****************************************************************************80
!  Modified:
!
!
!  Author:
!
!    Clement Joubert
!
!  Parameters:
!   Inputs
    real(8):: temp, alpha, beta, beta2, sal
    integer:: regdum !if using a linear probability model
!    Output

    real(8):: probexit 
!****************************************************************************80   

if (regdum.eq.0) then
    temp=exp(-alpha- beta*sal - beta2*sal*sal)  
    probexit = 1/(1+temp)  
elseif (regdum.eq.1) then
    temp=beta*(sal-alpha)+beta2
    probexit= dmin1(1.0d-0,dmax1(0.0d-0,temp))
endif

    
endfunction



subroutine Import_salary (nbins, year1, nyrs, nyrsdata, ntypes, &
                    &    ninter,muinter,siginter, &
                    &    nend, muend, sigend, &         !   end = end of period
                    &    nentry, muentry, sigentry, &
                    &    nexit, muexit, sigexit, &
                    &    alphagrowth, betagrowth, siggrowth, &
                    &    alphaexit, betaexit, beta2exit, &
                    &    noincrease,pathinput,entry_data,exit_data,sal_data)
implicit none


! arguments
integer:: nyrs,nyrsdata, ntypes, nstats,nbins,bin, year1
integer:: nend(nyrs,ntypes), nentry(nyrs,ntypes), nexit(nyrs,ntypes),ninter(nyrs,ntypes)
real(8):: muend(nyrs,ntypes), sigend(nyrs,ntypes),&
        & muentry(nyrs,ntypes), sigentry(nyrs,ntypes),&
        & muexit(nyrs,ntypes), sigexit(nyrs,ntypes), &
        & muinter(nyrs,ntypes), siginter(nyrs,ntypes), &
        & alphagrowth(nyrs,ntypes), betagrowth(nyrs,ntypes), siggrowth(nyrs,ntypes),sebetagrowth(nyrs,ntypes), &
        & noincrease(nyrs,ntypes), &
        & alphaexit(nyrs,ntypes), betaexit(nyrs,ntypes), beta2exit(nyrs,ntypes)

integer::   sal_data(nyrs,ntypes,nbins),&
          &               sal_data_int(nyrs,ntypes,nbins),&
          &                 entry_data(nyrs,ntypes,nbins),&
          &                  exit_data(nyrs,ntypes,nbins)
character*200:: pathinput,filename
character*2::ttchar


! internal
integer:: yy,tt,stat
character:: blah
    
    entry_data=0
    exit_data=0
    sal_data=0
    
    
    open(file=adjustl(trim(pathinput))//"simstats.txt",unit=1,position="rewind")
    read(1,*)
    do tt=1,ntypes
    !do tt=1,24
        do yy=1,nyrsdata
            write(ttchar,'(I2)') tt
            read(1,*) blah, nentry(yy,tt), muentry(yy,tt), sigentry(yy,tt), noincrease(yy,tt),  betagrowth(yy,tt),  sebetagrowth(yy,tt), alphagrowth(yy,tt), siggrowth(yy,tt),  betaexit(yy,tt), beta2exit(yy,tt), alphaexit(yy,tt), nend(yy,tt),muend(yy,tt),sigend(yy,tt), nexit(yy,tt), muexit(yy,tt), sigexit(yy,tt),ninter(yy,tt),muinter(yy,tt),siginter(yy,tt)
    enddo
    enddo
       
    open(file=adjustl(trim(pathinput))//"salhist.txt",unit=2)
    bin=0
    stat=0      
    do while(bin.le.nbins.and.stat.ge.0)
    read(2,iostat=stat,fmt=*) bin, tt, yy, entry_data(yy-year1+1,tt,bin), exit_data(yy-year1+1,tt,bin), sal_data(yy-year1+1,tt,bin)
    !write(*,*) bin
    enddo
    close(2)
    
endsubroutine Import_salary

subroutine salgrowth(nyrs,ntypes,binsize,nbins,sal_sim,sal_sim_int,idum,noincrease,alphagrowth,betagrowth,siggrowth,tt,yy)
implicit none

integer:: nyrs 
integer:: ntypes 
integer:: nbins
real(8):: binsize

! Import_salary variables
integer:: nend(nyrs,ntypes), nentry(nyrs,ntypes), nexit(nyrs,ntypes),ninter(nyrs,ntypes)
real(8):: muend(nyrs,ntypes),     sigend(nyrs,ntypes),&
        & muentry(nyrs,ntypes), sigentry(nyrs,ntypes),&
        & muexit(nyrs,ntypes),   sigexit(nyrs,ntypes), &
        & muint(nyrs,ntypes), sigint(nyrs,ntypes)
real(8):: alphagrowth(nyrs,ntypes), &
         & betagrowth(nyrs,ntypes), &
         &  siggrowth(nyrs,ntypes), &
         & noincrease(nyrs,ntypes)

! Output statistics
integer:: nend_sim(nyrs,ntypes),nenddif(nyrs,ntypes),nint_sim(nyrs,ntypes),nintdif(nyrs,ntypes)
real(8):: muend_sim(nyrs,ntypes),muenddif(nyrs,ntypes),&
        & sigend_sim(nyrs,ntypes), sigenddif(nyrs,ntypes), &
        & muint_sim(nyrs,ntypes),muintdif(nyrs,ntypes),&
        & sigint_sim(nyrs,ntypes), sigintdif(nyrs,ntypes)

! internal variables
integer:: ss,yy,tt,bin,hist,idum,newbin,i,yr,tp,binfreq
integer::   sal_sim(nyrs,ntypes,nbins),&
          &               sal_sim_int(nyrs,ntypes,nbins),&
          &                 entry_sim(nyrs,ntypes,nbins),&
          &                  exit_sim(nyrs,ntypes,nbins)
real(8)::median_sal(nbins),sal,raise,newsal,sig,mu,draw(1), temp(nbins)
real(8),allocatable:: raisedraw(:),saldraw(:),noincreasedraw(:)

 
        GROWTHLOOP: do bin=1,nbins
	    binfreq=sal_sim(yy-1,tt,bin)
            if (binfreq.gt.0) then
            
                allocate(raisedraw(binfreq))
                allocate(saldraw(binfreq))
                allocate(noincreasedraw(binfreq))
            
                call gasdev(idum,raisedraw,binfreq)
                call ran1(idum,saldraw,binfreq)
                call ran1(idum,noincreasedraw,binfreq)
            
                ! Consider each observation in each bin of last year's distribution
                BINSHIFTS: do ss=1,binfreq
                    ! Draw a salary from the bin
                    sal=(bin-1+saldraw(ss))*binsize
                    ! Draw salary growth and compute new salary
                    mu=alphagrowth(yy,tt)+sal*betagrowth(yy,tt)
                    sig=siggrowth(yy,tt)
                    raise=1+exp(raisedraw(ss)*sig+mu)
                    if (noincreasedraw(ss).lt.noincrease(yy,tt)) raise= 1       
                    newsal=sal*raise
                    ! Identify new salary bin for that simulated employee
                    newbin=min(nbins,1+int(newsal/binsize))
                    sal_sim_int(yy,tt,newbin)=sal_sim_int(yy,tt,newbin)+1
                    
                enddo BINSHIFTS
                
                deallocate(raisedraw)
                deallocate(saldraw)  
                deallocate(noincreasedraw) 
                 
            endif
           
         enddo GROWTHLOOP
         

endsubroutine salgrowth

subroutine fit(pathoutput,sal_sim,sal_data, sal_sim_int,entry_sim,entry_data,exit_sim,exit_data,nend,ninter,nexit,nentry,binsize,nbins,nyrs, ntypes, muend,sigend,muint,sigint,muentry,muexit,specification,tenure)
implicit none

character*200:: pathoutput 
integer:: nyrs 
integer:: ntypes 
integer:: nbins
real(8):: binsize
character*6:: specification
  
   ! Import_salary variables
integer:: nend(nyrs,ntypes), nentry(nyrs,ntypes), nexit(nyrs,ntypes),ninter(nyrs,ntypes), tenure
real(8):: muend(nyrs,ntypes),     sigend(nyrs,ntypes),&
        & muentry(nyrs,ntypes), sigentry(nyrs,ntypes),&
        & muexit(nyrs,ntypes),   sigexit(nyrs,ntypes), &
        & muint(nyrs,ntypes), sigint(nyrs,ntypes)
real(8):: alphagrowth(nyrs,ntypes), &
         & betagrowth(nyrs,ntypes), &
         &  siggrowth(nyrs,ntypes), &
         & noincrease(nyrs,ntypes)

! Output statistics
integer:: nend_sim(nyrs,ntypes),nenddif(nyrs,ntypes),nint_sim(nyrs,ntypes),nintdif(nyrs,ntypes),nentry_sim(nyrs,ntypes),nexit_sim(nyrs,ntypes),nentrydif(nyrs,ntypes),nexitdif(nyrs,ntypes)
real(8):: muend_sim(nyrs,ntypes),muenddif(nyrs,ntypes),&
        & sigend_sim(nyrs,ntypes), sigenddif(nyrs,ntypes), &
        & muint_sim(nyrs,ntypes),muintdif(nyrs,ntypes), &
        & sigint_sim(nyrs,ntypes), sigintdif(nyrs,ntypes), &
        & muexit_sim(nyrs,ntypes),muentry_sim(nyrs,ntypes), muentrydif(nyrs,ntypes), muexitdif(nyrs,ntypes)

! internal variables
integer:: ss,yy,tt,bin,hist,idum,newbin,i,yr,tp
integer::   sal_sim(nyrs,ntypes,nbins), sal_data(nyrs,ntypes,nbins),&
          &               sal_sim_int(nyrs,ntypes,nbins),&
          &                 entry_sim(nyrs,ntypes,nbins),entry_data(nyrs,ntypes,nbins),&
          &                  exit_sim(nyrs,ntypes,nbins),exit_data(nyrs,ntypes,nbins)
real(8)::median_sal(nbins),seed,sal,raise,newsal,sig,mu,draw(1), temp(nbins)
real(8),allocatable:: raisedraw(:),saldraw(:),noincreasedraw(:)

   
    ! Compute number of simulated individuals per type for each year
    nend_sim(:,:)=sum(sal_sim,3)
    nenddif=nend_sim-nend
    
    nint_sim(:,:)=sum(sal_sim_int,3)
    nintdif=nint_sim-ninter
    
    nentry_sim(:,:)=sum(entry_sim,3)
    nentrydif=nentry-nentry_sim
    
    nexit_sim(:,:) =sum(exit_sim,3)
    nexitdif=nexit-nexit_sim

    median_sal=(/(i,i=1,nbins)/)*binsize-binsize/2

    ! Compute means of the distribution of salaries before and after entries and exits   
    do yy=1,nyrs
    do tt=1,ntypes
        muend_sim(yy,tt)=dot_product(sal_sim(yy,tt,:),median_sal(:))
        muend_sim(yy,tt)=muend_sim(yy,tt)/nend_sim(yy,tt)
        muint_sim(yy,tt)=dot_product(sal_sim_int(yy,tt,:),median_sal(:))
        muint_sim(yy,tt)=muint_sim(yy,tt)/nint_sim(yy,tt)
        muexit_sim(yy,tt)=dot_product(exit_sim(yy,tt,:),median_sal(:))
        muexit_sim(yy,tt)=muexit_sim(yy,tt)/nexit_sim(yy,tt)
        muentry_sim(yy,tt)=dot_product(entry_sim(yy,tt,:),median_sal(:))
        muentry_sim(yy,tt)=muentry_sim(yy,tt)/nentry_sim(yy,tt)

    enddo
    enddo
    muenddif=muend_sim-muend
    muintdif=muint_sim-muint
    muentrydif=muentry_sim-muentry
    muexitdif=muexit_sim-muexit
    
    ! Compute standard deviations of the distribution of salaries before and after entries and exits
    do yy=1,nyrs
    do tt=1,ntypes
        temp(:)=(median_sal(:)*1.0d-0-muend_sim(yy,tt))
        temp(:)=temp(:)**2
        temp(:)=sal_sim(yy,tt,:)*temp(:)
        sigend_sim(yy,tt)=sum(temp(:))/nend_sim(yy,tt)
        sigend_sim(yy,tt)=sqrt(sigend_sim(yy,tt))
        
        temp(:)=(median_sal(:)*1.0d-0-muint_sim(yy,tt))
        temp(:)=temp(:)**2
        temp(:)=sal_sim_int(yy,tt,:)*temp(:)
        sigint_sim(yy,tt)=sum(temp(:))/nint_sim(yy,tt)
        sigint_sim(yy,tt)=sqrt(sigint_sim(yy,tt))
    enddo
    enddo
    sigenddif=sigend_sim-sigend
    sigintdif=sigint_sim-sigint
    
    
    open(unit=1,file=trim(adjustl(pathoutput))//"simulationfit"//specification//".txt")
    write(1,'(A9,20A7)')      "year ","type ",&
                                  & "mend_d ","mend_s ","mend_f ",  "send_f ", &
                                  & "nend_d ","nend_s ",&
                                  & "mout_d ", "mout_s ","mout_f ",&
                                  & "nout_d ", "nout_s ",& 
                                  & "min_d ", "min_s ", "min_f ", &
                                  & "nin_d ","nin_s ",&
                                  & "mint_d ","mint_s ", "mint_f "
    do tt=1,ntypes
    do yy=1,nyrs
        write(1,fmt='(2I7,4f7.1,2I7,3f7.1,2I7,3f7.1,2I7,3f7.1)')1986*(1-tenure)+yy,tt,&
                    & muend(yy,tt),muend_sim(yy,tt),muenddif(yy,tt), sigenddif(yy,tt)    ,&  
                    & nend(yy,tt),nend_sim(yy,tt),&
                    & muexit(yy,tt),muexit_sim(yy,tt),muexitdif(yy,tt),&
                    & nexit(yy,tt),nexit_sim(yy,tt),&
                    & muentry(yy,tt),muentry_sim(yy,tt),muentrydif(yy,tt),&
                    & nentry(yy,tt),nentry_sim(yy,tt), &
                    & muint(yy,tt),muint_sim(yy,tt),muintdif(yy,tt)
                  
    enddo
    enddo
    close(1)
    
end subroutine


function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
  real ( kind = 8 ) alnorm
  real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
  real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
  real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
  real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
  real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
  real ( kind = 8 ), parameter :: con = 1.28D+00
  real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
  real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
  real ( kind = 8 ), parameter :: ltone = 7.0D+00
  real ( kind = 8 ), parameter :: p = 0.398942280444D+00
  real ( kind = 8 ), parameter :: q = 0.39990348504D+00
  real ( kind = 8 ), parameter :: r = 0.398942280385D+00
  logical up
  logical upper
  real ( kind = 8 ), parameter :: utzero = 18.66D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 &
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
  end if

  return
end



SUBROUTINE ran1sub(idum,x) 
IMPLICIT NONE
INTEGER, PARAMETER :: K4B=selected_int_kind(9)
INTEGER(K4B), INTENT(INOUT) :: idum
REAL(8) :: x 
INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647,IQ=127773,IR=2836
REAL(8), SAVE :: am
INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
if (idum <= 0 .or. iy < 0) then
     am=nearest(1.0d-0,-1.0d-0)/IM
     iy=ior(ieor(888889999,abs(idum)),1)
     ix=ieor(777755555,abs(idum))
     idum=abs(idum)+1
end if
ix=ieor(ix,ishft(ix,13))
ix=ieor(ix,ishft(ix,-17))
ix=ieor(ix,ishft(ix,5))
k=iy/IQ
iy=IA*(iy-k*IQ)-IR*k
if (iy<0) iy=iy+IM
x=am*ior(iand(IM,ieor(ix,iy)),1)
END SUBROUTINE ran1sub


SUBROUTINE ran1(idum,vv,n)
 
implicit none

external ran1sub

REAL(8) vv(n),x
INTEGER kk,n,idum

vvloop: do kk = 1,n
     call ran1sub(idum,x)
     vv(kk) = x
end do vvloop

idum = idum+2

END SUBROUTINE ran1
 


!------------------------


       SUBROUTINE GASDEV(idum,vv,n)

! This function returns a normally distributed deviate with
! zero mean and unit variance, using ran1(IDUM) as the
! source of uniform deviates

       implicit none

       integer iset,n,k,gotoc,idum
       real(8) v1,v2,r,fac,gset,gasdev1
       real(8) temp1(2),vv(n),ix1,ix2

       external ran1

!       write(*,*) 'want',n,'random numbers'

       vvloop: do k = 1,n

       iset = 0
       v1 = 0.0d-0
       v2 = 0.0d-0
       r = 0.0d-0
       fac = 0.0d-0
       gset = 0.0d-0
       gasdev1 = 0.0d-0
       ix1=0.0d-0
       ix2=0.0d-0
       gotoc = 0

       if (iset.eq.0) then


1         call ran1(idum,temp1,2)

          v1 = (2.0d-0)*temp1(1)-1.0d-0 
          v2 = (2.0d-0)*temp1(2)-1.0d-0

!          write(*,*) 'v1',v1
!          write(*,*) 'v2',v2
!          write(*,*) 'idum',idum

          r = v1**2.0d-0+v2**2.0d-0
          if (r.ge.1.0d-0.or.r.eq.0.0d-0) then 
             gotoc = gotoc + 1
             if (gotoc.gt.n) then
                !write(*,*) 'error in gasdev'
             end if
             go to 1
          end if

          fac = dsqrt(-2.0d-0*dlog(r)/r)
          gset = v1*fac
          gasdev1 = v2*fac
          iset = 1
       else
          gasdev1 = gset
          iset = 0
       end if
       
       vv(k)=gasdev1
!       write(*,*) 'random normal number',k,vv(k)

       end do vvloop

       return
       end SUBROUTINE GASDEV





         








 
         








 

