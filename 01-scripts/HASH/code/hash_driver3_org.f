c Sample main driver program for using the focal mechanism subroutines.  
c Uses P polarities and S/P amplitude ratios.

      include 'param.inc'
      include 'rot.inc'

c variables for storing earthquake input information  
      integer icusp,icusp2          ! event ID
      real qlat,qlon,qdep           ! location
      real qmag                     ! magnitude
      integer iyr,imon,idy,ihr,imn  ! origin time, year, month, day, hour, minute
      real qsec                     ! origin time, seconds
      real seh, sez                 ! location error, horizontal & vertical 
      real rms                      ! residual travel time error 
      real terr                     ! origin time error 
      character*1 evtype            ! event type
      character*1 magtype           ! magnitude type
      character*1 locqual           ! location quality
      character*1 cns,cew           ! north/south and east/west codes
c
c variables for polarity information, input to HASH subroutines
      character*4 sname(npick0)                        ! station name
      character*3 scomp(npick0)                        ! station component
      character*2 snet(npick0)                         ! network code
      character*1 pickpol,pickonset                    ! polarity pick : U, u, +, D, d, or - ; onset, I or E
      integer p_pol(npick0),q_pol(npick0),spol                       ! polarity pick (+/-1), and reversal (+/-1)
      integer polpred(npick0,5)                        ! WSDW - predicted polarities for preferred solution
      real sppred(npick0,5)                         ! WSDW - predicted log(S/P) for preferred solution
      real sp_ratio(npick0),spin                       ! S/P ratio (log10)
      real p_azi_mc(npick0,nmc0),p_the_mc(npick0,nmc0) ! azimuth and takeoff angle for each trial
      integer index(nmc0)                              ! index into velocity models, for each trial
      real qlat2(nmc0), qlon2(nmc0)                    ! WSDW - add epicenter errors
      real qdep2(nmc0)                                 ! new depth, for each trail


      integer nmc                                      ! number of trails with different azimuth and take-off angles
      integer npol,nppl,nspr                           ! number of observations, P polarities, and S/P ratios
c
c variables for set of acceptable mechanisms, output of HASH subroutines
      integer nout2                                    ! number of acceptable mechanisms returned
      integer nmult                                    ! number of solutions (1 if no mulitples)
      real str_avg(5),dip_avg(5),rak_avg(5)            ! solution(s)
      real f1norm(3,nmax0),f2norm(3,nmax0)             ! normal vectors to the fault planes
      real strike2(nmax0),dip2(nmax0),rake2(nmax0)     ! strike, dip and rake
      real var_est(2,5),var_avg(5)                     ! variance of each plane, average
      real mfrac(5),stdr(5),mavg(5)                    ! fraction misfit polarities, station distribution, WSDW log10(S/P) misfit       
      real prob(5)                                     ! probability true mechanism is "close" to preferred solution(s)
*      character*1 qual(5),mflag                        ! solution quality rating, multiple flag
      character*1 qual(5)         ! solution quality rating
c
c control parameters
*      integer npolmin                                  ! minimum number of observations WSDW not used so eliminate
      integer max_agap,max_pgap                        ! maximum azimuthal & "plungal" gap
      integer nextra0,nextra1                                 ! WSDW minimum value of nextra
      real delmax                                      ! maximum station distance
      real dang,dang2                                  ! grid angle to sample focal sphere
      integer maxout                                   ! max number of acceptable mechanisms output
      real badfrac                                     ! assumed rate of polarity error (fraction)
      real cangle                                      ! definition of "close" == 45 degrees
      real ratmin                                      ! minimum allowed signal to noise ratio
      real qbadfac                                     ! assumed noise in amplitude ratios, log10 (0.3 for factor of 2)
c
c file names
      character*100 outfile1,corfile,fpfile
      character*100 stfile,plfile,ampfile
      character*100 outfile2,outfile3,outfile4
      
      degrad=180./3.1415927
      rad=1./degrad
      
      print *,'Enter station list file'       
      read (*,'(a)') stfile

      print *,'Enter station polarity reversal file'       
      read (*,'(a)') plfile

      print *,'Enter station correction file'       
      read (*,'(a)') corfile

      print *,'Enter name of amplitude input file'       
      read (*,'(a)') ampfile

      print *,'Enter name of P-polarity input file'       
      read (*,'(a)') fpfile

      print *,'Enter output file name for focal mechanisms'
      read (*,'(a)') outfile1
      open (13,file=outfile1)

      print *,'Enter output file name for acceptable planes'
      read (*,'(a)') outfile2
      open (14,file=outfile2)
      
      print *,'Enter output file for obs and preds'
      read (*,'(a)') outfile3
      open (15,file=outfile3)
      
      print *,'Enter output file for s/p ratio'
      read (*,'(a)') outfile4
      open (16,file=outfile4)

      print *,'Enter maximum azimuthal gap (e.g., 90)'
      read *,max_agap

      print *,'Enter maximum takeoff angle gap (e.g., 60)'
      read *,max_pgap

      print *,'Enter grid angle for focal mech search, in degrees 
     &  (max ',dang0,')'
      read *,dang
      dang2=max(dang0,dang) ! don't do finer than dang0

      print *,'Enter number of trials (e.g., 30)'
      read *,nmc

      print *,'Enter maxout for focal mech. output (e.g., 500)'
      read *,maxout

      print *,'Enter minimum allowed signal to noise ratio'
      read *,ratmin
      
      print *,'Enter fraction polarities assumed bad'
      read *,badfrac
      
      print *,'Number of additional impulsive polarity misfits'
      read *,nextra0

      print *,'Number of additional emergent polarity misfits'
      read *,nextra1

      print *,'Enter the assumed noise in amplitude ratios, log10  
     &  (e.g. 0.3 for a factor of 2)'
      read *,qbadfac

      print *,'Enter maximum allowed source-station distance, 
     &         in km (e.g., 120)'
      read *,delmax

      print *,'Enter angle for computing mechanisms probability, 
     &         in degrees (e.g., 45)'
      read *,cangle

      print *,'Enter probability threshold for multiples (e.g., 0.1)'
      read *,prob_max

c make tables of takeoff angles for various velocity models
      ntab=10
      call MK_TABLE(ntab)

c read in earthquake location, etc      ! SCEDC format - ** YOU MAY NEED TO CHANGE THE INPUT FORMAT **
      open (11,file=fpfile,status='old')
120   continue
      read (11,125,end=505) iyr,imon,idy,ihr,imn,qsec,ilatd,cns,qlatm,
     &                ilond,cew,qlonm,qdep,seh,sez,qmag,icusp
* Modify polarity file event info to a more compact form
125   format (i4,4i2,f5.2,i2,a1,f5.2,i3,a1,f5.2,f5.2,1x,
     &                f5.2,1x,f5.2,1x,f4.2,i16)     
*125   format (i4,4i2,f5.2,i2,a1,f5.2,i3,a1,f5.2,f5.2,3x,46x,
*     &                f5.2,1x,f5.2,40x,f4.2,6x,i16)     
      qlat=real(ilatd)+(qlatm/60.0)
      if (cns.eq.'S') then
        qlat=-qlat
      end if
* Fix sign so East is positive 
*      qlon=-(real(ilond)+(qlonm/60.0))
*      if (cew.eq.'E') then
*        qlon=-qlon
*      end if
      qlon=(real(ilond)+(qlonm/60.0))
      if (cew.eq.'W') then
        qlon=-qlon
      end if
      aspect=cos(qlat/degrad)
      if (sez.eq.0.) sez=1.
      terr=-9                  ! set parameters not given in input file
      rms=-9
      nppick=-9
      nspick=-9
      evtype='L'
      magtype='X'
      locqual='X'

c choose a new source location and velocity model for each trial 
      qdep2(1)=qdep
      qlon2(1)=qlon
      qlat2(1)=qlat
      index(1)=1
      do nm=2,nmc
        call RAN_NORM(val)
*@        qdep2(nm)=qdep
        qdep2(nm)=qdep+sez*val    ! randomly perturbed source depth
        if (qdep2(nm).lt.0.0) qdep2(nm)=1e-6
* Randomly perturb the lat and long too
        call RAN_NORM(val)
*@        qlat2(nm)=qlat
        qlat2(nm)=qlat+seh*val/1.414/111.2    
        call RAN_NORM(val)
*@        qlon2(nm)=qlon
        qlon2(nm)=qlon+seh*val/1.414/111.2/aspect    
        index(nm)=mod(nm,ntab)+1  ! index used to choose velocity model
      end do      
c         print*,'1 ',val,qdep2,qdep
c read in polarities       ! SCEDC format - ** YOU MAY NEED TO CHANGE THE INPUT FORMAT **
      k=1
130   continue
        read (11,135,end=140) sname(k),snet(k),scomp(k),
     &       pickonset,pickpol
135     format (a4,1x,a2,2x,a3,1x,a1,1x,a1)
        if (sname(k).eq.'    ')  goto 140 ! end of data for this event
        call GETSTAT_TRI(stfile,sname(k),scomp(k),snet(k),
     &               flat,flon,felv)   ! SCSN station information - ** YOU MAY NEED TO USE YOUR OWN SUBROUTINE **
        if (flat.eq.999.) go to 130
*        dx=(flon-qlon)*111.2*aspect
*        dy=(flat-qlat)*111.2
*        range=sqrt(dx**2+dy**2)
*        qazi=90.-atan2(dy,dx)*degrad
*        if (qazi.lt.0.) qazi=qazi+360.
        if (pickpol.eq.'U'.or.
     &                    pickpol.eq.'u'.or.
     &                    pickpol.eq.'+') then
          p_pol(k)=1
        else if (pickpol.eq.'D'.or.
     &                    pickpol.eq.'d'.or.
     &                    pickpol.eq.'-') then
          p_pol(k)=-1
        else
          goto 130
        end if
        if ((pickonset.ne.'I').and.               
     &      (pickonset.ne.'i')) then
          q_pol(k) = 1
        else
          q_pol(k) = 0
        end if
        call CHECK_POL(plfile,sname(k),iyr,imon,idy,ihr,ispol)  
                          ! SCSN station polarity reversal information - ** YOU MAY NEED TO USE YOUR OWN SUBROUTINE **        
        p_pol(k)=p_pol(k)*ispol
        sp_ratio(k)=0.
        do 105 nm=1,nmc  ! find azimuth and takeoff angle for each trial
          dx=(flon-qlon2(nm))*111.2*aspect
          dy=(flat-qlat2(nm))*111.2
          range=sqrt(dx**2+dy**2)
          qazi=90.-atan2(dy,dx)*degrad
          if (qazi.lt.0.) qazi=qazi+360.
          p_azi_mc(k,nm)=qazi
          call GET_TTS(index(nm),range,qdep2(nm),
     &                 p_the_mc(k,nm),iflag)
*          if (nm.eq.1) then
*            print*,'1 ',index(nm),range,qdep2(nm),p_the_mc(k,nm),iflag
*          end if
105     continue
        k=k+1
      goto 130
140   continue
      nppl=k-1

c read in amplitude ratios - find those for corresponding event ID   !  ** YOU MAY NEED TO CHANGE THE INPUT FORMAT **
      open (12,file=ampfile,status='old')
20    read (12,*,end=41) icusp2,nin
c      print*,'1 ',icusp2,icusp,nin
      if (icusp2.ne.icusp) then
        do i=1,nin
          read (12,*,end=41)
        end do
        goto 20
      end if
30    do 40 i=1,nin
        read (12,35,end=41) sname(k),scomp(k),snet(k),
     &     qns1,qns2,qpamp,qsamp
35      format (a4,1x,a3,1x,a2,17x,f10.3,1x,f10.3,
     &           1x,f10.3,1x,f10.3)
        call GETSTAT_TRI(stfile,sname(k),scomp(k),snet(k),
     &               flat,flon,felv)   ! SCSN station information - ** YOU MAY NEED TO USE YOUR OWN SUBROUTINE **
*            print*,'1 ',sname(1),scomp(1),snet(1),qns1,qns2,qpamp,qsamp
        if (flat.eq.999.) go to 40
*        dx=(flon-qlon)*111.2*aspect
*        dy=(flat-qlat)*111.2
*        range=sqrt(dx**2+dy**2)
*        qazi=90.-atan2(dy,dx)*degrad
*        if (qazi.lt.0.) qazi=qazi+360.
        call GET_COR(corfile,sname(k),scomp(k),snet(k),qcor)
        s2n1=abs(qpamp)/qns1
        s2n2=qsamp/qns2
        spin=qsamp/abs(qpamp)
c        print*,'2 ',spin,s2n1,s2n2,ratmin
        if (qcor.eq.-999.) goto 40
        if (qpamp.eq.0.) goto 40
        if ((s2n1.lt.ratmin).or.(s2n2.lt.ratmin)) goto 40
        sp_ratio(k)=log10(spin)-qcor
c        print*,'2 ',sp_ratio(k)
        p_pol(k)=0
        q_pol(k)=0
        do nm=1,nmc  ! find azimuth and takeoff angle for each trial
          dx=(flon-qlon2(nm))*111.2*aspect
          dy=(flat-qlat2(nm))*111.2
          range=sqrt(dx**2+dy**2)
          qazi=90.-atan2(dy,dx)*degrad
          if (qazi.lt.0.) qazi=qazi+360.
          p_azi_mc(k,nm)=qazi
          call GET_TTS(index(nm),range,qdep2(nm),
     &                 p_the_mc(k,nm),iflag)
        end do
        k=k+1
40    continue
41    close(12)
      npol=k-1
      nspr=npol-nppl
      print*,'1 ',npol,nppl,nspr
cc view polarity data
c      do k=1,npol
c        print *,k,p_azi_mc(k,1),p_the_mc(k,1),p_pol(k),sp_ratio(k)
c      end do
      
      do 256 k=1,npol 
          write(16,*) icusp,k,' ',sname(k),p_azi_mc(k,1),p_the_mc(k,1),
     &                  p_pol(k),sp_ratio(k)
256   continue

      if (nppl.lt.1) then
        print *,'*** warning - no p-wave polarity data for event',
     &            icusp
      end if
      if (nspr.lt.1) then
        print *,'*** warning - no s/p amplitude ratios for event',
     &            icusp
      end if

c determine maximum azimuthal and takeoff gap in polarity observations
c and stop if either gap is too big
      call GET_GAP(npol,p_azi_mc,p_the_mc,magap,mpgap)
      if ((magap.gt.max_agap).or.(mpgap.gt.max_pgap)) then
        if (magap.gt.maxagap) then
          print*,'Azimuthal gap too big - Skipping event ',icusp
        end if
        if (magap.gt.maxpgap) then
          print*,'Take-off gap too big - Skipping event ',icusp
        end if
        str_avg(1)=999
        dip_avg(1)=99
        rak_avg(1)=999
        var_est(1,1)=99
        var_est(2,1)=99
        mfrac(1)=0.99
        qual(1)='E'
        prob(1)=0.0
        nout2=0
        nmult=0
        goto 400
      end if
      
*      print *,'EVID =',icusp, '  NDATA =',npol,nppl,nspr

*      nmis0max=max(nint(nppl*badfrac),2) 
	  nmis0max=0                   
      nmis1max=0
*      qmismax=max(nspr*qbadfac,2.0)                    
*      qextra=max(nspr*qbadfac*0.5,2.0)
      qmismax=0                 
      qextra=0
      print*,'nmis0max qmismax ',nmis0max,qmismax
      

      call FOCALAMP_MC(p_azi_mc,p_the_mc,sp_ratio,p_pol,q_pol,
     &    npol,nmc,
     &    dang2,nmax0,nextra0,nmis0max,nextra1,nmis1max,qextra,
     &    qmismax,nf2,strike2,dip2,rake2,f1norm,f2norm)

*      print*,'nf2 ',nf2
      nout2=min(nmax0,nf2)  ! number mechs returned from sub
      nout1=min(maxout,nf2)  ! number mechs to return
*      print*,'nout1 nout2 ',nout1,nout2
      
c find the probable mechanism from the set of acceptable solutions          
      call MECH_PROB(nout2,f1norm,f2norm,cangle,prob_max,nmult,
     &        str_avg,dip_avg,rak_avg,prob,var_est)           

      do 390 imult=1,nmult
      
      var_avg(imult)=(var_est(1,imult)+var_est(2,imult))/2.
      print *,'cid = ',icusp,imult,'  mech = ',
     &          str_avg(imult),dip_avg(imult),rak_avg(imult)

c find misfit for prefered solution
      call GET_MISF_AMP(npol,p_azi_mc,p_the_mc,sp_ratio,
     &      p_pol,str_avg(imult),dip_avg(imult),rak_avg(imult),
     &      mfrac(imult),mavg(imult),stdr(imult),
     &      polpred(1,imult),sppred(1,imult))
      
c solution quality rating, completely ad-hoc - make up your own!
      if ((prob(imult).gt.0.8).and.(var_avg(imult).le.25)) then
        qual(imult)='A'
      else if ((prob(imult).gt.0.6).and.(var_avg(imult).le.35)) then
        qual(imult)='B'
      else if ((prob(imult).gt.0.5).and.(var_avg(imult).le.45)) then
        qual(imult)='C'
      else
        qual(imult)='D'
      end if

390   continue

400   continue
      if (nmult.eq.0) then
        goto 120
      end if
       
c output prefered mechanism  ** YOU MAY WISH TO CHANGE THE OUTPUT FORMAT **
      do i=1,nmult
      write (13,411) icusp,iyr,imon,idy,ihr,imn,qsec,evtype,
     &   qmag,magtype,qlat,qlon,qdep,locqual,rms,seh,sez,terr,
     &   nppick+nspick,nppick,nspick,
     &   nint(str_avg(i)),nint(dip_avg(i)),nint(rak_avg(i)),
     &   nint(var_est(1,i)),nint(var_est(2,i)),nppl,nint(mfrac(i)*100.),
     &   qual(i),nint(100*prob(i)),nint(100*stdr(i)),nspr,
     &   nint(mavg(i)*100.),nmult,magap,mpgap
      end do
411   format(i16,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3,1x,a1,1x,
     &  f5.3,1x,a1,1x,f9.5,1x,f10.5,1x,f7.3,1x,a1,1x,f7.3,1x,f7.3,
     &  1x,f7.3,1x,f7.3,3x,i4,1x,i4,1x,i4,1x,i4,1x,i3,1x,i4,3x,i2,
     &  1x,i2,1x,i3,1x,i2,1x,a1,1x,i3,1x,i2,1x,i3,1x,i3,1x,i1,1x,
     &  i3,1x,i2)

* Add more detailed output as per hash_driver2
c output set of acceptable mechanisms  ** YOU MAY WISH TO CHANGE THE OUTPUT FORMAT **      
      write (14,412) iyr,imon,idy,ihr,imn,qsec,qmag,
     &   qlat,qlon,qdep,sez,seh,npol,nout2,icusp,
     &   str_avg(1),dip_avg(1),rak_avg(1),var_est(1,1),var_est(2,1),
     &   mfrac(1),qual(1),prob(1),stdr(1)
412   format (i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3,2x,f3.1,1x,f9.4,1x,
     &   f10.4,1x,f6.2,1x,f8.4,1x,f8.4,1x,i5,1x,i5,1x,i16,1x,f7.1,1x,
     &   f6.1,1x,f7.1,1x,f6.1,1x,f6.1,1x,f7.3,2x,a1,1x,f7.3,1x,f4.2)
      do 500 ic=1,nout1
        write (14,550) strike2(ic),dip2(ic),rake2(ic),f1norm(1,ic),
     &      f1norm(2,ic),f1norm(3,ic),f2norm(1,ic),f2norm(2,ic),
     &      f2norm(3,ic)
500   continue
550   format (5x,3f9.2,6f9.4)

* Add observations and predictions 
c      write (15,421) icusp, npol,nppl, nspr, nmult
c421   format(i16,1x,i3,1x,i3,1x,i3,1x,i1)
c      do i = 1,nmult
c        write(15,422) sname(nm),p_azi_mc(nm,1), p_the_mc(nm,1),
c    &     q_pol(nm),p_pol(nm),polpred(nm,i),sp_ratio(nm),
c     &     sppred(nm,i)
c      end do
c      end do
c422   format (a4,1x,f7.1,1x,f6.1,1x,i2,1x,i2,1x,i2,1x,f7.2,1x,f7.2)

c output azimuth and take-off angle and polarity data % MZ added 040422
c      do 255 k=1,nppl 
c          write(15,*) icusp,' ',sname(k),p_azi_mc(k,1),p_the_mc(k,1),
c     &                  p_pol(k)
c255   continue

C WW copied
* Add observations and predictions 
      write (15,421) icusp, npol,nppl, nspr, nmult
421   format(i16,1x,i3,1x,i3,1x,i3,1x,i1)
      do i = 1,nmult
      do nm = 1,npol
        write(15,422) sname(nm),p_azi_mc(nm,1), p_the_mc(nm,1),
     &     p_pol(nm),polpred(nm,i),sp_ratio(nm),sppred(nm,i)
      end do
      end do
422   format (a4,1x,f7.1,1x,f6.1,1x,i2,1x,i2,1x,f7.2,1x,f7.2)
* wsdw end _________________________

      goto 120
      
505   continue
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      stop
      end
