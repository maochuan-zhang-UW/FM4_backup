c subroutine FOCALAMP_MC performs grid search to find focal mechanisms, using 
c            both P-polarity and S/P amplitude ratio information
*
* Modified by William Wilcock to weight all solutions that fit data equally
* and to include P-wave quality (0=impulsive; 1=emergent)
*
c  Inputs:  
c           p_azi_mc(npsta,nmc)  =  azimuth to station from event (deg. E of N)
c           p_the_mc(npsta,nmc)  =  takeoff angle (from vert, <90 upgoing, >90 downgoing)
c           sp_amp(npsta)  =  amplitude ratios
c           p_pol(nspta)   =  P polarities
*           p_qual(npsta)  = P polarity qualities (0 = impulsive; 1 = emergent)
c           npsta  =  number of observations
c           nmc    =  number of trials
c           dang   =  desired angle spacing for grid search
c           maxout =  maximum number of fault planes to return:
c                     if more are found, a random selection will be returned
*           nextra0 =  number of impulsive polarity additional misfits 
*                      allowed above minimum
c           ntotal0 =  total number of allowed impulsive polarity misfits
*           nextra1 =  number of emergent polarity misfits 
*                      allowed above minimum
*           ntotal1 =  total number of allowed emergent polarity misfits
c           qextra =  additional amplitude misfit allowed above minimum
c           qtotal =  total allowed amplitude misfit
c  Outputs: 
c           nf     =  number of fault planes found
c           strike(min(maxout,nf)) = strike
c           dip(min(maxout,nf))    = dip
c           rake(min(maxout,nf))   = rake
c           faults(3,min(maxout,nf)) = fault normal vector
c           slips(3,min(maxout,nf))  = slip vector
c

      subroutine FOCALAMP_MC(p_azi_mc,p_the_mc,sp_amp,p_pol,p_qual,
     &     npsta,nmc,dang,maxout,nextra0,ntotal0,nextra1,ntotal1,
     &     qextra,qtotal,nf,strike,dip,rake,faults,slips)
     
      include 'param.inc'
      include 'rot.inc'
      parameter (ntab=180) 

c input and output arrays
      dimension p_azi_mc(npick0,nmc0),p_the_mc(npick0,nmc0)
      real sp_amp(npsta)
      real p_a1(npick0),p_a2(npick0),p_a3(npick0)
      real faultnorm(3),slip(3),faults(3,nmax0),slips(3,nmax0)
      real strike(nmax0),dip(nmax0),rake(nmax0)
      integer p_pol(npsta),p_qual(npsta)
      save dangold,nrot,b1,b2,b3    
      save amptable,phitable,thetable

c coordinate transformation arrays
      real b1(3,ncoor),bb1(3)
      real b2(3,ncoor),bb2(3)
      real b3(3,ncoor),bb3(3)
c P and S amplitude arrays
      real amptable(2,ntab,2*ntab)
      real phitable(2*ntab+1,2*ntab+1)
      real thetable(2*ntab+1)
c misfit arrays
      real qmis(ncoor)
      integer nmis0(ncoor),nmis1(ncoor)
      integer nmis1min(npick0)
      integer irotgood(ncoor*nmax0)
      
      pi=3.1415927
      degrad=180./pi

      if (maxout.gt.nmax0) then
         maxout=nmax0
      end if

c Set up array with direction cosines for all coordinate transformations
      if (dang.eq.dangold) go to 8
      irot=0
      do 5 the=0.,90.1,dang
         rthe=the/degrad
         costhe=cos(rthe)
         sinthe=sin(rthe)
         fnumang=360./dang
         numphi=nint(fnumang*sin(rthe))
         if (numphi.ne.0) then
            dphi=360./float(numphi)
         else
            dphi=10000.
         end if
         do 4 phi=0.,359.9,dphi
            rphi=phi/degrad
            cosphi=cos(rphi)
            sinphi=sin(rphi)
            bb3(3)=costhe
            bb3(1)=sinthe*cosphi
            bb3(2)=sinthe*sinphi
            bb1(3)=-sinthe
            bb1(1)=costhe*cosphi
            bb1(2)=costhe*sinphi
            call CROSS(bb3,bb1,bb2)
            do 3 zeta=0.,179.9,dang
               rzeta=zeta/degrad
               coszeta=cos(rzeta)
               sinzeta=sin(rzeta)
               irot=irot+1
               if (irot.gt.ncoor) then
                  print *,'***FOCAL error: # of rotations too big'
                  return
               end if
               b3(3,irot)=bb3(3)
               b3(1,irot)=bb3(1)
               b3(2,irot)=bb3(2)
               b1(1,irot)=bb1(1)*coszeta+bb2(1)*sinzeta
               b1(2,irot)=bb1(2)*coszeta+bb2(2)*sinzeta                
               b1(3,irot)=bb1(3)*coszeta+bb2(3)*sinzeta
               b2(1,irot)=bb2(1)*coszeta-bb1(1)*sinzeta
               b2(2,irot)=bb2(2)*coszeta-bb1(2)*sinzeta                
               b2(3,irot)=bb2(3)*coszeta-bb1(3)*sinzeta
3           continue
4        continue
5     continue
      nrot=irot
      dangold=dang
      
      astep=1./real(ntab)
      do i=1,2*ntab+1
        bbb3=-1.+real(i-1)*astep
        thetable(i)=acos(bbb3)
        do j=1,2*ntab+1
          bbb1=-1.+real(j-1)*astep
          phitable(i,j)=atan2(bbb3,bbb1)
          if (phitable(i,j).lt.0.) then
            phitable(i,j)=phitable(i,j)+2.*pi
          end if
        end do
      end do

      do i=1,2*ntab
        phi=real(i-1)*pi*astep
        do j=1,ntab
          theta=real(j-1)*pi*astep
          amptable(1,j,i)=abs(sin(2*theta)*cos(phi))                
          s1=cos(2*theta)*cos(phi)  
          s2=-cos(theta)*sin(phi)
          amptable(2,j,i)=sqrt(s1*s1+s2*s2)
        end do
      end do

8     continue

c loop over multiple trials
      nfault = 0
      do 430 im=1,nmc 

c  Convert data to Cartesian coordinates//transforms spherical co-ordinates to cartesian
      do i=1,npsta
         call TO_CAR(p_the_mc(i,im),p_azi_mc(i,im),1.,
     &               p_a1(i),p_a2(i),p_a3(i))
      end do

c  find misfit for each solution and minimum misfit
         nmis0min=1e5
         do i=1,npsta+1
           nmis1min(i)=1e5
         end do
         qmis0min=1.0e5
         do 420 irot=1,nrot  
           qmis(irot)=0.
           nmis0(irot)=0
           nmis1(irot)=0
           do 400 ista=1,npsta
             p_b1= b1(1,irot)*p_a1(ista)
     &              +b1(2,irot)*p_a2(ista)
     &              +b1(3,irot)*p_a3(ista) 
             p_b3= b3(1,irot)*p_a1(ista)
     &              +b3(2,irot)*p_a2(ista)
     &              +b3(3,irot)*p_a3(ista) 
             if (sp_amp(ista).ne.0.) then
               p_proj1=p_a1(ista)-p_b3*b3(1,irot)
               p_proj2=p_a2(ista)-p_b3*b3(2,irot)
               p_proj3=p_a3(ista)-p_b3*b3(3,irot)
               plen=sqrt(p_proj1*p_proj1+p_proj2*p_proj2+
     &                    p_proj3*p_proj3)
               p_proj1=p_proj1/plen
               p_proj2=p_proj2/plen
               p_proj3=p_proj3/plen
               pp_b1=b1(1,irot)*p_proj1+b1(2,irot)*p_proj2
     &                +b1(3,irot)*p_proj3
               pp_b2=b2(1,irot)*p_proj1+b2(2,irot)*p_proj2
     &              +b2(3,irot)*p_proj3
               i=nint((p_b3+1.)/astep)+1
               theta=thetable(i)
               i=nint((pp_b2+1.)/astep)+1
               j=nint((pp_b1+1.)/astep)+1
               phi=phitable(i,j)
               i=nint(phi/(pi*astep))+1
               if (i.gt.2*ntab) i=1
               j=nint(theta/(pi*astep))+1
               if (j.gt.ntab) j=1
               p_amp=amptable(1,j,i)
               s_amp=amptable(2,j,i)
               if (p_amp.eq.0.0) then
                 sp_ratio=4.0
               else if (s_amp.eq.0.0) then
                 sp_ratio=-2.0
               else
* 1.7^3 = 4.9
                 sp_ratio=log10(5.088*s_amp/p_amp)
               end if
               qmis(irot)=qmis(irot)+abs(sp_amp(ista)-sp_ratio)
             end if
             if (p_pol(ista).ne.0) then
               prod=p_b1*p_b3
               ipol=-1
               if (prod.gt.0.) ipol=1 
               if (ipol.ne.p_pol(ista)) then
                 if (p_qual(ista).eq.0) then
                   nmis0(irot)=nmis0(irot)+1                       
                 else
                   nmis1(irot)=nmis1(irot)+1                       
                 end if
               end if
             end if
400         continue
            if (nmis0(irot).lt.nmis0min) then
              nmis0min=nmis0(irot)
              if (nmis1(irot).lt.nmis1min(nmis0(irot)+1)) then
                nmis1min(nmis0(irot)+1)=nmis1(irot)
              end if
            end if
            if (qmis(irot).lt.qmis0min) qmis0min=qmis(irot)
420      continue
 
         nmis0max=ntotal0
         nmis1max=ntotal1
         if (nmis0max.lt.nmis0min+nextra0) then
            nmis0max=nmis0min+nextra0
         end if
         if (nmis1max.lt.nmis1min(nmis0min+1)+nextra1) then
            nmis1max=nmis1min(nmis0min+1)+nextra1
         end if
         qmis0max=qtotal
         if (qmis0max.lt.qmis0min+qextra) then
            qmis0max=qmis0min+qextra
         end if

c loop over rotations - find those meeting fit criteria
         nstuck=0
425      nadd=0
         do irot=1,nrot        
           if (nstuck.eq.10) then
             print*,'fmamp_subs.f stuck in loop'
             stop
           end if
      ! print*,'9 nmis0(irot) nmis1(irot)',nmis0(irot),nmis1(irot)
            if ((nmis0(irot).le.nmis0max).and. 
     &          (nmis1(irot).le.nmis1max).and. 
     &            (qmis(irot).le.qmis0max)) then
              nfault = nfault+1
              irotgood(nfault)=irot
              nadd=nadd+1
            end if
         end do
        print*,'9.1 nadd nmc',nadd,nmc
         if (nadd.eq.0) then  ! if there are none that meet criteria
           nstuck = nstuck+1
           qmis0min=1.0e5     ! loosen the amplitude criteria
           do irot=1,nrot        
             if ((nmis0(irot).le.nmis0max).and.
     &           (nmis1(irot).le.nmis1max).and. 
     &           (qmis(irot).lt.qmis0min)) then
               qmis0min=qmis(irot)
             end if
           end do
           qmis0max=qtotal
           if (qmis0max.lt.qmis0min+qextra) then
              qmis0max=qmis0min+qextra
           end if
           goto 425
         end if

430     continue

c  Select output solutions  
        nf=0      
        if (nfault.le.maxout) then
          do i=1,nfault
            irot=irotgood(i)
            nf=nf+1
            faultnorm(1)=b3(1,irot)
            faultnorm(2)=b3(2,irot)
            faultnorm(3)=b3(3,irot)
            slip(1)=b1(1,irot)
            slip(2)=b1(2,irot)
            slip(3)=b1(3,irot)
            do m=1,3
              faults(m,nf)=faultnorm(m)
              slips(m,nf)=slip(m)
            end do
            call FPCOOR(s1,d1,r1,faultnorm,slip,2)
            strike(nf)=s1
            dip(nf)=d1
            rake(nf)=r1
          end do
        else
          do 441 i=1,99999
            fran=rand(0)
            iscr=nint(fran*float(nfault)+0.5)
            if (iscr.lt.1) iscr=1
            if (iscr.gt.nfault) iscr=nfault
            if (irotgood(iscr).le.0) goto 441
            irot=irotgood(iscr)
            irotgood(iscr)=-1
            nf=nf+1
            faultnorm(1)=b3(1,irot)
            faultnorm(2)=b3(2,irot)
            faultnorm(3)=b3(3,irot)
            slip(1)=b1(1,irot)
            slip(2)=b1(2,irot)
            slip(3)=b1(3,irot)
            do m=1,3
              faults(m,nf)=faultnorm(m)
              slips(m,nf)=slip(m)
            end do
            call FPCOOR(s1,d1,r1,faultnorm,slip,2)
            strike(nf)=s1
            dip(nf)=d1
            rake(nf)=r1
            if (nf.eq.maxout) go to 445
441       continue
445       continue
        end if  
         
      return 
      end

c ------------------------------------------------------------------- c
      
c GET_COR reads a file of station amplitude corrections
c
c   inputs:
c     stlfile - file with the stations and locations
c     snam - name of the station of interest, 4 characters
c     scom - station component, 3 characters
c     snet - network, 2 characters
c   outputs:
c     qcor - corrections to be subtracted from log(S/P)
c
c   input file format:
c     columns  format   variable
c     -------------------------
c     1-4        a4     station name
c     6-8        a3     station component (vertical)
c     10-11      a2     network code
c     13-19      f7.4   correction to be subtracted from log(S/P)
c
c
      subroutine GET_COR(stlfile,snam,scom,snet,qcor)
      parameter(nsta0=10000)
      character stlfile*100
      character*4 snam,stname(nsta0)
      character*3 scom,scom2,scompt(nsta0)
      character*2 snet,snetwk(nsta0)
      real corr_val(nsta0)
      logical firstcall
      save firstcall,stname,corr_val,nsta,scompt,snetwk
      data firstcall/.true./
      
c read in station list - in alphabetical order!
      if (firstcall) then
         firstcall=.false.
         open (19,file=stlfile)
         do i=1,nsta0
           read (19,11,end=12) stname(i),scompt(i),snetwk(i),
     &                              corr_val(i)
         end do
11       format (a4,1x,a3,1x,a2,1x,f7.4)
12       nsta=i-1
         close (19)
      end if  
      
      scom2=scom                             ! short-period stations are
      if (scom(1:1).eq."V") scom2(1:1)="E"   ! called both V and E     
      if (scom(1:1).eq."E") scom2(1:1)="V"           

c binary search for station name
      i1=1
      i2=nsta
      do it=1,30
         i=(i1+i2)/2
         if (snam.eq.stname(i)) then
           goto 40
         else if (i1.eq.i2) then
           goto 999
         else if (snam.lt.stname(i)) then
            i2=i-1
         else 
            i1=i+1
         end if
      end do
      print *,'station not found'
      goto 999
      
c search for proper component/network
40    i1=i
45    continue
        if (i1.gt.nsta) goto 50
        if (snam.ne.stname(i1)) goto 50
        if (scom(1:2).eq.scompt(i1)(1:2)) goto 900
        if (scom2(1:2).eq.scompt(i1)(1:2)) goto 900
        i1=i1+1
      goto 45
50    i1=i-1
55    continue
        if (i1.lt.1) goto 999
        if (snam.ne.stname(i1)) goto 999
        if (scom(1:2).eq.scompt(i1)(1:2)) goto 900
        if (scom2(1:2).eq.scompt(i1)(1:2)) goto 900
        i1=i1-1
      goto 55

900   qcor=corr_val(i1)
      return
999   print *,'GET_COR ***station not found ',snam,' ',scom,' ',snet,
     &  ' in file ',stlfile
      qcor=-999.
      return
      end

c ------------------------------------------------------------------- c
      

c subroutine GET_MISF_AMP finds the percent of misfit polarities and the
c                         average S/P ratio misfit for a given mechanism  
c    Inputs:    npol   = number of polarity observations
c               p_azi_mc(npol) = azimuths
c               p_the_mc(npol) = takeoff angles
c               sp_ratio(npol) = S/P ratio
c               p_pol(npol)  = polarity observations
c               str_avg,dip_avg,rak_avg = mechanism
c    Outputs:   mfrac = weighted fraction misfit polarities
c               mavg = average S/P misfit (log10)
c               stdr = station distribution ratio
* wsdw start _________________________
*               polpred(npol) = Predicted polarities
*               sppred(npol) = Predicted log10(S/P)
* wsdw end _________________________

      subroutine GET_MISF_AMP(npol,p_azi_mc,p_the_mc,sp_ratio,p_pol,
     &          str_avg,dip_avg,rak_avg,mfrac,mavg,stdr,
* wsdw start _________________________
     &          polpred,sppred)
* wsdw end _________________________

      dimension p_azi_mc(npol),p_the_mc(npol)
      real str_avg,dip_avg,rak_avg,M(3,3),a(3),b(3),sp_ratio(npol)
      real strike,dip,rake,mfrac,mavg,qcount,azi,toff,pol,wt,wo
      integer k,npol,p_pol(npol)
      real bb1(3),bb2(3),bb3(3)
* wsdw start _________________________
      integer polpred(npol)
      real sppred(npol)
* wsdw end _________________________
      
      rad=3.14159265/180.

      strike=str_avg*rad
      dip=dip_avg*rad
      rake=rak_avg*rad
      
      M(1,1)=-sin(dip)*cos(rake)*sin(2*strike)-sin(2*dip)*sin(rake)*
     & sin(strike)*sin(strike)
      M(2,2)=sin(dip)*cos(rake)*sin(2*strike)-sin(2*dip)*sin(rake)*
     & cos(strike)*cos(strike)
      M(3,3)=sin(2*dip)*sin(rake)
      M(1,2)=sin(dip)*cos(rake)*cos(2*strike)+0.5*sin(2*dip)*sin(rake)*
     & sin(2*strike)
      M(2,1)=M(1,2)
      M(1,3)=-cos(dip)*cos(rake)*cos(strike)-cos(2*dip)*sin(rake)*
     & sin(strike)
      M(3,1)=M(1,3)
      M(2,3)=-cos(dip)*cos(rake)*sin(strike)+cos(2*dip)*sin(rake)*
     & cos(strike)
      M(3,2)=M(2,3)
      call FPCOOR(strike,dip,rake,bb3,bb1,1)
      call CROSS(bb3,bb1,bb2)
      
      mfrac=0.
      qcount=0.
      stdr=0.
      scount=0.
      mavg=0.
      acount=0.
      
      do 600 k=1,npol
          call TO_CAR(p_the_mc(k),p_azi_mc(k),1.,p_a1,
     &                p_a2,p_a3)
          p_b1= bb1(1)*p_a1
     &              +bb1(2)*p_a2
     &              +bb1(3)*p_a3 
          p_b3= bb3(1)*p_a1
     &              +bb3(2)*p_a2
     &              +bb3(3)*p_a3
          p_proj1=p_a1-p_b3*bb3(1)
          p_proj2=p_a2-p_b3*bb3(2)
          p_proj3=p_a3-p_b3*bb3(3)
          plen=sqrt(p_proj1*p_proj1+p_proj2*p_proj2+
     &                    p_proj3*p_proj3)
          p_proj1=p_proj1/plen
          p_proj2=p_proj2/plen
          p_proj3=p_proj3/plen
          pp_b1=bb1(1)*p_proj1+bb1(2)*p_proj2
     &              +bb1(3)*p_proj3
          pp_b2=bb2(1)*p_proj1+bb2(2)*p_proj2
     &              +bb2(3)*p_proj3
          phi=atan2(pp_b2,pp_b1)
          theta=acos(p_b3)
          p_amp=abs(sin(2*theta)*cos(phi))     
          wt=sqrt(p_amp)
          if (p_pol(k).ne.0) then
            azi=rad*p_azi_mc(k)
            toff=rad*p_the_mc(k)        
            a(1)=sin(toff)*cos(azi)
            a(2)=sin(toff)*sin(azi)
            a(3)=-cos(toff)
            do in=1,3
              b(in)=0
              do jn=1,3
                 b(in)=b(in)+M(in,jn)*a(jn)
              end do
            end do
            if ((a(1)*b(1)+a(2)*b(2)+a(3)*b(3)).lt.0) then
              pol=-1
            else
             pol=1
            end if
* wsdw start _________________________
            polpred(k)=pol
* wsdw end _________________________
            if ((pol*p_pol(k)).lt.0) then
              mfrac=mfrac+wt
            end if
            qcount=qcount+wt
            stdr=stdr+wt
            scount=scount+1.0
* wsdw start _________________________
          else
            polpred(k)=0
* wsdw end _________________________
          end if
          if (sp_ratio(k).ne.0.) then
            s1=cos(2*theta)*cos(phi)  
            s2=-cos(theta)*sin(phi)
            s_amp=sqrt(s1*s1+s2*s2)
* wsdw start _________________________
* 1.7^3 = 4.9
            sp_rat=log10(5.088*s_amp/p_amp)               
            sppred(k) = sp_rat
* wsdw end _________________________
            mavg=mavg+abs(sp_ratio(k)-sp_rat)
            acount=acount+1.0
            stdr=stdr+wt
            scount=scount+1.0
* wsdw start _________________________
          else
            sppred(k)=0.0
* wsdw end _________________________
          end if
600    continue
       mfrac=mfrac/qcount
       if (qcount.eq.0.0) mfrac=0.0
       mavg=mavg/acount
       if (acount.eq.0.0) mavg=0.0
       stdr=stdr/scount
       if (scount.eq.0.0) stdr=0.0
       
       return 
       end
