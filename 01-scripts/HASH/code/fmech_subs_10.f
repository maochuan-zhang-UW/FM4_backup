c subroutine FOCALMC performs grid search to find acceptable focal mechanisms,
c                    for multiple trials of ray azimuths and takeoff angles.
c                    Acceptable mechanisms are those with less than "ntotal"
c                    misfit polarities, or the minimum plus "nextra" if this
c                    is greater.
c
c  Inputs:  
c           p_azi_mc(npsta,nmc)  =  azimuth to station from event (deg. E of N)
c           p_the_mc(npsta,nmc)  =  takeoff angle (from vert, up=0, <90 upgoing, >90 downgoing)
c           p_pol(npsta)  =  first motion, 1=up, -1=down
c           p_qual(npsta) =  quality, 0=impulsive, 1=emergent
c           npsta  =  number of first motions
c           nmc    =  number of trials
c           dang   =  desired angle spacing for grid search
c           maxout =  maximum number of fault planes to return:
c                     if more are found, a random selection will be returned
c           nextra =  number of additional misfits allowed above minimum
c           ntotal =  total number of allowed misfits
c  Outputs: 
c           nf     =  number of fault planes found
c           strike(min(maxout,nf)) = strike
c           dip(min(maxout,nf))    = dip
c           rake(min(maxout,nf))   = rake
c           faults(3,min(maxout,nf)) = fault normal vector
c           slips(3,min(maxout,nf))  = slip vector
c
c
      subroutine FOCALMC(p_azi_mc,p_the_mc,p_pol,p_qual,npsta,nmc,
     &    dang,maxout,nextra,ntotal,nf,strike,dip,
     &    rake,faults,slips)
     
      include 'param.inc'
      include 'rot_10.inc'
     
c input and output arrays
      dimension p_azi_mc(npick0,nmc0),p_the_mc(npick0,nmc0)
      integer p_pol(npsta),p_qual(npsta)
      real p_a1(npick0),p_a2(npick0),p_a3(npick0)
      real faultnorm(3),slip(3),faults(3,maxout),slips(3,maxout)
      real strike(maxout),dip(maxout),rake(maxout)
      save dangold,nrot,b1,b2,b3

c coordinate transformation arrays
      real b1(3,ncoor),bb1(3)
      real b2(3,ncoor),bb2(3)
      real b3(3,ncoor),bb3(3)

c fit arrays
      integer fit(2,ncoor),nmiss01min(0:npick0)
      integer irotgood(ncoor2),irotgood2(ncoor2)
      real fran

      pi=3.1415927
      degrad=180./pi

c Set up array with direction cosines for all coordinate transformations
      if (dang.eq.dangold) go to 8
      print*,'dang ',dang
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
8     continue
      print*,'nrot', nrot

c loop over multiple trials
      nfault=0
      do 430 im=1,nmc 
      print*,'nmc im ',nmc,im

c  Convert data to Cartesian coordinates
      do 40 i=1,npsta
         call TO_CAR(p_the_mc(i,im),p_azi_mc(i,im),1.,
     &               p_a1(i),p_a2(i),p_a3(i))
40    continue

c  find misfit for each solution and minimum misfit
         nmiss0min=999
         nmissmin=999
c eeeh         do 390 i=0,nsta
         do 390 i=0,npsta
390         nmiss01min(i)=999
         do 420 irot=1,nrot  
            nmiss=0
            nmiss0=0           
            do 400 ista=1,npsta
               p_b1= b1(1,irot)*p_a1(ista)
     &              +b1(2,irot)*p_a2(ista)
     &              +b1(3,irot)*p_a3(ista) 
               p_b3= b3(1,irot)*p_a1(ista)
     &              +b3(2,irot)*p_a2(ista)
     &              +b3(3,irot)*p_a3(ista) 
               prod=p_b1*p_b3
               ipol=-1
               if (prod.gt.0.) ipol=1    ! predicted polarization
               if (ipol.ne.p_pol(ista)) then
                  nmiss=nmiss+1                       
                  if (p_qual(ista).eq.0) nmiss0=nmiss0+1
               end if
400         continue
            fit(1,irot)=nmiss0           ! misfit impulsive polarities
            fit(2,irot)=nmiss            ! total misfit polarities
            if (nmiss0.lt.nmiss0min) nmiss0min=nmiss0
            if (nmiss.lt.nmissmin) nmissmin=nmiss
            if (nmiss.lt.nmiss01min(nmiss0)) then
               nmiss01min(nmiss0)=nmiss
            end if
420      continue

c choose fit criteria
         if (nmiss0min.eq.0) then 
            nmiss0max=ntotal
            nmissmax=ntotal
         else
            nmiss0max=ntotal
            nmissmax=npsta
         end if
         if (nmiss0max.lt.nmiss0min+nextra) then
            nmiss0max=nmiss0min+nextra
         end if
         if (nmissmax.lt.nmiss01min(nmiss0min)+nextra) then
            nmissmax=nmiss01min(nmiss0min)+nextra
         end if

c loop over rotations - find those meeting fit criteria
         do 440 irot=1,nrot        
            nmiss0=fit(1,irot)
            nmiss=fit(2,irot)
            if (nmiss0.gt.nmiss0max.or.nmiss.gt.nmissmax) go to 440
            nfault=nfault+1
            print*,'nfault irot ',nfault,irot
            print*,'  nmiss nmiss0 nmissmax nmiss0max '
     c              ,nmiss, nmiss0, nmissmax, nmiss0max
            if (nfault.le.ncoor2) then
              irotgood(nfault)=irot
              irotgood2(nfault)=irot
            end if
440     continue

430     continue

c  Select output solutions        
        nreturn=nfault
        if (nfault.le.maxout) go to 445
        nreturn=maxout
        print*,'nreturn maxout',nreturn,maxout
        j=0
        do 441 i=1,99999
           fran=rand(0)
           iscr=nint(fran*float(nfault)+0.5)
           if (iscr.lt.1) iscr=1
           if (iscr.gt.nfault) iscr=nfault
           if (irotgood(iscr).lt.0) go to 441
           j=j+1
           irotgood2(j)=irotgood(iscr)
           irotgood(iscr)=-irotgood(iscr)
           if (j.eq.maxout) go to 445
441     continue

445     nf=0
        print*,'$$ nreturn ', nreturn
        do 450 i=1,nreturn
            irot=irotgood2(i)
            nf=nf+1
            print*,i,nf
            faultnorm(1)=b3(1,irot)
            faultnorm(2)=b3(2,irot)
            faultnorm(3)=b3(3,irot)
            slip(1)=b1(1,irot)
            slip(2)=b1(2,irot)
            slip(3)=b1(3,irot)
            do 447 m=1,3
              faults(m,nf)=faultnorm(m)
              slips(m,nf)=slip(m)
447         continue
            call FPCOOR(s1,d1,r1,faultnorm,slip,2)
            strike(nf)=s1
            dip(nf)=d1
            rake(nf)=r1
*            if (nf.eq.nf0) go to 452
450      continue 
452      continue
*         nf=nfault            
         print*,'$ nf ', nf
*         nf=nreturn            
         print*,'$$ nf ', nf
c
      return 
      end
