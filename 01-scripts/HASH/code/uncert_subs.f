c subroutine MECH_PROB determines the average focal mechanism of a set
c   of mechanisms, after removing outliers, and returns the probability
c   of the mechanism being within a given cutoff angle of the result,
c   also checks outliers for multiple solutions and returns any
c
c  Inputs:  nf     =  number of fault planes
c           norm1(3,nf) = normal to fault plane
c           norm2(3,nf) = slip vector
c           cangle  =  cutoff angle 
c           prob_max = cutoff percent for multiples
c  Output:  nsltn = number of solutions, up to 5
c           str_avg(5)   = average strike
c           dip_avg(5)   = average dip
c           rak_avg(5)   = average rake
c           prob(5)      = percent of mechs within cutoff angle
c                          of average mechanism
c           rms_diff(2,5)  = RMS angular difference of all planes to average 
c                          plane (1=fault plane, 2=auxiliary plane)
c

      subroutine MECH_PROB(nf,norm1in,norm2in,cangle,prob_max,nsltn,
     &             str_avg,dip_avg,rak_avg,prob,rms_diff)

      include 'param.inc'
      
      integer nf
      real str_avg(5),dip_avg(5),rak_avg(5),rota(nmax0)
      real norm1(3,nmax0),norm2(3,nmax0),temp1(3),temp2(3)
      real norm1in(3,nf),norm2in(3,nf),ln_norm1,ln_norm2
      real norm1_avg(3),norm2_avg(3),slipol,rms_diff(2,5)
      real stv(2),udv(3),dd_rad,di_rad,a,b,prob(5)

      pi=3.1415927
      degrad=180./3.1415927
      
c if there is only one mechanism, return that mechanism
   
      if (nf.le.1) then
        do i=1,3
          norm1_avg(i)=norm1in(i,1)
          norm2_avg(i)=norm2in(i,1)
        end do
        call fpcoor(str_avg(1),dip_avg(1),rak_avg(1),
     &              norm1_avg,norm2_avg,2)
        prob(1)=1.
        rms_diff(1,1)=0.
        rms_diff(2,1)=0.
        nsltn=1
        return
      end if
      
c otherwise find the prefered mechanism and any multiples

      do j=1,nf
        do i=1,3
          norm1(i,j)=norm1in(i,j)
          norm2(i,j)=norm2in(i,j)
        end do
      end do
      print*,'8.2 nf',nf
      nfault=nf
      nc=nf
      do 380 imult=1,5
        if (nc.lt.1) goto 385

c find the average mechanism by repeatedly throwing out the mechanism
c with the largest angular difference from the average, and finding the average of 
c the remaining mechanisms - stop when all mechanisms are within cangle of average
   
      do 77 count=1,nf 
        call MECH_AVG(nc,norm1,norm2,norm1_avg,norm2_avg) 
        do 50 i=1,nc      
          d11=norm1(1,i)*norm1_avg(1)+norm1(2,i)*norm1_avg(2)+
     &                              norm1(3,i)*norm1_avg(3)
          d12=norm1(1,i)*norm2_avg(1)+norm1(2,i)*norm2_avg(2)+
     &                              norm1(3,i)*norm2_avg(3)
          d21=norm2(1,i)*norm1_avg(1)+norm2(2,i)*norm1_avg(2)+
     &                              norm2(3,i)*norm1_avg(3)
          d22=norm2(1,i)*norm2_avg(1)+norm2(2,i)*norm2_avg(2)+
     &                              norm2(3,i)*norm2_avg(3)
          if ((d11*d22.ge.d21*d12).or.(abs(d11).gt.0.99999).or.
     &                              (abs(d22).gt.0.99999)) then
            if (abs(d11).gt.0.000001) then
              isgn=sign(1.,d11)
            else
              isgn=sign(1.,d22)
            end if
            do j=1,3
              temp1(j)=isgn*norm1(j,i)
              temp2(j)=isgn*norm2(j,i)
            end do
          else
            if (abs(d21).gt.0.000001) then
              isgn=sign(1.,d21)
            else
              isgn=sign(1.,d12)
            end if
            do j=1,3
              temp1(j)=isgn*norm2(j,i)
              temp2(j)=isgn*norm1(j,i)
            end do
          end if
          call MECH_ROT(temp1,norm1_avg,temp2,norm2_avg,rota(i))
50      continue
        maxrot=0.
        do i=1,nc
          if (abs(rota(i)).gt.maxrot) then
            maxrot=abs(rota(i))
            imax=i
          end if
        end do
        if (maxrot.le.cangle) goto 78
        nc=nc-1
        do i=1,3
          temp1(i)=norm1(i,imax)
          temp2(i)=norm2(i,imax)
        end do
        do j=imax,nc            
          do i=1,3
            norm1(i,j)=norm1(i,j+1)
            norm2(i,j)=norm2(i,j+1)
          end do
        end do
        do i=1,3
          norm1(i,nc+1)=temp1(i)
          norm2(i,nc+1)=temp2(i)
        end do
77    continue

78    continue
      a=nc
      b=nfault
      prob(imult)=a/b 
      print*,'8.4 a b prob',a,b,prob

      if ((imult.gt.1).and.(prob(imult).lt.prob_max)) goto 385

      do j=1,nf-nc   ! set up for next round
        do i=1,3
          norm1(i,j)=norm1(i,j+nc)
          norm2(i,j)=norm2(i,j+nc)
        end do
      end do
      nc=nf-nc 
      nf=nc

c determine the RMS observed angular difference between the average 
c normal vectors and the normal vectors of each mechanism

      rms_diff(1,imult)=0.
      rms_diff(2,imult)=0.
      do 80 i=1,nfault
        d11=norm1in(1,i)*norm1_avg(1)+norm1in(2,i)*norm1_avg(2)+
     &                              norm1in(3,i)*norm1_avg(3)
        d12=norm1in(1,i)*norm2_avg(1)+norm1in(2,i)*norm2_avg(2)+
     &                              norm1in(3,i)*norm2_avg(3)
        d21=norm2in(1,i)*norm1_avg(1)+norm2in(2,i)*norm1_avg(2)+
     &                              norm2in(3,i)*norm1_avg(3)
        d22=norm2in(1,i)*norm2_avg(1)+norm2in(2,i)*norm2_avg(2)+
     &                              norm2in(3,i)*norm2_avg(3)
        if (d11.ge.1.) d11=1.
        if (d11.le.-1.) d11=-1.
        if (d12.ge.1.) d12=1.
        if (d12.le.-1.) d12=-1.
        if (d21.ge.1.) d21=1.
        if (d21.le.-1.) d21=-1.
        if (d22.ge.1.) d22=1.
        if (d22.le.-1.) d22=-1.
        a11=acos(abs(d11))
        a22=acos(sign(1.,d11)*d22)
        a21=acos(abs(d21))
        a12=acos(sign(1.,d21)*d12)
        if ((d11*d22.ge.d21*d12).or.(abs(d11).gt.0.99999).or.
     &                              (abs(d22).gt.0.99999)) then
          rms_diff(1,imult)=rms_diff(1,imult)+a11*a11
          rms_diff(2,imult)=rms_diff(2,imult)+a22*a22
        else
          rms_diff(1,imult)=rms_diff(1,imult)+a21*a21
          rms_diff(2,imult)=rms_diff(2,imult)+a12*a12
        end if
80    continue
      rms_diff(1,imult)=degrad*sqrt(rms_diff(1,imult)/nfault)
      rms_diff(2,imult)=degrad*sqrt(rms_diff(2,imult)/nfault)

c the average normal vectors may not be exactly orthogonal (although
c usually they are very close) - find the misfit from orthogonal and 
c adjust the vectors to make them orthogonal - adjust the more poorly 
c constrained plane more
 
      if ((rms_diff(1,imult)+rms_diff(2,imult)).lt.0.0001) goto 120

      maxmisf=0.01
      fract1=rms_diff(1,imult)/(rms_diff(1,imult)+rms_diff(2,imult))
90    do 115 count=1,100  
        dot1=norm1_avg(1)*norm2_avg(1)+norm1_avg(2)
     &     *norm2_avg(2)+norm1_avg(3)*norm2_avg(3)
        misf=90.-acos(dot1)*degrad
        if (abs(misf).le.maxmisf) goto 120
        theta1=misf*fract1/degrad
        theta2=misf*(1.-fract1)/degrad
        do j=1,3
          temp=norm1_avg(j)
          norm1_avg(j)=norm1_avg(j)-norm2_avg(j)*sin(theta1)
          norm2_avg(j)=norm2_avg(j)-temp*sin(theta2)
        end do
        ln_norm1=0
        ln_norm2=0
        do j=1,3
          ln_norm1=ln_norm1+norm1_avg(j)*norm1_avg(j)
          ln_norm2=ln_norm2+norm2_avg(j)*norm2_avg(j)
        end do
        ln_norm1=sqrt(ln_norm1)
        ln_norm2=sqrt(ln_norm2)
        do i=1,3
          norm1_avg(i)=norm1_avg(i)/ln_norm1
          norm2_avg(i)=norm2_avg(i)/ln_norm2
        end do
115   continue

120   continue  

      call fpcoor(str_avg(imult),dip_avg(imult),rak_avg(imult),
     &            norm1_avg,norm2_avg,2)

380   continue

385   nsltn=imult-1  ! only use ones with high enough probability
          
      return 
      end

c ------------------------------------------------------------ c

c subroutine MECH_AVG determines the average focal mechanism of a set
c   of mechanisms
c
c  Inputs:  nf     =  number of fault planes
c           norm1(3,nf) = normal to fault plane
c           norm2(3,nf) = slip vector
c  Output:  norm1_avg(3)    = normal to avg plane 1
c           norm2_avg(3)    = normal to avg plane 2
c
c    Written  10/4/2000 by Jeanne Hardebeck                              
c    Modified 5/14/2001 by Jeanne Hardebeck                              
c
      subroutine MECH_AVG(nf,norm1,norm2,norm1_avg,norm2_avg)
           
      real dot1,fract1
      real misf,maxmisf,avang1,avang2
      real norm1(3,nf),norm2(3,nf)
      real norm1_avg(3),norm2_avg(3),ln_norm1,ln_norm2
      real theta1,theta2
      integer nf

      pi=3.1415927
      degrad=180./3.1415927
      
c if there is only one mechanism, return that mechanism
   
      if (nf.le.1) then
        do i=1,3
          norm1_avg(i)=norm1(i,1)
          norm2_avg(i)=norm2(i,1)
        end do
        goto 120
      end if
            
c find the average normal vector for each plane - determine which
c nodal plane of each event corresponds to which running sum by 
c computing the dot product of both nodal planes with the fault  
c plane of the first mechanism and taking the one which is closer

      do j=1,3
        norm1_avg(j)=norm1(j,1)
        norm2_avg(j)=norm2(j,1)
      end do
      do 50 i=2,nf
        d11=norm1(1,i)*norm1(1,1)+norm1(2,i)*norm1(2,1)+
     &                              norm1(3,i)*norm1(3,1)
        d12=norm1(1,i)*norm2(1,1)+norm1(2,i)*norm2(2,1)+
     &                              norm1(3,i)*norm2(3,1)
        d21=norm2(1,i)*norm1(1,1)+norm2(2,i)*norm1(2,1)+
     &                              norm2(3,i)*norm1(3,1)
        d22=norm2(1,i)*norm2(1,1)+norm2(2,i)*norm2(2,1)+
     &                              norm2(3,i)*norm2(3,1)
        if ((d11*d22.ge.d21*d12).or.(abs(d11).gt.0.99999).or.
     &                              (abs(d22).gt.0.99999)) then
          if (abs(d11).gt.0.000001) then
            isgn=sign(1.,d11)
          else
            isgn=sign(1.,d22)
          end if
          do j=1,3
            norm1_avg(j)=norm1_avg(j)+isgn*norm1(j,i)
            norm2_avg(j)=norm2_avg(j)+isgn*norm2(j,i)
          end do
        else
          if (abs(d21).gt.0.000001) then
            isgn=sign(1.,d21)
          else
            isgn=sign(1.,d12)
          end if
          do j=1,3
            norm1_avg(j)=norm1_avg(j)+isgn*norm2(j,i)
            norm2_avg(j)=norm2_avg(j)+isgn*norm1(j,i)
          end do
        end if
50    continue
 
      ln_norm1=0
      ln_norm2=0
      do j=1,3
        ln_norm1=ln_norm1+norm1_avg(j)*norm1_avg(j)
        ln_norm2=ln_norm2+norm2_avg(j)*norm2_avg(j)
      end do
      ln_norm1=sqrt(ln_norm1)
      ln_norm2=sqrt(ln_norm2)
      do i=1,3
        norm1_avg(i)=norm1_avg(i)/ln_norm1
        norm2_avg(i)=norm2_avg(i)/ln_norm2
      end do

c determine the RMS observed angular difference between the average 
c normal vectors and the normal vectors of each mechanism

      avang1=0.
      avang2=0.

      do 80 i=1,nf
        d11=norm1(1,i)*norm1_avg(1)+norm1(2,i)*norm1_avg(2)+
     &                              norm1(3,i)*norm1_avg(3)
        d12=norm1(1,i)*norm2_avg(1)+norm1(2,i)*norm2_avg(2)+
     &                              norm1(3,i)*norm2_avg(3)
        d21=norm2(1,i)*norm1_avg(1)+norm2(2,i)*norm1_avg(2)+
     &                              norm2(3,i)*norm1_avg(3)
        d22=norm2(1,i)*norm2_avg(1)+norm2(2,i)*norm2_avg(2)+
     &                              norm2(3,i)*norm2_avg(3)
        if (d11.ge.1.) d11=1.
        if (d11.le.-1.) d11=-1.
        if (d12.ge.1.) d12=1.
        if (d12.le.-1.) d12=-1.
        if (d21.ge.1.) d21=1.
        if (d21.le.-1.) d21=-1.
        if (d22.ge.1.) d22=1.
        if (d22.le.-1.) d22=-1.
        a11=acos(abs(d11))
        a22=acos(sign(1.,d11)*d22)
        a21=acos(abs(d21))
        a12=acos(sign(1.,d21)*d12)
        if (d11*d22.ge.d21*d12) then
          avang1=avang1+a11*a11
          avang2=avang2+a22*a22
        else
          avang1=avang1+a21*a21
          avang2=avang2+a12*a12
        end if
80    continue
      avang1=sqrt(avang1/nf)
      avang2=sqrt(avang2/nf)

c the average normal vectors may not be exactly orthogonal (although
c usually they are very close) - find the misfit from orthogonal and 
c adjust the vectors to make them orthogonal - adjust the more poorly 
c constrained plane more
 
      if ((avang1+avang2).lt.0.0001) goto 120

      maxmisf=0.01
      fract1=avang1/(avang1+avang2)
90    do 115 count=1,100  
        dot1=norm1_avg(1)*norm2_avg(1)+norm1_avg(2)
     &     *norm2_avg(2)+norm1_avg(3)*norm2_avg(3)
        misf=90.-acos(dot1)*degrad
        if (abs(misf).le.maxmisf) goto 120
        theta1=misf*fract1/degrad
        theta2=misf*(1.-fract1)/degrad
        do j=1,3
          temp=norm1_avg(j)
          norm1_avg(j)=norm1_avg(j)-norm2_avg(j)*sin(theta1)
          norm2_avg(j)=norm2_avg(j)-temp*sin(theta2)
        end do
        ln_norm1=0
        ln_norm2=0
        do j=1,3
          ln_norm1=ln_norm1+norm1_avg(j)*norm1_avg(j)
          ln_norm2=ln_norm2+norm2_avg(j)*norm2_avg(j)
        end do
        ln_norm1=sqrt(ln_norm1)
        ln_norm2=sqrt(ln_norm2)
        do i=1,3
          norm1_avg(i)=norm1_avg(i)/ln_norm1
          norm2_avg(i)=norm2_avg(i)/ln_norm2
        end do
115   continue

120   continue      
      return 
      end
      

c ------------------------------------------------------------ c


c subroutine MECH_ROT finds the minimum rotation angle between two
c mechanisms
c
c  Inputs:  P1(3) = normal to fault plane 1  OR  P-axis 1
c           T1(3) = slip vector 1            OR  T-axis 1
c           P2(3) = normal to fault plane 2  OR  P-axis 2
c           T2(3) = slip vector 2            OR  T-axis 2
c  Output:  rota  = rotation angle
c
      subroutine MECH_ROT(P1,P2,T1,T2,rota)
      
      real P1(3),P2(3),T1(3),T2(3),B1(3),B2(3)
      real rota,phi(3),n(3,3),scale(3),R(3),qdot(3)
      real theta(3),n1(3),n2(3),phi1
      
      pi=3.1415927
      degrad=180./3.1415927

      call cross(P1,T1,B1)
      call cross(P2,T2,B2)

      phi1=P1(1)*P2(1)+P1(2)*P2(2)+P1(3)*P2(3)
      phi(1)=acos(phi1)
      phi1=T1(1)*T2(1)+T1(2)*T2(2)+T1(3)*T2(3)
      phi(2)=acos(phi1)
      phi1=B1(1)*B2(1)+B1(2)*B2(2)+B1(3)*B2(3)
      phi(3)=acos(phi1)
      
c if the mechanisms are very close, return 0
      if ((phi(1).lt.1e-4).and.(phi(2).lt.1e-4).and.
     &    (phi(3).lt.1e-4)) then
        rota=0.0
c if one vector is the same, it is the rotation axis
      else if (phi(1).lt.1e-4) then
        rota=degrad*phi(2)
      else if (phi(2).lt.1e-4) then
        rota=degrad*phi(3)
      else if (phi(3).lt.1e-4) then
        rota=degrad*phi(1)
      else
c find difference vectors - the rotation axis must be orthogonal
c to all three of these vectors
        do i=1,3
          n(i,1)=P1(i)-P2(i)
          n(i,2)=T1(i)-T2(i)
          n(i,3)=B1(i)-B2(i)
        end do
        do j=1,3
          scale(j)=sqrt(n(1,j)*n(1,j)+n(2,j)*n(2,j)+n(3,j)*n(3,j))
          do i=1,3
            n(i,j)=n(i,j)/scale(j)
          end do
        end do
        qdot(3)=n(1,1)*n(1,2)+n(2,1)*n(2,2)+n(3,1)*n(3,2)
        qdot(2)=n(1,1)*n(1,3)+n(2,1)*n(2,3)+n(3,1)*n(3,3)
        qdot(1)=n(1,2)*n(1,3)+n(2,2)*n(2,3)+n(3,2)*n(3,3)
c use the two largest difference vectors, as long as they aren't orthogonal 
        iout=0
        do i=1,3
          if (qdot(i).gt.0.9999) iout=i
        end do
        if (iout.eq.0) then
          qmins=10000.
          do i=1,3
            if (scale(i).lt.qmins) then
              qmins=scale(i)
              iout=i
            end if
          end do
        end if
        k=1
        do j=1,3
          if (j.ne.iout) then
            if (k.eq.1) then
              do i=1,3
                n1(i)=n(i,j)
              end do
              k=2
            else
              do i=1,3
                n2(i)=n(i,j)
              end do
            end if
          end if
        end do
c  find rotation axis by taking cross product
         call CROSS(n1,n2,R)
         scaleR=sqrt(R(1)*R(1)+R(2)*R(2)+R(3)*R(3))
         do i=1,3
           R(i)=R(i)/scaleR
         end do
c find rotation using axis furthest from rotation axis
         theta(1)=acos(P1(1)*R(1)+P1(2)*R(2)+P1(3)*R(3))
         theta(2)=acos(T1(1)*R(1)+T1(2)*R(2)+T1(3)*R(3))
         theta(3)=acos(B1(1)*R(1)+B1(2)*R(2)+B1(3)*R(3))
         qmindif=1000.
         do i=1,3
           if (abs(theta(i)-pi/2.0).lt.qmindif) then
             qmindif=abs(theta(i)-pi/2.0)
             iuse=i
           end if
         end do
         rota=(cos(phi(iuse))-cos(theta(iuse))*cos(theta(iuse)))/
     &          (sin(theta(iuse))*sin(theta(iuse)))
         if (rota.gt.1.0) then
           rota=1.0
         end if
         if (rota.lt.-1.0) then
           rota=-1.0
         end if
         rota=degrad*acos(rota)
       end if
       return
       end      
