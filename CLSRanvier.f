C This program calculates the time evolution of the membrane potential
C at a node of Ranvier with left-shifted NaV currents.
C shft = the left shift of the affected channels in mV
C LS = the fraction of affected channels at a node
      program singlenodeHH
      implicit none
C Initializes variables
      integer imax,zi
      real*8 dt,tmax,V,t,RTF,m,h,n
      real*8 mLS,hLS,dVdt,Ii
      real*8 Ki,Nai,Ko,Nao
      real*8 INat,Inal,Ikl,Ikn,Ikp,Inapmp,Ileak,ILS
      real*8 LS,shft
      common/imax/imax
      common/I/INat,Inal,Ikl,Ikn,Ikp,Inapmp,Ileak,ILS,Ii
      common/RTF/RTF
      common/LS/LS
      common/shft/shft
      common/dt/dt

C Initializes the output files
      open(20,file='node1.dat')
      open(21,file='node2.dat')
      open(22,file='node3.dat')
      open(23,file='node4.dat')
      open(24,file='node5.dat')
      open(25,file='node6.dat')
      open(26,file='node7.dat')
      open(27,file='node8.dat')
      open(28,file='node9.dat')
      write(20,*) 't',' V',' It'
      write(21,*) 'm',' h',' n'
      write(22,*) ' mLS',' hLS'
      write(23,*) 'INat',' INal',' IKl'
      write(24,*) 'Ikn',' dVdt'
      write(25,*) ' Ileak',' Inapump',' Ikpump'
      write(26,*) 'Nai',' Ki'
      write(27,*) 'Nao',' Ko'
      write(28,*) 'Ena',' Ek'


C the do loop (using 4 as its end marker) is controled by the dummy-variable zi.
C it determines the phases of the simulation.
C phase 1 (zi = 1) is usually an equilibration phase
C imax = number of steps between output to file
C tmax = maximal time value for a phase (ms)
C Ii = injected current
C LS = relative amount of left-shifted channels
C shft = left-shift of channels (in mV)

      do 4 zi=1,3
      if (zi.eq.1) then
            imax=1000
            tmax=200.0d0
            Ii=0.0d0
            LS=0.00d0
            shft=0.0d0
            call init(t,m,h,n,mLS,hLS,V,Nao,Nai,Ko,Ki)
         elseif (zi.eq.2) then
            imax=1000
            tmax=400.0d0
            Ii=-12.0d0
         else
            imax=1000
            tmax=900.0d0
            Ii=0.0d0
      endif


C this while loop controls the progress of the simulation,
C time evolves as long as time is smaller than the set value of tmax
      dowhile (t.lt.tmax)

      
C propagate fait avancer de imax pas de longueur dt (variable?)
C the propagate subroutine calculates the time evolution of NaV channel
C parameters (m, h as well as mLS and hLS for shifted channels), n for
C K channels, Nai and Ki (interior concentration), Nao and Ko (exterior concentration),
C time t and membrane potential V.
C the propagate subroutine advances for a number of imax steps.
      call propagate(m,h,mLS,hLS,n,V,Nai,Nao,Ki,Ko,t)

C the update subroutine recalculates the values of currents with
C changed values of parameters.
C update is called here simply to output correct current values.
      call update(m,h,mLS,hLS,n,V,Nai,Nao,Ki,Ko)


C Output to files.
      write(20,*) t,V,INat+Inal+Ikl+Ikn+Ikp
     * +Inapmp+Ii+Ileak+ILS
      write(21,*) m,h,n
      write(22,*) mLS,hLS
      write(23,*) ILS+Inat,Inal,Ikl
      write(24,*) Ikn,dVdt(INat+Inal+Ikl+
     * Ikn+Ikp+Inapmp+Ii+Ileak+ILS)
      write(25,*) Ileak,Inapmp,Ikp
      write(26,*) Nai,Ki
      write(27,*) Nao,Ko
      write(28,*) RTF*dlog(Nao/Nai),RTF*dlog(Ko/Ki)

      enddo
4     continue
      stop
      end



C The propagate subroutine uses an RK4 scheme. The update subroutine
C is called each step to use correct values for currents.
      subroutine propagate(m,h,mLS,hLS,n,V,Nai,Nao,Ki,Ko,t)
      implicit none
      integer i,imax,count
      real*8 m,h,mLS,hLS,n,V,Nai,Nao,Ki,Ko,t
      real*8 mk1,mk2,mk3,mk4,hk1,hk2,hk3,hk4
      real*8 mLSk1,mLSk2,mLSk3,mLSk4,hLSk1,hLSk2,hLSk3,hLSk4
      real*8 nk1,nk2,nk3,nk4,vk1,vk2,vk3,vk4
      real*8 naik1,naik2,naik3,naik4,naok1,naok2,naok3,naok4
      real*8 kik1,kik2,kik3,kik4,kok1,kok2,kok3,kok4
      real*8 dmdt,dhdt,dndt,dCdt,dVdt
      real*8 Inat,Inal,Ikl,Ikn,Ikp,Inapmp,Ii,Ileak,Ils
      real*8 shft,dt,voli,volo,test1,test2,test3
      common/imax/imax
      common/vol/volo
      common/volume/voli
      common/shft/shft
      common/dt/dt
      common/I/INat,Inal,Ikl,Ikn,Ikp,Inapmp,Ileak,Ils,Ii
      count=0
      do 1 i=1,imax
2     continue
      mk1=dt*dmdt(m,V)
      hk1=dt*dhdt(h,V)
      mLSk1=dt*dmdt(mLS,V+shft)
      hLSk1=dt*dhdt(hLS,V+shft)
      nk1=dt*dndt(n,V)
      call update(m,h,mLS,hLS,n,V,Nai,Nao,Ki,Ko)
      Vk1=dt*dVdt(INat+Inal+Ikl+Ikn+Ikp+Inapmp+Ii+Ileak+ILS)
      naik1=dt*dCdt(Inat+Inal+Inapmp+ILS)
      kik1=dt*dCdt(Ikl+Ikn+Ikp)
      Naok1=-1.0d0*naik1*voli/volo
      kok1=-1.0d0*kik1*voli/volo


      mk2=dt*dmdt(m+0.5d0*mk1,V+0.5d0*Vk1)
      hk2=dt*dhdt(h+0.5d0*hk1,V+0.5d0*Vk1)
      mLSk2=dt*dmdt(mLS+0.5d0*mLSk1,V+0.5d0*Vk1+shft)
      hLSk2=dt*dhdt(hLS+0.5d0*hLSk1,V+0.5d0*Vk1+shft)
      nk2=dt*dndt(n+0.5d0*nk1,V+0.5d0*Vk1)
      call update(m+0.5d0*mk1,h+0.5d0*hk1,mLS+0.5d0*mLSk1,
     * hLS+0.5d0*hLSk1,n+0.5d0*nk1,V+0.5d0*Vk1,
     *Nai+0.5d0*Naik1,Nao+0.5d0*Naok1,Ki+0.5d0*Kik1,Ko+0.5d0*kok1)
      Vk2=dt*dVdt(INat+Inal+Ikl+Ikn+Ikp+Inapmp+Ii+Ileak+ILS)
      naik2=dt*dCdt(Inat+Inal+Inapmp+ILS)
      kik2=dt*dCdt(Ikl+Ikn+Ikp)
      Naok2=-1.0d0*naik2*voli/volo
      kok2=-1.0d0*kik2*voli/volo

      mk3=dt*dmdt(m+0.5d0*mk2,V+0.5d0*Vk2)
      hk3=dt*dhdt(h+0.5d0*hk2,V+0.5d0*Vk2)
      mLSk3=dt*dmdt(mLS+0.5d0*mLSk2,V+0.5d0*Vk2+shft)
      hLSk3=dt*dhdt(hLS+0.5d0*hLSk2,V+0.5d0*Vk2+shft)
      nk3=dt*dndt(n+0.5d0*nk2,V+0.5d0*Vk2)
      call update(m+0.5d0*mk2,h+0.5d0*hk2,mLS+0.5d0*mLSk2,
     * hLS+0.5d0*hLSk2,n+0.5d0*nk2,V+0.5d0*Vk2,
     *nai+0.5d0*naik2,nao+0.5d0*naok2,Ki+0.5d0*kik2,Ko+0.5d0*kok2)
      Vk3=dt*dVdt(INat+Inal+Ikl+Ikn+Ikp+Inapmp+Ii+Ileak+ILS)
      naik3=dt*dCdt(Inat+Inal+Inapmp+ILS)
      kik3=dt*dCdt(Ikl+Ikn+Ikp)
      Naok3=-1.0d0*naik3*voli/volo
      kok3=-1.0d0*kik3*voli/volo

      mk4=dt*dmdt(m+mk3,V+Vk3)
      hk4=dt*dhdt(h+hk3,V+Vk3)
      mLSk4=dt*dmdt(mLS+mLSk3,V+Vk3+shft)
      hLSk4=dt*dhdt(hLS+hLSk3,V+Vk3+shft)
      nk4=dt*dndt(n+nk3,V+Vk3)
      call update(m+mk3,h+hk3,mLS+mLSk3,hLS+hLSk3,n+nk3,
     *V+dt*Vk3,nai+naik3,nao+naok3,Ki+kik3,Ko+kok3)
      Vk4=dt*dVdt(INat+Inal+Ikl+Ikn+Ikp+Inapmp+Ii
     * +Ileak+ILS)
      naik4=dt*dCdt(Inat+Inal+Inapmp+ILS)
      kik4=dt*dCdt(Ikl+Ikn+Ikp)
      Naok4=-1.0d0*naik4*voli/volo
      kok4=-1.0d0*kik4*voli/volo


C this section ensures a rapid simulation by allowing the timestep
C to vary is values are changing rapidly or slowly.
      test1=dabs((mk1+2.0d0*mk2+2.0d0*mk3+mk4)/6.0d0/m)
      test2=dabs((Vk1+2.0d0*Vk2+2.0d0*Vk3+Vk4)/6.0d0/V)
      test3=dabs(((mLSk1+2.0d0*mLSk2+2.0d0*mLSk3+mLSk4)/6.0d0))/mLS
      if (((test1.ge.1.0d-3).or.
     *(test2.ge.1.0d-3).or.(test3.ge.1.0d-3)).and.(dt.ge.1.0d-4)) then
      if (count.ge.2) goto 3
      dt=dt/10.0d0
      count=count+1
      goto 2
      elseif ((test1.le.1.0d-5).and.(test2.le.1.0d-5).and.
     * (test3.le.1.0d-5).and.(dt.lt.1.0d-2)) then
      if (count.ge.2) goto 3
      dt=dt*10.0d0
      count=count+1
      goto 2
      endif
3     continue
      count=0

C The values of the parameters are calculated with the correct RK4 weight.
      m=m+(mk1+2.0d0*mk2+2.0d0*mk3+mk4)/6.0d0
      h=h+(hk1+2.0d0*hk2+2.0d0*hk3+hk4)/6.0d0
      mLS=mLS+(mLSk1+2.0d0*mLSk2+2.0d0*mLSk3+mLSk4)/6.0d0
      hLS=hLS+(hLSk1+2.0d0*hLSk2+2.0d0*hLSk3+hLSk4)/6.0d0
      n=n+(nk1+2.0d0*nk2+2.0d0*nk3+nk4)/6.0d0
      V=V+(Vk1+2.0d0*Vk2+2.0d0*Vk3+Vk4)/6.0d0
      Nai=Nai+(naik1+2.0d0*naik2+2.0d0*naik3+naik4)/6.0d0
      Ki=Ki+(Kik1+2.0d0*Kik2+2.0d0*Kik3+Kik4)/6.0d0
      Nao=Nao+(naok1+2.0d0*naok2+2.0d0*naok3+naok4)/6.0d0
      Ko=Ko+(Kok1+2.0d0*Kok2+2.0d0*Kok3+Kok4)/6.0d0
      t=t+dt
1     continue

      return
      end
      
      
      
C The update subroutine calculates the values for the currents (which are common variables)
C given a set of parameters.
      subroutine update(m,h,mLS,hLS,n,V,Nai,Nao,Ki,Ko)
      implicit none
      real*8 m,h,mLS,hLS,n,V,Nai,Nao,Ki,Ko
      real*8 INat,Inal,Ikl,Ikn,Ikp,Inapmp,Ileak,Ii,ILS
      real*8 gna,gnal,gkl,gkn,Imaxp
      real*8 Apump,gleak,LS,Ek,Ena,RTF,Eleak
      common/up1/gna,gnal,gkl,gkn,Imaxp,gleak
      common/I/INat,Inal,Ikl,Ikn,Ikp,Inapmp,Ileak,ILS,Ii
      common/RTF/RTF
      common/LS/LS
      common/Eleak/Eleak
      Ena=RTF*dlog(Nao/Nai)
      INat=(1.0d0-LS)*gna*m**3*h*(V-Ena)
      ILS=LS*gna*mLS**3*hLS*(V-Ena)
      Inal=gnal*(V-Ena)

      Ek=RTF*dlog(Ko/Ki)
      Ikl=gkl*(V-Ek)
      Ikn=gkn*n**4*(V-Ek)

      Ileak=gleak*(V-Eleak)

      Apump=Imaxp*(1.0d0+3.5d0/Ko)**(-2)*(1.0d0+10.0d0/Nai)**(-3)
      Ikp=-2.0d0*Apump
      Inapmp=3.0d0*Apump

      return
      end


C this function calculates the rate of change of membrane potential.
      function dVdt(Itot)
      implicit none
      real*8 Itot,C,dVdt
      common/Cap/C
      dVdt=-1.0d0*Itot/C
      return
      end

C this function calculates the rate of change of one species of ion.
      function dCdt(Itot)
      implicit none
      real*8 dCdt,Itot,F,voli,surf
      common/F/F
      common/volume/voli
      common/surf/surf
      dCdt=-1.0d-6*Itot*surf/F/voli
      return
      end

C rate of change for the n variable.
      function dndt(n,V)
      implicit none
      real*8 dndt,n,V,alpha,beta
      alpha=0.01d0*(V+55.0d0)/(1.0d0-dexp(-1.0d0*(V+55.0d0)/10.0d0))
      beta=0.1250d0*dexp(-1.0d0*(V+65.0d0)/80.0d0)
      dndt=alpha*(1.0d0-n)-beta*n
      return
      end


C rate of change of the m variable.
      function dmdt(m,V)
      implicit none
      real*8 dmdt,m,V,alpha,beta
      alpha=0.1d0*(V+40.0d0)/(1.0d0-dexp(-1.0d0*(V+40.0d0)/10.0d0))
      beta=4.0d0*dexp(-1.0d0*(V+65.0d0)/18.0d0)
      dmdt=alpha*(1.0d0-m)-beta*m
      return
      end

C rate of change of the h variable.
      function dhdt(h,V)
      implicit none
      real*8 dhdt,h,V,alpha,beta
      alpha=0.07d0*dexp(-1.0d0*(V+65.0d0)/20.0d0)
      beta=1.0d0/(1.0d0+dexp(-1.0d0*(V+35.0d0)/10.0d0))
      dhdt=alpha*(1.0d0-h)-beta*h
      return
      end


C this subroutine simply sets ths initial values of most parameters.
      subroutine init(t,m,h,n,mLS,hLS,V,Nao,Nai,Ko,Ki)
      implicit none
      real*8 t,m,h,n,mLS,hLS,V,Nao,Nai,Ko,Ki
      real*8 gna,gkn,Imaxp,gnal,gkl,gleak,voli,volo,C,F
      real*8 Temp,dt,RTF,Eleak,surf,R
      real*8 INat,Inal,Ikl,Ikn,Ikp,Inapmp,Ileak,ILS,Ii
      common/cap/C
      common/F/F
      common/RTF/RTF
      common/dt/dt
      common/vol/volo
      common/volume/voli
      common/up1/gna,gnal,gkl,gkn,Imaxp,gleak
      common/Eleak/Eleak
      common/surf/surf
      common/I/INat,Inal,Ikl,Ikn,Ikp,Inapmp,Ileak,ILS,Ii
      F=96485.3399d0
      R=8.314472d0
      Temp=293.150d0
      RTF=R*temp/F*1000.0d0
      
      gna=120.0d0
      gkn=36.0d0
C      gleak=0.25d0
      gleak=0.50d0
      Eleak=-59.9d0
      C=1.0d0
      
      gkl=0.1d0
      Imaxp=1.0d0
      gnal=1.0d0
      
C Volume in m3
C surface in cm2
C      voli=1.0d100
      voli=3.0d-15
C      volo=1.0d100
      volo=voli
C      volo=0.15*1.0d-15
C      volo=0.150d0*voli
      surf=6.0d-8

      Ko=6.0d0
      Nao=154.0d0
      Ki=150.0d0
      Nai=20.0d0
      dt=1.0d-3
      t=0.0d0

      m=0.095d0
      h=0.414d0
      n=0.398d0
      V=-59.9d0
      mls=m
      hls=h

      write(*,*) m,h,m**3*h
      write(*,*) 'mLS=',mls,', hls=',hls,' mls**3hls=',mls**3*hls
      call update(m,h,mLS,hLS,n,V,Nai,Nao,Ki,Ko)
      gnal=(-3.0d0/2.0d0*(Ikl+Ikn)-(Inat+ILS))/Inal
      call update(m,h,mLS,hLS,n,V,Nai,Nao,Ki,Ko)
      Imaxp=-1.0d0*(Inat+Inal+ILS)/Inapmp
      call update(m,h,mLS,hLS,n,V,Nai,Nao,Ki,Ko)
      write(*,*) 'Inattotal: ',Inat+Inal+ILS,'=',Inat,'+',Inal,'+',ILS
      write(*,*) 'Iktot: ',Ikl+Ikn,'=',Ikl,'+',Ikn
      write(*,*) 'Ina/Iktot: ',(Inat+Inal+ILS)/(Ikl+Ikn)
      write(*,*) 'A= ', (-3.0d0/2.0d0*(Ikl+Ikn)-(Inat+ILS))/Inal
      write(*,*) 'Inapump: ',Inapmp,' Ikp: ',Ikp
      write(*,*) 'Inapmp/Ikp: ',Inapmp/Ikp
      write(*,*) '-Inat/Inp: ', -1.0d0*(Inat+Inal+ILS)/Inapmp
      write(*,*) 'Ileak= ',Ileak,' Itotal=',INat+Inal+Ikl+
     * Ikn+Ikp+Inapmp+Ii+ILS
      write(*,*) 'Ek= ',RTF*dlog(Ko/Ki),
     * ' Ena= ',RTF*dlog(Nao/Nai)
      write(*,*) 'gnal:',gnal,' Imaxp:',Imaxp
C      pause
      return
      end
