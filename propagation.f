C this program calculates the time evolution of the membrane voltage over 
C a series of nodes connected by simple resistors. The simulation is
C carried 30 times (30 values of shift) for 6 different values of the fraction
C of affected channels (LS variable in this program). No
C stimulating current is used on the initial node, 
C so spontaneous firing of the nodes is studied.
C a single node (node 5) has shifted channels.
      program propagation
      implicit none
C variables are declared
      integer i,ii,spikeg,zi
      real*8 m(0:10),h(0:10),n(0:10),V(0:10),t,tmax,kappa,Iext,V0,temps
      real*8 kappa1,shft,LS,mLS,hLS
      common/C2/kappa,Iext
      common/LS/LS
      common/shft/shft
C files are initialized.
      open(11,file='propLS1.dat')
      open(12,file='propLS2.dat')
      open(13,file='propLS3.dat')
      open(14,file='propLS4.dat')
      open(15,file='propLS5.dat')
      open(16,file='propLS6.dat')

C zi is a dummy variable controlling the simulation, which is done over
C 6 different phases. Each phase uses a different value of LS (the ratio
C of shifted channels).
      do 123 zi=0,5
C text is wirtten in the appropriate file.
      write(zi+11,*) 'shft',' LS',' rate'
C the ratio of shifted channels is set.
      LS=0.2d0*zi

C this loop controls the simulation which is done 30 times with a fixed
C LS (previously set with the zi loop). The ii loop controls the value
C of the left-shift itself.
      do 911 ii=0,30
      shft=1.00d0*ii

C kappa1 is the interaction parameters between nodes.
      kappa1=0.14d0
C the initC and initVar subroutines initializes the parameters
      call initC(t,spikeg)
      call initVar(m,h,n,V,mls,hls)

C the i loop controls the time evolution of all simulations.
      do 99 i=1,2
C first, nodes are isolated (kappa = 0) and allowed to evolve without
C interactions or injected current for 100 ms.
      if (i.eq.1) then
      kappa=0.0d0
      Iext=0.0d0
      tmax=100.0d0
C then, interactions are turned on and the simulation continues until tmax.
      elseif (i.eq.2) then
      spikeg=0
      kappa=kappa1
      tmax=1100.0d0
      temps=tmax-t
      endif


C simulation continues until tmax.
      dowhile (t.lt.tmax)
C potential of the node 5 is recorded before the timestep
      V0=V(5)

C propagate updates 1 timestep.
      call propagate(m,h,n,V,t,mls,hls)

C if the potential of node 5 crosses 15 mV, a spike is recorded.
      if ((V0.lt.15.0d0).and.(V(5).ge.15.0d0)) then
      spikeg=spikeg+1
      endif


      enddo
99    continue

C number of spikes and parameters are written to file.
      write(zi+11,*) shft,LS,spikeg*1.0d0/temps
C999   continue
911   continue
123   continue
      stop
      end


C this subroutine propagates parameters 1 timestep using Euler's method.
      subroutine propagate(m,h,n,V,t,mls,hls)
      implicit none
      integer i,nnode
      real*8 m(0:10),h(0:10),n(0:10),V(0:10),t,inter(0:10)
      real*8 dt,dmdt,dhdt,dndt,dVdt,dV2dt,C,shft,mLS,hLS,V5
      common/C3/C,dt
      common/node/nnode
      common/shft/shft
C interaction currents are calculated.
      call interactions(V,inter)
C potential of node 5 is recorded.
      V5=V(5)
C all parameters are updated 1 timestep. Note that potential of node
C 5 is also updated here, but not correctly. It is corrected below.
      do 1 i=0,(Nnode-1)
      m(i)=dmdt(m(i),V(i))*dt+m(i)
      h(i)=dhdt(h(i),V(i))*dt+h(i)
      n(i)=dndt(n(i),V(i))*dt+n(i)
      V(i)=dt*(dVdt(m(i),h(i),n(i),V(i))+inter(i))/C+V(i)
1     continue
C left-shifted parameters are updated.
      mLS=dmdt(mLS,V(5)+shft)*dt+mLS
      hLS=dhdt(hLS,V(5)+shft)*dt+hLS
C potential of node 5 is correctly updated.
      V(5)=dt*(dV2dt(m(5),h(5),n(5),V5,mls,hls)+inter(5))/C+V5
      t=t+dt
      return
      end


C this subroutine calculates the interaction between nodes.
      subroutine interactions(V,inter)
      implicit none
      integer i,nnode
      real*8 V(0:10),inter(0:10),kappa,Iext
      common/C2/kappa,Iext
      common/node/nnode
      inter(0)=kappa*(V(1)-V(0))+Iext
      inter(9)=kappa*(V(8)-V(9))
      do 2 i=1,(nnode-2)
      inter(i)=kappa*(V(i-1)+V(i+1)-2.0d0*V(i))
2     continue
      return
      end
      
      function dmdt(m,V)
      implicit none
      real*8 dmdt,m,V,alpha,beta
      alpha=0.1d0*(V+40.0d0)/(1.0d0-dexp(-1.0d0*(V+40.0d0)/10.0d0))
      beta=4.0d0*dexp(-1.0d0*(V+65.0d0)/18.0d0)
      dmdt=alpha*(1.0d0-m)-beta*m
      return
      end
      
      function dhdt(h,V)
      implicit none
      real*8 dhdt,h,V,alpha,beta
      alpha=0.07d0*dexp(-1.0d0*(V+65.0d0)/20.0d0)
      beta=1.0d0/(1.0d0+dexp(-1.0d0*(V+35.0d0)/10.0d0))
      dhdt=alpha*(1.0d0-h)-beta*h
      return
      end
      
      function dndt(n,V)
      implicit none
      real*8 dndt,n,V,alpha,beta
      alpha=0.01d0*(V+55.0d0)/(1.0d0-dexp(-1.0d0*(V+55.0d0)/10.0d0))
      beta=0.1250d0*dexp(-1.0d0*(V+65.0d0)/80.0d0)
      dndt=alpha*(1.0d0-n)-beta*n
      return
      end
      
C rate of change for potential.
      function dVdt(m,h,n,V)
      implicit none
      real*8 dVdt,m,h,n,V,gna,Ena,gk,Ek,gleak,El
      common/C1/gna,Ena,gk,Ek,gleak,El
      dVdt=-gna*m**3*h*(V-Ena)-gk*n**4*(V-Ek)-gleak*(V-El)
      return
      end

C rate of change of potential in a node with shifted channels.
      function dV2dt(m,h,n,V,mls,hls)
      implicit none
      real*8 dV2dt,m,h,n,V,gna,Ena,gk,Ek,gleak,El,mLS,hLS,LS
      common/C1/gna,Ena,gk,Ek,gleak,El
      common/LS/LS
      dV2dt=-gna*m**3*h*(V-Ena)*(1.0d0-Ls)-gk*n**4*(V-Ek)-gleak*(V-El)
     *-gna*mLS**3*hLS*(V-Ena)*LS
      return
      end
      
      subroutine initC(t,spikeg)
      implicit none
      integer nnode,spikeg
      real*8 C,Ena,Ek,El,gleak,gna,gk,dt,t
      common/C1/gna,Ena,gk,Ek,gleak,El
      common/C3/C,dt
      common/node/nnode
      nnode=10
      C=1.0d0
      Ena=50.0d0
      Ek=-77.0d0
      El=-54.4d0
      gleak=0.25d0
      gna=120.0d0
      gk=36.0d0
      dt=2.0d-3
      t=0.0d0
      spikeg=0
      return
      end
      
      subroutine initVar(m,h,n,V,mLS,hLS)
      implicit none
      integer i,nnode
      real*8 m(0:10),h(0:10),n(0:10),V(0:10),mLS,hLS
      common/node/nnode
      do 3 i=0,(nnode-1)
      m(i)=0.095d0
      h(i)=0.414d0
      n(i)=0.328d0
      V(i)=-66.3d0
3     continue
      mLS=m(5)
      hLS=h(5)
      return
      end
