
c**************************************************************************
c                           C.F.T. PROGRAME
c**************************************************************************
c   * A computer program for the calculation of plane or axisymmetric
c     steady two-dimensional recirculating flows
c**************************************************************************
c     subroutine main
      implicit double precision (a-h,o-z)
c
chapter  0  0  0  0  0  0  0  0  preliminaries  0  0  0  0  0  0  0  0  0
c
      dimension hedu(6),hedv(6),hedp(6),hedt(6),hedk(6),hedd(6),hedm(6)
     1    ,heda(6),hedb(6)
      common
     1/uvel/resoru,nswpu,urfu,dxepu(32),dxpwu(32),sewu(32)
     1/vvel/resorv,nswpv,urfv,dynpv(32),dypsv(32),snsv(32),rcv(32)
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/ten/resork,nswrk,urfk
     1/tdis/resore,nswpd,urfe
     1/var/u(22,22),v(22,22),p(22,22),pp(22,22),te(22,22),ed(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/kase t1/uin,tein,edin,flowin,alamda,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
     1/turb/gen(22,22),cd,cmu,c1,c2,cappa,elog,pred,prte
     1/wallf/yplusn(22),xplusw(22),taun(22),tauw(22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22),
     1       sp(22,22)
      logical incalu,incalv,incalp,inpro,incalk,incald
      great=1.e30
      istep=2
      niter=0
      it=22
      jt=22
      nswpu=3
      nswpv=3
      nswpp=5
      nswpk=3
      nswpd=3
c     read(70,010)hedu,hedv,hedp,hedt,hedk,hedd,hedm,heda,hedb
c 010 format(6a6)
c
chapter  1  1  1  1  parameters and control indices  1  1  1  1  1  1  1
c
c-----grid
      ni=16
      nj=16
      nim1=ni-1
      njm1=nj-1
      indcos=2
      jstep=5
      istp1=istep+1
      istm1=istep-1
      jstp1=jstep+1
      jstm1=jstep-1
      rlarge=0.0124
      altot=16.0*rlarge
      epsx=1.2
      sumx=.5*epsx**(ni-4)+(epsx**(ni-3)-1.)/(epsx-1.)+0.5
      dx=altot/sumx
      x(1)=-0.5*dx
      x(2)=-x(1)
      do   100  i=3,nim1
      x(i)=x(i-1)+dx
  100 dx=epsx*dx
      x(ni)=x(nim1)-x(ni-2)+x(nim1)
      al1=0.5*(x(istep)+x(istm1))
      al2=altot-al1
      dy=rlarge/float(nj-2)
      y(1)=-0.5*dy
      do   101  j=2,nj
  101 y(j)=y(j-1)+dy
      rsmall=0.5*(y(jstep)+y(jstp1))
c-----dependent   variable   selection
      incalu=.true.
      incalp=.true.
      incalv=.true.
      incalk=.true.
      incald=.true.
      inpro=.true.
c-----fiuid   properties
      densit=1000.0
      viscos=1.002e-3
      prandt=1.0
c-----turbulence  constants
      cmu=0.09
      cd=1.00
      c1=1.44
      c2=1.92
      cappa=0.4187
      elog=9.793
      pred=cappa*cappa/(c2-c1)/(cmu**0.5)
      prte=1.0
c-----boundary  values
      uin=0.5
      ularge=uin*(rsmall/rlarge)**2
      turbin=0.03
      tein=turbin*uin**2
      alamda=0.005
      edin=tein**1.5/(alamda*rlarge)
      reyn=(uin*rsmall*densit)/viscos
c-----pressure   calculation
      ipref=2
      jpref=2
c-----program   control and monitor
      maxit=200
      imon=6
      jmon=6
      urfu=0.5
      urfv=0.5
      urfp=1.0
      urfe=0.7
      urfk=0.7
      urfvis=0.7
      indpri=190
      sormax=0.01
c
chater  2  2  2  2  2  2  initial   operations  2  2  2  2  2  2  2  2
c
c-----calculate geometrical quantities and set variables to zero
      call init
c-----initialise variable  fields
      flowin=0.0
      arden=0.0
      do    200   j=2,jstep
      u(2,j)=uin
      te(1,j)=tein
      ed(1,j)=edin
      arden=0.5*(den(1,j)+den(2,j))*r(j)*sns(j)
  200 flowin=flowin+arden*u(2,j)
      jfin=jstp1
      do   202  i=2,ni
      if (i.ge.istp1) jfin=nj
      factor=(yv(jstp1)*rv(jstp1))/(yv(jfin)*rv(jfin))
      jend=jfin-1
      do   202   j=2,jend
      te(i,j)=tein
      ed(i,j)=edin
  202 u(i,j)=uin*factor
      do   203  i=2,nim1
  203 yplusn(i)=11.0
      do   204   j=jstep,nj
      xplusw(j)=11.0
  204 if (j.eq.jstep) xplusw(j)=0.0
      call props
c-----calculate  source-normalisation factors
      snormm=flowin
      snormu=flowin*uin
c-----initial  output
      write(*,210)
      write(*,220) uin
      re=uin*rsmall*2.0*densit/viscos
      write(*,230) re
      rsorl=rsmall/rlarge
      write(*,240) rsorl
      write(*,250) viscos
      write(*,260) densit
      if (incalu) call print(2,2,ni,nj,it,jt,xu,y,u,hedu)
      if (incalv) call print(2,2,ni,nj,it,jt,x,yv,v,hedv)
      if (incalp) call print(2,2,ni,nj,it,jt,x,y,p,hedp)
      if (incalk) call print(2,2,ni,nj,it,jt,x,y,te,hedk)
      if (incald) call print(2,2,ni,nj,it,jt,x,y,ed,hedd)
c
chapter  3  3  3  3  3  3  iteration loop  3  3  3  3  3  3  3  3  3  3
c
      write(*,310) imon,jmon
  300 niter=niter+1
c-----update  main  dependent  variables
      if (incalu) call calcu
      if (incalv) call calcv
      if (incalp) call calcp
      if (incalk) call calcte
      if (incald) call calced
c-----update fluid properties
      if (inpro) call props
c-----normalise  residual-source  terms
      resoru=resoru/snormu
      resorv=resorv/snormu
      resort=0.0
      resorm=resorm/snormm
c-----intermediate   output
      dummy=0.0
      write(*,311) niter,resoru,resorv,resorm,resort,resork,resore
     1             ,u(imon,jmon),v(imon,jmon),p(imon,jmon),dummy,
     1                te(imon,njm1),ed(imon,njm1)
      if (mod(niter,indpri).ne.0) go to 301
      if (incalu) call print(2,2,ni,nj,it,jt,xu,y,u,hedu)
      if (incalv) call print(2,2,ni,nj,it,jt,x,yv,v,hedv)
      if (incalp) call print(2,2,ni,nj,it,jt,x,y,p,hedp)
      if (incalk) call print(2,2,ni,nj,it,jt,x,y,te,hedk)
      if (incald) call print(2,2,ni,nj,it,jt,x,y,ed,hedd)
      write(*,310) imon,jmon
  301 continue
c-----termination  tests
      source=dmax1(resoru,resurv,resorm)
      if (niter.eq.20.and.source.gt.1.0e4*sormax) go to 302
      if (niter.eq.maxit) go to 302
      if (source.gt.sormax) go to 300
  302 continue
c
chapter  4  4  4  4  4  final operations and output  4  4  4  4  4  4  4
c
      if (incalu) call print(2,2,ni,nj,it,jt,xu,y,u,hedu)
      if (incalv) call print(2,2,ni,nj,it,jt,x,yv,v,hedv)
      if (incalp) call print(2,2,ni,nj,it,jt,x,y,p,hedp)
      if (incalk) call print(2,2,ni,nj,it,jt,x,y,te,hedk)
      if (incald) call print(2,2,ni,nj,it,jt,x,y,ed,hedd)
c-----calculation of non dimensional turbulence  energy and length scale
      if (.not.incalk) go to 500
      do  400  i=2,nim1
      do  400  j=2,njm1
      su(i,j)=te(i,j)*den(i,j)/abs(taun(i))
  400 sp(i,j)=te(i,j)**1.5/ed(i,j)/rlarge
      call print(2,2,ni,nj,it,jt,x,y,su,heda)
      call print(2,2,ni,nj,it,jt,x,y,sp,hedb)
c-----calculation of shear-stress coefficient  along  large duct wall
  500 continue
      write(*,402)
      do  401  i=istep,nim1
      ssc=taun(i)/(1.0*densit*ularge*ularge)
      xud=xu(i)*rlarge/2.
      write(*,403) i,xud,ssc
  401 continue
      stop
c-----format   statements
  210 format(1h0,47x,'kase t1 turbulent flow through a sudden enlargemen
     1t'////)
  220 format(//1h0,15x,'inlet fluid velocity ',t60,1h=,3x,1pe11.3)
  230 format(1h0,15x,'reynolds number ',t60,1h=,3x,1pe11.3)
  240 format(1h0,15x,'diameter ratio ',t60,1h=,3x,1pe11.3)
  250 format(1h0,15x,'laminar viscosity ',t60,1h=,3x,1pe11.3)
  260 format(1h0,15x,'fluid density ',t60,1h=,3x,1pe11.3)
 310  format(1h0,'iter   ','i---------------absolute residual source sum
     1s---------------i   i-------field values at monitoring location','
     2(',i2,',',i2,')','--------i' / 2x,'no,',7x,'umom',6x,'vmom',6x,'ma
     3ss',6x,'ener',6x,'tkin',6x,'disp',10x,'u',9x,'v',9x,'p',9x,'t',9x,
     4'k',9x,'o'/)
  311 format(1h ,i3,4x,1p6e10.3,3x,1p6e10.3)
  402 format(///5x,1hi,7x,5hxu(i),6x,10hs,s,coeff,)
  403 format(/5x,i5,2(1pe11.3))
      end
c************************************************************************
      subroutine init
      implicit double precision(a-h,o-z)
chapter  0  0  0  0  0  0  0  0 preliminaries  0  0  0  0  0  0  0  0  0
c
      common
     1/uvel/resoru,nswpu,urfu,dxepu(32),dxpwu(32),sewu(32)
     1/vvel/resorv,nswpv,urfv,dynpv(32),dypsv(32),snsv(32),rcv(32)
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/ten/resork,nswrk,urfk
     1/tdis/resore,nswpd,urfe
     1/var/u(22,22),v(22,22),p(22,22),pp(22,22),te(22,22),ed(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/kase t1/uin,tein,edin,flowin,alamda,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
     1/turb/gen(22,22),cd,cmu,c1,c2,cappa,elog,pred,prte
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22)
     1     ,sp(22,22)
c
chapter  1  1  1  1  1  geometrical quantities  1  1  1  1  1  1  1  1
c
      do   100  j=1,nj
      r(j)=y(j)
  100 if (indcos.eq.1) r(j)=1.0
      dxpw(1)=0.0
      dxep(ni)=0.0
      do  101  i=1,nim1
      dxep(i)=x(i+1)-x(i)
 101  dxpw(i+1)=dxep(i)
      dyps(1)=0.0
      dynp(nj)=0.0
      do  102  j=1,njm1
      dynp(j)=y(j+1)-y(j)
 102  dyps(j+1)=dynp(j)
      sew(1)=0.0
      sew(ni)=0.0
      do  103  i=2,nim1
 103  sew(i)=0.5*(dxep(i)+dxpw(i))
      sns(1)=0.0
      sns(nj)=0.0
      do  104  j=2,njm1
 104  sns(j)=0.5*(dynp(j)+dyps(j))
      xu(1)=0.0
      do  105 i=2,ni
  105 xu(i)=0.5*(x(i)+x(i-1))
      dxpwu(1)=0.0
      dxpwu(2)=0.0
      dxepu(1)=0.0
      dxepu(ni)=0.0
      do  106  i=2,nim1
      dxepu(i)=xu(i+1)-xu(i)
 106  dxpwu(i+1)=dxepu(i)
      sewu(1)=0.0
      sewu(2)=0.0
      do  107  i=3,nim1
 107  sewu(i)=0.5*(dxepu(i)+dxpwu(i))
      yv(1)=0.0
      rv(1)=0.0
      do  108  j=2,nj
      rv(j)=0.5*(r(j)+r(j-1))
      rcv(j)=0.5*(rv(j)+rv(j-1))
 108  yv(j)=0.5*(y(j)+y(j-1))
      dypsv(1)=0.0
      dypsv(2)=0.0
      dynpv(nj)=0.0
      do  109  j=2,njm1
      dynpv(j)=yv(j+1)-yv(j)
 109  dypsv(j+1)=dynpv(j)
      snsv(1)=0.0
      snsv(2)=0.0
      snsv(nj)=0.0
      do  110  j=3,njm1
  110 snsv(j)=0.5*(dynpv(j)+dypsv(j))
c
chapter  2  2  2  2  2  2  set variable to zero  2  2  2  2  2  2  2
c
      do  200  i=1,ni
      do  200  j=1,nj
      u(i,j)=0.0
      v(i,j)=0.0
      p(i,j)=0.0
      pp(i,j)=0.0
      te(i,j)=1.0e-10
      ed(i,j)=0.001
      den(i,j)=densit
      vis(i,j)=viscos
      du(i,j)=0.0
      dv(i,j)=0.0
      su(i,j)=0.0
      sp(i,j)=0.0
  200 continue
      return
      end
c*********************************************************************
      subroutine props
      implicit double precision(a-h,o-z)
chapter  0  0  0  0  0  0  0  0 preliminaries  0  0  0  0  0  0  0  0  0
c
      common
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/var/u(22,22),v(22,22),p(22,22),pp(22,22),te(22,22),ed(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/turb/gen(22,22),cd,cmu,c1,c2,cappa,elog,pred,prte
c
chapter  1  1  1  1  1   viscosity  1  1  1  1  1  1  1  1  1  1  1  1
c
      do 100 i=2,nim1
      do 100 j=2,njm1
      visold=vis(i,j)
      if (ed(i,j).eq.0.) go to 102
      vis(i,j)=den(i,j)*te(i,j)**2*cmu/ed(i,j)+viscos
      go to 101
 102  vis(i,j)=viscos
c-----under-relax viscosity
  101 vis(i,j)=urfvis*vis(i,j)+(1.-urfvis)*visold
  100 continue
      return
      end
c**********************************************************************
      subroutine calcu
      implicit double precision(a-h,o-z)
chapter  0  0  0  0  0  0  0  0 preliminaries  0  0  0  0  0  0  0  0  0
c
      common
     1/uvel/resoru,nswpu,urfu,dxepu(32),dxpwu(32),sewu(32)
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/var/u(22,22),v(22,22),p(22,22),pp(22,22),te(22,22),ed(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22),
     1      sp(22,22)
     1/kase t1/uin,tein,edin,flowin,alamda,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
c
chapter  1  1  1  1  1  assembly of coefficients   1  1  1  1  1  1  1
c
      do 100 i=3,nim1
      do 101 j=2,njm1
c-----compute areas and volume
      arean=rv(j+1)*sewu(i)
      areas=rv(j)*sewu(i)
      areaew=r(j)*sns(j)
      vol=r(j)*sewu(i)*sns(j)
c-----calculate convection coefficients
      gn=0.5*(den(i,j+1)+den(i,j))*v(i,j+1)
      gnw=0.5*(den(i-1,j)+den(i-1,j+1))*v(i-1,j+1)
      gs=0.5*(den(i,j-1)+den(i,j))*v(i,j)
      gsw=0.5*(den(i-1,j)+den(i-1,j-1))*v(i-1,j)
      ge=0.5*(den(i+1,j)+den(i,j))*u(i+1,j)
      gp=0.5*(den(i,j)+den(i-1,j))*u(i,j)
      gw=0.5*(den(i-1,j)+den(i-2,j))*u(i-1,j)
      cn=0.5*(gn+gnw)*arean
      cs=0.5*(gs+gsw)*areas
      ce=0.5*(ge+gp)*areaew
      cw=0.5*(gp+gw)*areaew
c-----calculate diffusion coefficients
      visn=0.25*(vis(i,j)+vis(i,j+1)+vis(i-1,j)+vis(i-1,j+1))
      viss=0.25*(vis(i,j)+vis(i,j-1)+vis(i-1,j)+vis(i-1,j-1))
      dn=visn*arean/dynp(j)
      ds=viss*areas/dyps(j)
      de=vis(i,j)*areaew/dxepu(i)
      dw=vis(i-1,j)*areaew/dxpwu(i)
c-----calculate coefficients of source terms
      smp=cn-cs+ce-cw
      cp=dmax1(0.0d0,smp)
      cpo=cp
c-----assemble main coefficients
      an(i,j)=dmax1(abs(0.5*cn),dn)-0.5*cn
      as(i,j)=dmax1(abs(0.5*cs),ds)+0.5*cs
      ae(i,j)=dmax1(abs(0.5*ce),de)-0.5*ce
      aw(i,j)=dmax1(abs(0.5*cw),dw)+0.5*cw
      du(i,j)=areaew
      su(i,j)=cpo*u(i,j)+du(i,j)*(p(i-1,j)-p(i,j))
      sp(i,j)=-cp
  101 continue
  100 continue
c
chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
c
      call modu
c
chapter  3  final coeff. assembly and residual source calculations  3  3
c
      resoru=0.0
      do 300 i=3,nim1
      do 301 j=2,njm1
      ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
      du(i,j)=du(i,j)/ap(i,j)
      resor=an(i,j)*u(i,j+1)+as(i,j)*u(i,j-1)+ae(i,j)*u(i+1,j)
     1    +aw(i,j)*u(i-1,j)-ap(i,j)*u(i,j)+su(i,j)
      vol=r(j)*sew(i)*sns(j)
      sorvol=great*vol
      if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
      resoru=resoru+abs(resor)
c-----under-relaxation
      ap(i,j)=ap(i,j)/urfu
      su(i,j)=su(i,j)+(1.-urfu)*ap(i,j)*u(i,j)
      du(i,j)=du(i,j)*urfu
 301  continue
 300  continue
c
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
c
      do  400  n=1,nswpu
 400  call lisolv(3,2,ni,nj,it,jt,u)
      return
      end
c**********************************************************************
      subroutine calcv
      implicit double precision(a-h,o-z)
chapter  0  0  0  0  0  0  0  0 preliminaries  0  0  0  0  0  0  0  0  0
c
      common
     1/vvel/resorv,nswpv,urfv,dynpv(32),dypsv(32),snsv(32),rcv(32)
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/var/u(22,22),v(22,22),p(22,22),pp(22,22),te(22,22),ed(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22),
     1      sp(22,22)
     1/kase t1/uin,tein,edin,flowin,alamda,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
c
chapter  1  1  1  1  1  assembly of coefficients  1  1  1  1  1  1  1  1
c
      do 100 i=2,nim1
      do 101 j=3,njm1
c-----compute areas and volume
      arean=rcv(j+1)*sew(i)
      areas=rcv(j)*sew(i)
      areaew=rv(j)*snsv(j)
      vol=rv(j)*sew(i)*snsv(j)
c-----calculate convection coefficients
      gn=0.5*(den(i,j+1)+den(i,j))*v(i,j+1)
      gp=0.5*(den(i,j)+den(i,j-1))*v(i,j)
      gs=0.5*(den(i,j-1)+den(i,j-2))*v(i,j-1)
      ge=0.5*(den(i+1,j)+den(i,j))*u(i+1,j)
      gse=0.5*(den(i,j-1)+den(i+1,j-1))*u(i+1,j-1)
      gw=0.5*(den(i,j)+den(i-1,j))*u(i,j)
      gsw=0.5*(den(i,j-1)+den(i-1,j-1))*u(i,j-1)
      cn=0.5*(gn+gp)*arean
      cs=0.5*(gs+gp)*areas
      ce=0.5*(ge+gse)*areaew
      cw=0.5*(gw+gsw)*areaew
c-----calculate diffusion coefficients
      vise=0.25*(vis(i,j)+vis(i+1,j)+vis(i,j-1)+vis(i+1,j-1))
      visw=0.25*(vis(i,j)+vis(i-1,j)+vis(i,j-1)+vis(i-1,j-1))
      dn=vis(i,j)*arean/dynpv(j)
      ds=vis(i,j-1)*areas/dypsv(j)
      de=vise*areaew/dxep(i)
      dw=visw*areaew/dxpw(i)
c-----calculate coefficients of source terms
      smp=cn-cs+ce-cw
      cp=dmax1(0.0d0,smp)
      cpo=cp
c-----assemble main coefficients
      an(i,j)=dmax1(abs(0.5*cn),dn)-0.5*cn
      as(i,j)=dmax1(abs(0.5*cs),ds)+0.5*cs
      ae(i,j)=dmax1(abs(0.5*ce),de)-0.5*ce
      aw(i,j)=dmax1(abs(0.5*cw),dw)+0.5*cw
      dv(i,j)=0.5*(arean+areas)
      su(i,j)=cpo*v(i,j)+dv(i,j)*(p(i,j-1)-p(i,j))
      sp(i,j)=-cp
      if (indcos.eq.2) sp(i,j)=sp(i,j)-vis(i,j)*vol/rv(j)**2
 101  continue
 100  continue
      arean=rcv(3)*sew(2)
      areas=rcv(2)*sew(2)
      dv(2,2)=0.5*(arean+areas)
c
chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
c
      call modv
c
chapter  3  final coeff. assembly and residual source calculations  3  3
c
      resorv=0.0
      do 300 i=2,nim1
      do 301 j=3,njm1
      ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
      dv(i,j)=dv(i,j)/ap(i,j)
      resor=an(i,j)*v(i,j+1)+as(i,j)*v(i,j-1)+ae(i,j)*v(i+1,j)
     1    +aw(i,j)*v(i-1,j)-ap(i,j)*v(i,j)+su(i,j)
      vol=r(j)*sew(i)*sns(j)
      sorvol=great*vol
      if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
      resorv=resorv+abs(resor)
c-----under-relaxation
      ap(i,j)=ap(i,j)/urfv
      su(i,j)=su(i,j)+(1.-urfv)*ap(i,j)*v(i,j)
      dv(i,j)=dv(i,j)*urfv
  301 continue
  300 continue
c
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
c
      do 400 n=1,nswpv
 400  call lisolv(2,3,ni,nj,it,jt,v)
      return
      end
c*************************************************************************
      subroutine calcp
      implicit double precision(a-h,o-z)
chapter  0  0  0  0  0  0  0  0 preliminaries  0  0  0  0  0  0  0  0  0
c
      common
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/var/u(22,22),v(22,22),p(22,22),pp(22,22),te(22,22),ed(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22),
     1      sp(22,22)
      resorm=0.0
c
chapter  1  1  1  1  1  assembly of coefficients   1  1  1  1  1  1  1
c
      do 100 i=2,nim1
      do 101 j=2,njm1
c-----compute areas and volume
      arean=rv(j+1)*sew(i)
      areas=rv(j)*sew(i)
      areaew=r(j)*sns(j)
      vol=r(j)*sns(j)*sew(i)
c-----calculate coefficients
      denn=0.5*(den(i,j)+den(i,j+1))
      dens=0.5*(den(i,j)+den(i,j-1))
      dene=0.5*(den(i,j)+den(i+1,j))
      denw=0.5*(den(i,j)+den(i-1,j))
      an(i,j)=denn*arean*dv(i,j+1)
      as(i,j)=dens*areas*dv(i,j)
      ae(i,j)=dene*areaew*du(i+1,j)
      aw(i,j)=denw*areaew*du(i,j)
c-----calculate source term
      cn=denn*v(i,j+1)*arean
      cs=dens*v(i,j)*areas
      ce=dene*u(i+1,j)*areaew
      cw=denw*u(i,j)*areaew
      smp=cn-cs+ce-cw
      sp(i,j)=0.0
      su(i,j)=-smp
c-----compute sum of absolute mass sources
      resorm=resorm+abs(smp)
 101  continue
 100  continue
c
chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
c
      call modp
c
chapter  3  final coeff. assembly and residual source calculations  3  3
c
      do 300 i=2,nim1
      do 301 j=2,njm1
  301 ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
  300 continue
c
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
c
chapter  5  5  5  5  correct velocities and pressure  5  5  5  5  5  5  5
      do 400 n=1,nswpp
  400 call lisolv(2,2,ni,nj,it,jt,pp)
c
c-----velocities
      do 500 i=2,nim1
      do 501 j=2,njm1
      if (i.ne.2) u(i,j)=u(i,j)+du(i,j)*(pp(i-1,j)-pp(i,j))
      if (j.ne.2) v(i,j)=v(i,j)+dv(i,j)*(pp(i,j-1)-pp(i,j))
  501 continue
  500 continue
c-----pressures (with provision for under-relaxation)
      ppref=pp(ipref,jpref)
      do 502 i=2,nim1
      do 503 j=2,njm1
      p(i,j)=p(i,j)+urfp*(pp(i,j)-ppref)
      pp(i,j)=0.0
  503 continue
  502 continue
      return
      end
c*************************************************************************
      subroutine calcte
      implicit double precision(a-h,o-z)
chapter  0  0  0  0  0  0  0  0 preliminaries  0  0  0  0  0  0  0  0  0
c
       common
     1/ten/resork,nswpk,urfk
     1/var/u(22,22),v(22,22),p(22,22),pp(22,22),te(22,22),ed(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/kase t1/uin,tein,edin,flowin,alamda,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
     1/turb/gen(22,22),cd,cmu,c1,c2,cappa,elog,pred,prte
     1/wallf/yplusn(22),xplusw(22),taun(22),tauw(22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22),
     1      sp(22,22)
     1/susp/sukd(22,22),spkd(22,22)
c
chapter  1  1  1  1  1  assembly of coefficients   1  1  1  1  1  1  1
c
      prte=1.0
      do 100 i=2,nim1
      do 101 j=2,njm1
c-----compute areas and volume
      arean=rv(j+1)*sew(i)
      areas=rv(j)*sew(i)
      areaew=r(j)*sns(j)
      vol=r(j)*sns(j)*sew(i)
c-----calculate convection coefficients
      gn=0.5*(den(i,j)+den(i,j+1))*v(i,j+1)
      gs=0.5*(den(i,j)+den(i,j-1))*v(i,j)
      ge=0.5*(den(i,j)+den(i+1,j))*u(i+1,j)
      gw=0.5*(den(i,j)+den(i-1,j))*u(i,j)
      cn=gn*arean
      cs=gs*areas
      ce=ge*areaew
      cw=gw*areaew
c-----calculate diffusion coefficients
      gamn=0.5*(vis(i,j)+vis(i,j+1))/prte
      gams=0.5*(vis(i,j)+vis(i,j-1))/prte
      game=0.5*(vis(i,j)+vis(i+1,j))/prte
      gamw=0.5*(vis(i,j)+vis(i-1,j))/prte
      dn=gamn*arean/dynp(j)
      ds=gams*areas/dyps(j)
      de=game*areaew/dxep(i)
      dw=gamw*areaew/dxpw(i)
c-----source terms
      smp=cn-cs+ce-cw
      cp=dmax1(0.0d0,smp)
      cpo=cp
      dudx=(u(i+1,j)-u(i,j))/sew(i)
      dvdy=(v(i,j+1)-v(i,j))/sns(j)
      dudy=((u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))/4.-(u(i,j)+u(i+1,j)+
     1u(i,j-1)+u(i+1,j-1))/4.)/sns(j)
      dvdx=((v(i,j)+v(i,j+1)+v(i+1,j)+v(i+1,j+1))/4.-(v(i,j)+v(i,j+1)+v(
     1i-1,j)+v(i-1,j+1))/4.)/sew(i)
      gen(i,j)=(2.*(dudx**2+dvdy**2)+(dudy+dvdx)**2)*vis(i,j)
      if (j.gt.2) vdr=v(i,j)/rv(j)
      if (j.eq.2) vdr=0.0
      if (indcos.eq.2) gen(i,j)=gen(i,j)+vis(i,j)*0.5*(vdr+v(i,j+1)/
     1                 rv(j+1))**2
c-----assemble main coefficients
      an(i,j)=dmax1(abs(0.5*cn),dn)-0.5*cn
      as(i,j)=dmax1(abs(0.5*cs),ds)+0.5*cs
      ae(i,j)=dmax1(abs(0.5*ce),de)-0.5*ce
      aw(i,j)=dmax1(abs(0.5*cw),dw)+0.5*cw
      su(i,j)=cpo*te(i,j)
      sukd(i,j)=su(i,j)
      su(i,j)=su(i,j)+gen(i,j)*vol
      sp(i,j)=-cp
      spkd(i,j)=sp(i,j)
      sp(i,j)=sp(i,j)-cd*cmu*den(i,j)**2*te(i,j)*vol/vis(i,j)
  101 continue
  100 continue
c
chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
c
      call modte
c
chapter  3  final coeff. assembly and residual source calculations  3  3
c
      resork=0.0
      do 300 i=2,nim1
      do 301 j=2,njm1
      ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
      resor=an(i,j)*te(i,j+1)+as(i,j)*te(i,j-1)+ae(i,j)*te(i+1,j)
     1    +aw(i,j)*te(i-1,j)-ap(i,j)*te(i,j)+su(i,j)
      vol=r(j)*sew(i)*sns(j)
      sorvol=great*vol
      if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
      resork=resork+abs(resor)
c-----under-relaxation
      ap(i,j)=ap(i,j)/urfk
      su(i,j)=su(i,j)+(1.-urfk)*ap(i,j)*te(i,j)
 301  continue
 300  continue
c
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
      do 400 n=1,nswpk
  400 call lisolv(2,2,ni,nj,it,jt,te)
      return
      end
c***********************************************************************
      subroutine calced
      implicit double precision(a-h,o-z)
chapter  0  0  0  0  0  0  0  0 preliminaries  0  0  0  0  0  0  0  0  0
c
      common
     1/tdis/resore,nswpd,urfe
     1/var/u(22,22),v(22,22),p(22,22),pp(22,22),te(22,22),ed(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/susp/sukd(22,22),spkd(22,22)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/kase t1/uin,tein,edin,flowin,alamda,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
     1/turb/gen(22,22),cd,cmu,c1,c2,cappa,elog,pred,prte
     1/wallf/yplusn(22),xplusw(22),taun(22),tauw(22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22),
     1      sp(22,22)
c
chapter  1  1  1  1  1  assembly of coefficients   1  1  1  1  1  1  1
c
      do 100 i=2,nim1
      do 101 j=2,njm1
c-----compute areas and volume
      arean=rv(j+1)*sew(i)
      areas=rv(j)*sew(i)
      areaew=r(j)*sns(j)
      vol=r(j)*sns(j)*sew(i)
c-----calculate convection coefficients
      gn=0.5*(den(i,j)+den(i,j+1))*v(i,j+1)
      gs=0.5*(den(i,j)+den(i,j-1))*v(i,j)
      ge=0.5*(den(i,j)+den(i+1,j))*u(i+1,j)
      gw=0.5*(den(i,j)+den(i-1,j))*u(i,j)
      cn=gn*arean
      cs=gs*areas
      ce=ge*areaew
      cw=gw*areaew
c-----calculate diffusion coefficients
      gamn=0.5*(vis(i,j)+vis(i,j+1))/pred
      gams=0.5*(vis(i,j)+vis(i,j-1))/pred
      game=0.5*(vis(i,j)+vis(i+1,j))/pred
      gamw=0.5*(vis(i,j)+vis(i-1,j))/pred
      dn=gamn*arean/dynp(j)
      ds=gams*areas/dyps(j)
      de=game*areaew/dxep(i)
      dw=gamw*areaew/dxpw(i)
c-----source terms
      smp=cn-cs+ce-cw
      cp=dmax1(0.0d0,smp)
      cpo=cp
c-----assemble main coefficients
      an(i,j)=dmax1(abs(0.5*cn),dn)-0.5*cn
      as(i,j)=dmax1(abs(0.5*cs),ds)+0.5*cs
      ae(i,j)=dmax1(abs(0.5*ce),de)-0.5*ce
      aw(i,j)=dmax1(abs(0.5*cw),dw)+0.5*cw
      su(i,j)=cpo*ed(i,j)
      sukd(i,j)=su(i,j)
      su(i,j)=su(i,j)+c1*cmu*gen(i,j)*vol*den(i,j)*te(i,j)/vis(i,j)
      sp(i,j)=-cp
      spkd(i,j)=sp(i,j)
      sp(i,j)=sp(i,j)-c2*den(i,j)*ed(i,j)*vol/te(i,j)
 101  continue
 100  continue
c
chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
c
      call moded
c
chapter  3  final coeff. assembly and residual source calculations  3  3
c
      resore=0.0
      do 300 i=2,nim1
      do 301 j=2,njm1
      ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
      resor=an(i,j)*ed(i,j+1)+as(i,j)*ed(i,j-1)+ae(i,j)*ed(i+1,j)
     1    +aw(i,j)*ed(i-1,j)-ap(i,j)*ed(i,j)+su(i,j)
      vol=r(j)*sew(i)*sns(j)
      sorvol=great*vol
      if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
      resore=resore+abs(resor)
c-----under-relaxation
      ap(i,j)=ap(i,j)/urfe
      su(i,j)=su(i,j)+(1.-urfe)*ap(i,j)*ed(i,j)
 301  continue
 300  continue
c
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
c
      do 400 n=1,nswpd
 400  call lisolv(2,2,ni,nj,it,jt,ed)
      return
      end
c***********************************************************************
      subroutine lisolv(istart,jstart,ni,nj,it,jt,phi)
      implicit double precision(a-h,o-z)
chapter  0  0  0  0  0  0  0  0 preliminaries  0  0  0  0  0  0  0  0  0
c
      dimension phi(it,jt),a(32),b(32),c(32),d(32)
      common
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22),
     1      sp(22,22)
      nim1=ni-1
      njm1=nj-1
      jstm1=jstart-1
      a(jstm1)=0.0
c-----commence w-e sweep
      do 100 i=istart,nim1
      c(jstm1)=phi(i,jstm1)
c-----commence s-n traverse
      do 101 j=jstart,njm1
c-----assemble toma coefficients
      a(j)=an(i,j)
      b(j)=as(i,j)
      c(j)=ae(i,j)*phi(i+1,j)+aw(i,j)*phi(i-1,j)+su(i,j)
      d(j)=ap(i,j)
c-----calculate coefficients of recorrence formula
      term=1./(d(j)-b(j)*a(j-1))
      a(j)=a(j)*term
 101  c(j)=(c(j)+b(j)*c(j-1))*term
c-----obtain new phi's
      do 102 jj=jstart,njm1
      j=nj+jstm1-jj
 102  phi(i,j)=a(j)*phi(i,j+1)+c(j)
 100  continue
      return
      end
c*************************************************************************
      subroutine promod
      implicit double precision(a-h,o-z)
chapter  0  0  0  0  0  0  0  0 preliminaries  0  0  0  0  0  0  0  0  0
c
      common
     1/uvel/resoru,nswpu,urfu,dxepu(32),dxpwu(32),sewu(32)
     1/vvel/resorv,nswpv,urfv,dynpv(32),dypsv(32),snsv(32),rcv(32)
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/var/u(22,22),v(22,22),p(22,22),pp(22,22),te(22,22),ed(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/kase t1/uin,tein,edin,flowin,alamda,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
     1/susp/sukd(22,22),spkd(22,22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22)
     1     ,sp(22,22)
     1/turb/gen(22,22),cd,cmu,c1,c2,cappa,elog,pred,prte
     1/wallf/yplusn(22),xplusw(22),taun(22),tauw(22)
c
chapter  1  1  1  1  1  geometrical quantities  1  1  1  1  1  1  1  1
c
      entry modpro
c-----no modifications for this problem
      return
c
chpter  2  2  2  2  2  2  2 u  momentum  2  2  2  2  2  2  2  2  2  2  2
c
      entry modu
c-----top wall
      cdterm=cmu**0.25
      yp=yv(nj)-y(njm1)
      j=njm1
      do 210 i=3,nim1
      sqrtk=sqrt(0.5*(te(i,j)+te(i-1,j)))
      denu=0.5*(den(i,j)+den(i-1,j))
      yplusa=0.5*(yplusn(i)+yplusn(i-1))
      if (yplusa.le.11.63)  go to 211
      tmult=denu*cdterm*sqrtk*cappa/log(elog*yplusa)
      go to 212
 211  tmult=viscos/yp
 212  taun(i)=-tmult*u(i,j)
      sp(i,j)=sp(i,j)-tmult*sewu(i)*rv(nj)
 210  an(i,j)=0.0
      taun(2)=taun(3)
      taun(ni)=taun(nim1)
c-----side wall
      do 213 j=jstp1,njm1
 213  aw(3,j)=0.0
c-----symmetry  axis
      do 203 i=1,ni
 203  as(i,2)=0.0
c-----outlet
      ardent=0.0
      flow=0.0
      do 204 j=2,njm1
      arden=0.5*(den(njm1,j)+den(nim1-1,j))*r(j)*sns(j)
      ardent=ardent+arden
 204  flow=flow+arden*u(nim1,j)
      uinc=(flowin-flow)/ardent
      do 205 j=2,njm1
 205  u(ni,j)=u(nim1,j)+uinc
      return
c
chapter  3  3  3  3  3  3  3  3  v  momentum  3  3  3  3  3  3  3  3  3
c
      entry modv
c-----side wall
      cdterm=cmu**0.25
      xp=x(istep)-xu(istep)
      i=istep
      do 310 j=jstp1,njm1
      sqrtk=sqrt(0.5*(te(i,j)+te(i,j-1)))
      denv=0.5*(den(i,j)+den(i,j-1))
      xplusa=0.5*(xplusw(j)+xplusw(j-1))
      if (xplusa.le.11.63)  go to 311
      tmult=denv*cdterm*sqrtk*cappa/log(elog*xplusa)
      go to 312
 311  tmult=viscos/xp
 312  tauw(j)=-tmult*v(i,j)
      sp(i,j)=sp(i,j)-tmult*snsv(j)*rv(j)
 310  aw(i,j)=0.0
      tauw(jstep)=tauw(jstp1)
      tauw(nj)=tauw(njm1)
c-----top wall
      do 313 i=istep,nim1
 313  an(i,njm1)=0.0
c-----symmetry  axis
      do 302 i=2,nim1
  302 as(i,3)=0.0
      return
c
chapter  4  4  4  4  4  4  4  pressure correction  4  4  4  4  4  4  4
c
      entry modp
      return
c
chapter  5  5  5  5  5  5  5  5  thermal energy  5  5  5  5  5  5  5  5
c
      entry modt
c-----no modifications for this problem
      return
c
chapter  6  6  6  6  6  6  6  turbulent kinetic energy   6  6  6  6  6
c
      entry modte
c-----top wall
      cdterm=cmu**0.25
      yp=yv(nj)-y(njm1)
      j=njm1
      do  610  i=2,nim1
      denu=den(i,j)
      sqrtk=sqrt(te(i,j))
      vol=r(j)*sns(j)*sew(i)
      gencou=0.5*(abs(taun(i+1)*u(i+1,j))+abs(taun(i)*u(i,j)))/yp
      yplusn(i)=denu*sqrtk*cdterm*yp/viscos
      dudy=((u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))/4.-(u(i,j)+u(i+1,j)+
     1u(i,j-1)+u(i+1,j-1))/4.)/sns(j)
      genres=gen(i,j)-vis(i,j)*dudy**2
      gen(i,j)=genres+gencou
      if (yplusn(i).le.11.63)  go to 611
      diterm=den(i,j)*(cmu**0.75)*sqrtk*log(elog*yplusn(i))/(cappa*yp)
      go to 612
 611  continue
      diterm=den(i,j)*(cmu**0.75)*sqrtk*yplusn(i)/yp
 612  continue
      su(i,j)=gen(i,j)*vol+sukd(i,j)
      sp(i,j)=-diterm*vol+spkd(i,j)
 610  an(i,j)=0.0
c-----side wall
      xp=x(istep)-xu(istep)
      i=istep
      do 620 j=jstp1,njm1
      denv=den(i,j)
      sqrtk=sqrt(te(i,j))
      vol=r(j)*sns(j)*sew(i)
      xplusw(j)=denv*sqrtk*cdterm*xp/viscos
      gencou=0.5*(abs(tauw(j+1)*v(i,j+1))+abs(tauw(j)*v(i,j)))/xp
      dvdx=((v(i,j)+v(i,j+1)+v(i+1,j)+v(i+1,j+1))/4.-(v(i,j)+v(i,j+1)+v(
     1i-1,j)+v(i-1,j+1))/4.)/sew(i)
      genres=gen(i,j)-vis(i,j)*dvdx**2
      gen(i,j)=genres+gencou
      if (xplusw(j).le.11.63) go to 621
      diterm=den(i,j)*(cmu**0.75)*sqrtk*log(elog*xplusw(j))/(cappa*xp)
      go to 622
 621  continue
      diterm=den(i,j)*(cmu**0.75)*sqrtk*xplusw(j)/xp
 622  continue
      su(i,j)=sukd(i,j)+gen(i,j)*vol+su(i,j)
      sp(i,j)=spkd(i,j)-diterm*vol+sp(i,j)
 620  aw(i,j)=0.0
c-----symmetry  axis
      j=2
      do  630  i=2,nim1
      dudy=((u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))/4.-(u(i,j)+u(i+1,j)+
     1u(i,j-1)+u(i+1,j-1))/4.)/sns(j)
      vol=r(j)*sns(j)*sew(i)
      gen(i,j)=gen(i,j)-vis(i,j)*dudy**2
      su(i,j)=sukd(i,j)+gen(i,j)*vol
 630  as(i,2)=0.0
      return
c
chapter  7  7  7  7  7  7  7  7  dissipation  7  7  7  7  7  7  7  7  7
c
      entry moded
c-----top wall
      yp=yv(nj)-y(njm1)
      j=njm1
      term=(cmu**0.75)/(cappa*yp)
      do 710 i=2,nim1
      su(i,j)=great*term*te(i,j)**1.5
 710  sp(i,j)=-great
c-----side wall
      xp=x(istep)-xu(istep)
      i=istep
      term=(cmu**.75)/(cappa*xp)
      njm2=nj-2
      do 720 j=jstp1,njm2
      su(i,j)=great*term*te(i,j)**1.5
 720  sp(i,j)=-great
c-----symmetry axis
      do 730 i=2,nim1
 730  as(i,2)=0.0
      return
      end
c%options ansi77
c***********************************************************************
      subroutine print(istart,jstart,ni,nj,it,jt,x,y,phi,head)
      implicit double precision(a-h,o-z)
      character*6 hi,hy
      dimension phi(it,jt),x(it),y(jt),head(6),store(50)
      character*28 f
      character*4 f4(14)
      data f/'(1h  ,a6,i3, 14i 10, 7x, a6)'/
      data f4/' 1i ',' 2i ',' 3i ',' 4i ',' 5i ',' 6i ',
     1        ' 7i ',' 8i ',' 9i ','10i ','11i ','12i ',
     2        '13i ','14i '/
      data hi,hy/'  i = ',' y =  '/
      iskip=1
      jskip=1
      write(*,110) head
      ista=istart-15
 100  continue
      ista=ista+15
      iend=ista+14
      iend=min0(ni,iend)
      f(13:16)=f4(iend-ista)
      write(*,f) hi,  (i,i=ista,iend,iskip), hy
      write(*,112)
      do 101 jj=jstart,nj,jskip
      j=jstart+nj-jj
      do 120 i=ista,iend
      a=phi(i,j)
      if (abs(a).lt.1.e-20) a=0.0
 120  store(i)=a
 101  write(*,113) j,(store(i),i=ista,iend,iskip),y(j)
      write(*,114)  (x(i),i=ista,iend,iskip)
c-----------------------------------------------------------
      if (iend.lt.ni) go to 100
      return
 110  format(1h0,30(2h*-),7x,6a6,7x,30(2h-*))
 111  format(1h0,6h  i = ,i3,14i10,7x,' y = ')
 112  format(3h  j)
 113  format(1h ,i3,1p15e10.2,0pf7.3)
 114  format(6h0x=   ,f7.3,14f10.3)
      end
