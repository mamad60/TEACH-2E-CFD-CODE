
c**********************************************************************
c                              teach - t
c**********************************************************************
c     A computer program for the calculation of plane or axisymmetric
c     steady two-dimensional recirculating flows
c**********************************************************************
c     subroutine main
      implicit double precision (a-h,o-z)
       dimension hedu(6),hedv(6),hedp(6),hedt(6)
      common
     1/uvel/resoru,nswpu,urfu,dxepu(32),dxpwu(32),sewu(32)
     1/vvel/resorv,nswpv,urfv,dynpv(32),dypsv(32),snsv(32),rcv(32)
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/var/u(22,22),v(22,22),w(22,22),p(22,22),pp(22,22),ut(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/kase t1/uin,tein,edin,flowin,alamda,win,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
     1/wallf/yplusn(22),xplusw(22),taun(22),tauw(22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22),
     1       sp(22,22)
      logical incalu,incalv,incalp,inpro,incald
      open(6,file='TEACHLAM.ANS')
      great=1.e30
      istep=2
      niter=0
      it=22
      jt=22
      nswpu=3
      nswpv=3
      nswpp=5
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
      incalv=.true.
      incalp=.true.
      inpro=.true.
c-----fiuid   properties
      densit=917.000
      viscos=2.900e-1
      prandt=1.0
         uin=.5
c-----pressure   calculation
      ipref=2
      jpref=2
c-----program   control and monitor
      maxit=200
      urfu=0.5
         urfv=0.5
         urfp=1.0
         urfvis=0.7
         sormax=0.01
c-----calculate geometrical quantities and set variables to zero
      call init
c-----initialise variable  fields
      flowin=0.0
      arden=0.0
      do    200   j=2,jstep
      u(2,j)=uin
         arden=0.5*(den(1,j)+den(2,j))*r(j)*sns(j)
 200   flowin=flowin+arden*u(2,j)
      call props
c-----calculate  source-normalisation factors
      snormm=flowin
      snormu=flowin*uin
  300 niter=niter+1
c-----update  main  dependent  variables
      if (incalu) call calcu
      if (incalv) call calcv
      if (incalp) call calcp
c-----update fluid properties
      if (inpro) call props
c-----normalise  residual-source  terms
      resoru=resoru/snormu
      resorv=resorv/snormu
      resort=0.0
      resorm=resorm/snormm
c-----termination  tests
      source=dmax1(resoru,resurv,resorm)
      if (niter.ge.maxit) go to 302
      if (niter.lt.20.or.source.gt.sormax) go to 300
  302 continue
c-------intermediate  output for tecplot -------------------------------
         write(6,5) nI-1,nJ-1
5      format(2x,
     $'title="SHAHB AHMADI"      ',/,
     $2x,'VariabLes = \    ',/,
     $2x,'"x        " \    ',/,
     $2x,'"y        " \    ',/,
     $2x,'"u        " \    ',/,
     $2x,'"v        " \    ',/,
     $2x,'"p        " \    ',/,
     $2x,'Zone t ="  "\    ',/,
     $2x,'i= ',i3,/,
     $2x,'j= ',i3,/,
     $2x,'f=block           ')
        do 111 m=2,nj
        write(6,7000)(x(i),i=2,ni)
7000    format(15(1x,1pe11.3))
111     continue
        write(6,*)
        do 112 k=ni,2,-1
        write(6,7001)(y(k),j=2,nj)
7001    format(15(1x,1pe11.3))
112     continue
        write(6,*)
c       write(6,*)  '--------------- U-VELOCITY OF NODES ---------------'
        write(6,7004)((u(i,j),i=2,ni),j=nj,2,-1)
7004     format(15(1x,1pe11.3))
        write(6,*)
c        write(6,*)  '---------------- V-VELOCITY OF NODES -------------'
         write(6,7005)((v(i,j),i=2,ni),j=nj,2,-1)
7005     format(15(1x,1pe11.3))
        write(6,*)
c        write(6,*)  '---------------- PRESSURE OF NODES --------------'
         write(6,7006)((p(i,j),i=2,ni),j=nj,2,-1)
7006     format(15(1x,1pe11.3))
      stop
      end
c***********************************************************************
      subroutine init
      implicit double precision(a-h,o-z)
      common
     1/uvel/resoru,nswpu,urfu,dxepu(32),dxpwu(32),sewu(32)
     1/vvel/resorv,nswpv,urfv,dynpv(32),dypsv(32),snsv(32),rcv(32)
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/tdis/resore,nswpd,urfe
     1/var/u(22,22),v(22,22),p(22,22),pp(22,22),w(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/kase t1/uin,win,tein,edin,flowin,alamda,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
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
      common
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/var/u(22,22),v(22,22),p(22,22),pp(22,22),w(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
      do 100 i=2,nim1
      do 100 j=2,njm1
      visold=vis(i,j)
      vis(i,j)=viscos
c-----under-relax viscosity
  100 continue
      return
      end
c**********************************************************************
      subroutine calcu
      implicit double precision(a-h,o-z)
      common
     1/uvel/resoru,nswpu,urfu,dxepu(32),dxpwu(32),sewu(32)
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/var/u(22,22),v(22,22),w(22,22),p(22,22),pp(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22),
     1      sp(22,22)
     1/kase t1/uin,win,tein,edin,flowin,alamda,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
chapter  1  1  1  1  1  assembly of coefficients   1  1  1  1  1  1  1
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
      cp=dmax1(0.0D0,smp)
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
      call modu
chapter  3  final coeff. assembly and residual source calculations  3  3
      resoru=0.0
      do 300 i=3,nim1
      do 301 j=2,njm1
      ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
      du(i,j)=du(i,j)/ap(i,j)
      resor=an(i,j)*u(i,j+1)+as(i,j)*u(i,j-1)+ae(i,j)*u(i+1,j)
     1    +aw(i,j)*u(i-1,j)-ap(i,j)*u(i,j)+su(i,j)
      vol=r(j)*sew(i)*sns(j)
      sorvol=great*vol
      if (-sp(i,j).gt.0.5*sorvol) resor=0.0
      resoru=resoru+abs(resor)
c-----under-relaxation
      ap(i,j)=ap(i,j)/urfu
      su(i,j)=su(i,j)+(1.-urfu)*ap(i,j)*u(i,j)
      du(i,j)=du(i,j)*urfu
 301  continue
 300  continue
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
      do  400  n=1,nswpu
 400  call lisolv(3,2,ni,nj,it,jt,u)
      return
      end
c**********************************************************************
      subroutine calcv
      implicit double precision(a-h,o-z)
      common
     1/vvel/resorv,nswpv,urfv,dynpv(32),dypsv(32),snsv(32),rcv(32)
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/var/u(22,22),v(22,22),w(22,22),p(22,22),pp(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22),
     1      sp(22,22)
     1/kase t1/uin,win,tein,edin,flowin,alamda,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
chapter  1  1  1  1  1  assembly of coefficients  1  1  1  1  1  1  1  1
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
     1 +vol*((den(i,j)+den(i,j-1))*(w(i,j)+w(i,j-1))**2)/(8*rv(j))
      sp(i,j)=-cp
      if (indcos.eq.2) sp(i,j)=sp(i,j)-2*vis(i,j)*vol/rv(j)**2
 101  continue
 100  continue
      arean=rcv(3)*sew(2)
      areas=rcv(2)*sew(2)
      dv(2,2)=0.5*(arean+areas)
      call modv
chapter  3  final coeff. assembly and residual source calculations  3  3
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
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
      do 400 n=1,nswpv
 400  call lisolv(2,3,ni,nj,it,jt,v)
      return
      end
c***********************************************************************
        subroutine calcp
      implicit double precision(a-h,o-z)
      common
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/var/u(22,22),v(22,22),w(22,22),p(22,22),pp(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22),
     1      sp(22,22)
      resorm=0.0
chapter  1  1  1  1  1  assembly of coefficients   1  1  1  1  1  1  1
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
      call modp
chapter  3  final coeff. assembly and residual source calculations  3  3
      do 300 i=2,nim1
      do 301 j=2,njm1
  301 ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
  300 continue
chapter  5  5  5  5  correct velocities and pressure  5  5  5  5  5  5
      do 400 n=1,nswpp
  400 call lisolv(2,2,ni,nj,it,jt,pp)
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
c***********************************************************************
      subroutine lisolv(istart,jstart,ni,nj,it,jt,phi)
      implicit double precision(a-h,o-z)
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
      common
     1/uvel/resoru,nswpu,urfu,dxepu(32),dxpwu(32),sewu(32)
     1/vvel/resorv,nswpv,urfv,dynpv(32),dypsv(32),snsv(32),rcv(32)
     1/pcor/resorm,nswpp,urfp,du(22,22),dv(22,22),ipref,jpref
     1/var/u(22,22),v(22,22),w(22,22),p(22,22),pp(22,22)
     1/all/it,jt,ni,nj,nim1,njm1,great
     1/geom/indcos,x(32),y(32),dxep(32),dxpw(32),dynp(32),dyps(32),
     1      sns(32),sew(32),xu(32),yv(32),r(32),rv(32)
     1/flupr/urfvis,viscos,densit,prandt,den(22,22),vis(22,22)
     1/kase t1/uin,win,tein,edin,flowin,alamda,
     2       rsmall,rlarge,al1,al2,jstep,istep,jstp1,jstm1,istp1,istm1
     1/susp/sukd(22,22),spkd(22,22)
     1/coef/ap(22,22),an(22,22),as(22,22),ae(22,22),aw(22,22),su(22,22)
     1     ,sp(22,22)
     1/wallf/yplusn(22),xplusw(22),taun(22),tauw(22)
chapter  1  1  1  1  1  geometrical quantities  1  1  1  1  1  1  1  1
      entry modpro
      return
chpter  2  2  2  2  2  2  2 u  momentum  2  2  2  2  2  2  2  2  2  2  2
      entry modu
c-----top wall
      yp=yv(nj)-y(njm1)
      j=njm1
      do 210 i=3,nim1
      denu=0.5*(den(i,j)+den(i-1,j))
        tmult=viscos/yp
        ueff=sqrt(wavg*wavg+u(i,j)*u(i,j))
        if(u(i,j) .lt.0.) ueff=-ueff
        taun(i)=-tmult*ueff
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
      arden=0.5*(den(nim1,j)+den(nim1-1,j))*r(j)*sns(j)
      ardent=ardent+arden
 204  flow=flow+arden*u(nim1,j)
      uinc=(flowin-flow)/ardent
      do 205 j=2,njm1
 205  u(ni,j)=u(nim1,j)+uinc
      return
chapter  3  3  3  3  3  3  3  3  v  momentum  3  3  3  3  3  3  3  3  3
      entry modv
c-----side wall
      xp=x(istep)-xu(istep)
      i=istep
      do 310 j=jstp1,njm1
      denv=0.5*(den(i,j)+den(i,j-1))
      tmult=viscos/xp
        veff=sqrt(wavg*wavg+v(i,j)*v(i,j))
        if(v(i,j).lt.0.) veff=-veff
        tauw(j)=-tmult*veff
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
chapter  4  4  4  4  4  4  4  pressure correction  4  4  4  4  4  4  4
      entry modp
      return
      end
c*************************************************************************
