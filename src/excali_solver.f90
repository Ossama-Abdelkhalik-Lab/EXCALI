! EXCALI (EXact Cauchy-Analytic Lambert Integration): Lambert's problem for
! every conic and revolution count from the zero-height limit of Exact v2.2's
! thin-rectangle ratio, root = c + Q[((x-c)f)']/Q[f'], f = 1/(T-t), by real
! Gauss-Legendre quadrature. Zero revolutions: a table or closed-form starter
! and one fixed rectangle, no iteration and no fallback. Multiple revolutions:
! a table starter and one certified rectangle, or fixed ratios over regions
! chosen from the geometry.
! Geometry and time functions adapted from ConfluentContour's contour_simple.
module excali_solver
   use iso_fortran_env, only: real64, int64
   use ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_fma
   use iso_c_binding, only: c_double
   implicit none
   private
   public :: excali_one, excali_solve, excali_batch
   integer,parameter :: dp=real64
   real(dp),parameter :: pi=acos(-1.0_dp),sqrt2=sqrt(2.0_dp),rsqrt2=1/sqrt2
   ! The bit-pattern helpers (pow2m, kexp, rlog, rcbrt) assume IEEE-754
   ! binary64; on any other real format they use the standard intrinsics.
   logical,parameter  :: ieee64=radix(1._dp)==2.and.digits(1._dp)==53.and. &
      maxexponent(1._dp)==1024.and.minexponent(1._dp)==-1021.and.bit_size(0_int64)==64
   integer,parameter  :: nb=32          ! excali_batch chunk size
   ! Thin-rectangle quadrature in its zero-height limit. On [c-L,c+L],
   ! L=2**(-kr), with f=1/(T-t) real, Gauss-Legendre nodes x and weights w:
   !   root = c + Q[((x-c)f)'] / Q[f'],   Q[g'] = g(c+L)-g(c-L) - L*sum w g'(x),
   ! the quadrature error of an exact derivative. It is exact for a simple
   ! pole anywhere inside and converges like (L/2R)**(2n) for singularities
   ! at distance R. Node trig comes from angle addition with compile-time
   ! half-angle constants (positive nodes; the last entry is the endpoint).
   integer,parameter  :: kmax=62
   integer :: j_,k_
   ! Stage 1: 8-point Gauss-Legendre (positive half, then the endpoint 1).
   integer,parameter  :: n1=8,h1=n1/2
   real(dp),parameter :: g8x(h1+1)=[.1834346424956498_dp,.5255324099163290_dp,.7966664774136267_dp,.9602898564975363_dp,1._dp]
   real(dp),parameter :: g8w(h1)=[.3626837833783620_dp,.3137066458778873_dp,.2223810344533745_dp,.1012285362903763_dp]
   real(dp),parameter :: g3x(3)=[0._dp,.7745966692414834_dp,1._dp],g3w(2)=[8._dp/9,5._dp/9]
   ! Full node lists (nodes, then the endpoints +1,-1) with their weights,
   ! and cos/sin (cosh/sinh) of the half offsets L*u/2 for each L=2**(-k).
   ! They are protected module data, initialized at compile time: indexing a
   ! parameter array at run time makes some gfortran versions copy it to the
   ! stack.
   real(dp),parameter :: u8(n1+2)=[g8x(1:h1),-g8x(1:h1),1._dp,-1._dp],w8(n1+2)=[g8w,g8w,0._dp,0._dp]
   real(dp),parameter :: u3(5)=[0._dp,g3x(2),-g3x(2),1._dp,-1._dp],w3(5)=[g3w(1),g3w(2),g3w(2),0._dp,0._dp]
   real(dp),protected :: c8f(n1+2,0:kmax)=reshape([((cos(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   real(dp),protected :: s8f(n1+2,0:kmax)=reshape([((sin(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   real(dp),protected :: ch8f(n1+2,0:kmax)=reshape([((cosh(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   real(dp),protected :: sh8f(n1+2,0:kmax)=reshape([((sinh(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   real(dp),protected :: c3f(5,0:kmax)=reshape([((cos(.5_dp**(k_+1)*u3(j_)),j_=1,5),k_=0,kmax)],[5,kmax+1])
   real(dp),protected :: s3f(5,0:kmax)=reshape([((sin(.5_dp**(k_+1)*u3(j_)),j_=1,5),k_=0,kmax)],[5,kmax+1])
   real(dp),protected :: ch3f(5,0:kmax)=reshape([((cosh(.5_dp**(k_+1)*u3(j_)),j_=1,5),k_=0,kmax)],[5,kmax+1])
   real(dp),protected :: sh3f(5,0:kmax)=reshape([((sinh(.5_dp**(k_+1)*u3(j_)),j_=1,5),k_=0,kmax)],[5,kmax+1])
   ! Reciprocal factorials 1/j!, so the Taylor series multiply by constants.
   real(dp),parameter :: rf(16)=[(1/gamma(j_+1._dp),j_=1,16)]
   ! 32-point Gauss-Legendre rule (positive nodes and weights on [-1,1]) for the
   ! multi-revolution minimum, and on v in [0,pi]: gv, gw*pi/2, sin/cos(gv/2).
   integer,parameter  :: ng=32,hg=ng/2
   real(dp),parameter :: ggx(hg)=[0.04830766568773832_dp,0.1444719615827965_dp,0.23928736225213706_dp, &
      0.33186860228212767_dp,0.42135127613063533_dp,0.5068999089322294_dp,0.5877157572407623_dp, &
      0.6630442669302152_dp,0.7321821187402897_dp,0.7944837959679424_dp,0.84936761373257_dp, &
      0.8963211557660521_dp,0.9349060759377397_dp,0.9647622555875064_dp,0.9856115115452683_dp,0.9972638618494816_dp]
   real(dp),parameter :: ggw(hg)=[0.09654008851472758_dp,0.0956387200792746_dp,0.09384439908080439_dp, &
      0.09117387869576371_dp,0.08765209300440369_dp,0.08331192422694662_dp,0.07819389578707012_dp, &
      0.0723457941088482_dp,0.0658222227763616_dp,0.05868409347853546_dp,0.05099805926237596_dp, &
      0.042835898022226926_dp,0.03427386291302163_dp,0.02539206530926273_dp,0.016274394730906486_dp,0.007018610009470362_dp]
   ! 24-point rule for the region formulas (ratio_root).
   integer,parameter  :: nz=24,hz=nz/2
   real(dp),parameter :: gzx(hz)=[0.06405689286260563_dp,0.1911188674736163_dp,0.3150426796961634_dp, &
      0.4337935076260451_dp,0.5454214713888396_dp,0.6480936519369755_dp,0.7401241915785544_dp,0.820001985973903_dp, &
      0.8864155270044011_dp,0.9382745520027328_dp,0.9747285559713095_dp,0.9951872199970213_dp]
   real(dp),parameter :: gzw(hz)=[0.1279381953467523_dp,0.1258374563468285_dp,0.12167047292780352_dp, &
      0.11550566805372571_dp,0.10744427011596573_dp,0.09761865210411409_dp,0.08619016153195343_dp, &
      0.07334648141108027_dp,0.059298584915436534_dp,0.04427743881741942_dp,0.028531388628933684_dp,0.012341229799986784_dp]
   real(dp),parameter :: gvp(ng)=[(pi/2*(1+ggx(j_)),j_=1,hg),(pi/2*(1-ggx(j_)),j_=1,hg)]
   real(dp),protected :: gv(ng)=gvp,gw(ng)=[ggw*(pi/2),ggw*(pi/2)]
   real(dp),protected :: gsv(ng)=sin(gvp/2),gcv(ng)=cos(gvp/2)
   ! (1-cos(v))/X in X=v**2 (highest first), for the near-parabolic coordinate.
   real(dp),parameter :: cser(9)=[1._dp/6402373705728000._dp,-1._dp/20922789888000._dp,1._dp/87178291200._dp, &
      -1._dp/479001600._dp,1._dp/3628800._dp,-1._dp/40320._dp,1._dp/720._dp,-1._dp/24._dp,.5_dp]
   ! Parabolic series of (v-sin(v)cos(v))/sin(v)**3 in v**2 (highest first),
   ! used for |v|<0.3 on zero-revolution arcs where the closed form cancels.
   real(dp),parameter :: pser(10)=[115665628927._dp/2601853536522240000._dp, &
      3821869001._dp/10751460894720000._dp,133541._dp/48117888000._dp,763967._dp/36324288000._dp, &
      1393481._dp/9081072000._dp,1181._dp/1108800,29._dp/4200,17._dp/420,1._dp/5,2._dp/3]
   real(dp),parameter :: pser2(10)=pser/sqrt(2._dp)
   ! Zero-revolution starter tables (tools/gen_tables.py).
   include 'excali_tables.inc'
   ! Transfer frame normalized by the larger radius (see geometry).
   type :: frame
      real(dp) :: scale,rn1,rn2,e1(3),e2(3),t1v(3),t2v(3),cs,ss,s1,s2,p,q,a,b,y0,yend,ypi,t
   end type

   interface
      pure real(c_double) function cbrt_c(x) bind(c,name='cbrt')
         import::c_double
         real(c_double),value,intent(in)::x
      end function
      pure real(c_double) function log1p_c(x) bind(c,name='log1p')
         import::c_double
         real(c_double),value,intent(in)::x
      end function
   end interface

contains

   subroutine excali_solve(r1,r2,tof,mu,direction,nrev,branch,v1,v2,status,plane_normal)
      real(dp),intent(in)::r1(3),r2(3),tof,mu
      integer,intent(in)::direction,nrev,branch
      real(dp),intent(out)::v1(3),v2(3)
      integer,intent(out)::status
      real(dp),optional,intent(in)::plane_normal(3)
      real(dp)::sc,tt,vscale
      integer::half
      v1=ieee_value(0._dp,ieee_quiet_nan);v2=v1;status=-1
      if(.not.all(ieee_is_finite(r1)).or..not.all(ieee_is_finite(r2)))return
      if(.not.ieee_is_finite(mu).or..not.ieee_is_finite(tof).or.mu<=0.or.tof<=0)return
      sc=max(norm2(r1),norm2(r2));if(sc<=0.or..not.ieee_is_finite(sc))return
      vscale=sqrt(mu/sc);tt=tof*vscale/sc
      if(.not.(vscale>tiny(sc).and.vscale<huge(sc)).or..not.(tt>tiny(sc).and.tt<huge(sc)))then
         ! Extreme units: scale through logarithms.
         tt=exp(log(tof)+.5_dp*log(mu)-1.5_dp*log(sc))
         vscale=exp(.5_dp*(log(mu)-log(sc)))
      endif
      if(.not.ieee_is_finite(tt).or..not.ieee_is_finite(vscale).or.tt<=0.or.vscale<=0)return
      call excali_one(r1/sc,r2/sc,tt,direction,v1,v2,status,half,nrev,branch,plane_normal)
      if(status==0.or.status==1)then
         v1=v1*vscale;v2=v2*vscale
      endif
   end subroutine

   subroutine excali_one(r1,r2,tof,direction,v1,v2,status,half,nrev,branch,plane_normal)
   real(dp),intent(in)          :: r1(3),r2(3),tof
   integer,intent(in)           :: direction
   real(dp),intent(out)         :: v1(3),v2(3)
   integer,intent(out)          :: status,half
   integer,optional,intent(in)  :: nrev,branch
   real(dp),optional,intent(in) :: plane_normal(3)
   integer :: m(1),br(1),st(1),hf(1)
   real(dp) :: w1(3,1),w2(3,1)
   m=0;br=-1
   if(present(nrev))m=nrev
   if(present(branch))br=branch
   call solve_chunk(1,r1,r2,[tof],[direction],w1,w2,st,hf,m,br,plane_normal)
   v1=w1(:,1);v2=w2(:,1);status=st(1);half=hf(1)
   end subroutine excali_one

   ! Multiple revolutions off the tables (or whose rectangle fails its
   ! certificate): fixed ratios over regions chosen from the geometry, in a
   ! coordinate free of singularities near the root. Returns the root in
   ! coordinate k with its half-angle functions; status 0/1 solved, 5 below
   ! the minimum time, -4 failed.
   subroutine general_root(g,m,br,k,x,hs,hc,trig,status)
   type(frame),intent(in) :: g
   integer,intent(in)     :: m,br
   integer,intent(out)    :: k,status
   real(dp),intent(out)   :: x,hs,hc
   logical,intent(out)    :: trig
   real(dp) :: rn1,rn2,ss,a,b,y0,yend,ypi,t,lo,hi,tmin,d(0:3),c,dc,tref,s3,s0
   rn1=g%rn1;rn2=g%rn2;ss=g%ss;a=g%a;b=g%b;y0=g%y0;yend=g%yend;ypi=g%ypi;t=g%t
   trig=.false.;status=-4;x=0;hs=0;hc=1;k=1;s0=0
   if(m==0)return             ! zero revolutions all take the table path
   ! Multiple revolutions. Every ellipse through both points has a>=s/2,
   ! so T > m*2*pi*(s/2)**1.5 (within 31% of the minimum time) rejects
   ! most infeasible times from the geometry alone.
   status=5
   dc=sqrt((rn1-rn2)**2+4*rn1*rn2*ss*ss);c=(rn1+rn2+dc)/2     ! chord, semiperimeter
   if(t<m*pi*c*sqrt(c)/sqrt2*(1-1e-12_dp))return
   select case(mrev_global())
    case(0)
      return
    case(1)
      status=-4;tref=tmin
      if(branch_fixed())status=0
    case default
      status=-4                                   ! certificate failed (not observed)
   end select

   contains

      ! Multiple revolutions from two fixed 32-node zero-height ratios over the
      ! whole interval [0,pi] (or the scaled interval below).
      ! 1. v* is the only zero of D=T'*sin(v)**4 on (0,pi); D is pole-free with
      !    D(0)=-3m*pi*y0**1.5/sqrt2, D(pi)=3(m+1)*pi*ypi**1.5/sqrt2. Then
      !    tmin=T(v*) is exact to rounding (T'(v*)=0) and decides feasibility.
      ! 2. rho=sign(v-v*)*sqrt(1-(tmin/T)**(1/3)) rises from -1 to 1 through the
      !    fold and both poles; branch br is its level set
      !    rho=br*sqrt(1-(tmin/t)**(1/3)). T-tmin comes from sum-to-product
      !    differences about v*, so rho keeps full precision near v*.
      ! Each ratio is certified by its second moment. Result: 0 when done
      ! (status set), 1 when v* and tmin are known and the branch is outside
      ! the scaled interval (branch_fixed), 2 when the v* certificate fails.
      integer function mrev_global() result(r)
         real(dp)::tv(ng),t1(ng),ff(ng),fp(ng),xn(ng),wn(ng),sxn(ng),cxn(ng),dd(0:3),s,co,dv,dvp,fa,fb, &
                   vs,err,hss,hcs,ss3,fs,ys,qs,sd,cd,sm,cm,dyv,dsv,dphi,ds3,dfv,yv,qv,dtv,u,ut,rt,rt2, &
                   rho,rho2,den,du,sn,cs2,pv,pt,t2s,lo2,hi2,cen,hl,ve,q, &
                   xn2(ng+2),sx2(ng+2),cx2(ng+2),tv2(ng+2),t12(ng+2),ff2(ng+2),fp2(ng+2),rr(ng+2)
         integer::j,nj
         logical::whole
         r=2
         ! Nodes: the whole interval while the sqrt(y) branch points stay away
         ! from it (|b/(y0+b)|<=0.95); otherwise v* lies between the branch
         ! points' scale and 1, and the interval [ve/2,2ve] around the large-m
         ! minimum ve (cos ve = b/(y0+b+sqrt(y0*ypi))), or its mirror about pi,
         ! keeps every singularity ~1.7 half-widths from the centre.
         whole=abs(b)<=.95_dp*abs(y0+b)
         if(whole)then
            lo2=0;hi2=pi;xn=gv;wn=gw;sxn=gsv;cxn=gcv
         else
            ! For small m the a*sqrt(y) term moves v*: as y0/b -> 0, T'=0
            ! balances 3m*pi*e2/v**3 against 2 (e2=2*y0/b), so v* is at most
            ! (1.5m*pi*e2)**(1/3); towards pi it tends to ~2/((m+1)*pi) instead.
            q=y0+b+sqrt(y0*ypi)
            if(b>0)then
               ve=2*asin(sqrt((y0+sqrt(y0*ypi))/(2*q)))          ! 1-cos ve, exact
               ve=min(ve,rcbrt(3*m*pi*y0/b))
               lo2=ve/2;hi2=min(2*ve,(ve+pi)/2)
            else
               ve=2*asin(sqrt((ypi+sqrt(y0*ypi))/(2*q)))         ! pi - ve
               ve=max(ve,2/((m+1)*pi))
               lo2=max(pi-2*ve,(pi-ve)/2);hi2=pi-ve/2
            endif
            cen=(lo2+hi2)/2;hl=(hi2-lo2)/2
            xn=[cen+hl*ggx,cen-hl*ggx];wn=[ggw*hl,ggw*hl]
            do j=1,ng
               call sincos_half(xn(j),sxn(j),cxn(j))
            enddo
         endif
         do j=1,ng
            call time_derivs(xn(j),1,m,y0,yend,a,b,dd,sxn(j),cxn(j))
            tv(j)=dd(0);t1(j)=dd(1)
            s=2*sxn(j)*cxn(j);co=(cxn(j)-sxn(j))*(cxn(j)+sxn(j))
            dv=dd(1)*s**4;dvp=(dd(2)*s+4*dd(1)*co)*s**3
            ff(j)=1/dv;fp(j)=-dvp*ff(j)*ff(j)
         enddo
         if(whole)then
            fa=-sqrt2/(3*m*pi*y0*sqrt(y0));fb=sqrt2/(3*(m+1)*pi*ypi*sqrt(ypi))
         else
            ! The ends' signs certify that v* is inside.
            call sincos_half(lo2,hss,hcs);call time_derivs(lo2,1,m,y0,yend,a,b,dd,hss,hcs)
            s=2*hss*hcs;fa=1/(dd(1)*s**4)
            call sincos_half(hi2,hss,hcs);call time_derivs(hi2,1,m,y0,yend,a,b,dd,hss,hcs)
            s=2*hss*hcs;fb=1/(dd(1)*s**4)
            if(.not.(fa<0.and.fb>0))return
         endif
         call global_ratio(ng,ff,fp,fa,fb,lo2,hi2,xn,wn,(lo2+hi2)/2,vs,err)
         if(.not.(vs>lo2.and.vs<hi2.and.err<=1e-13_dp*(hi2-lo2)))return
         call sincos_half(vs,hss,hcs)
         call time_derivs(vs,1,m,y0,yend,a,b,dd,hss,hcs)
         tmin=dd(0);t2s=dd(2)
         r=0;status=5
         if(t<tmin*(1-16*epsilon(t)))return
         if(abs(t-tmin)<=16*epsilon(t)*tmin)then
            status=1;k=1;x=vs;hs=hss;hc=hcs;trig=.true.;return
         endif
         status=-4
         ! From here a failure keeps v* and tmin for branch_fixed.
         r=1;x=vs;d=dd;hs=hss;hc=hcs
         lo=lo2;hi=hi2
         if(whole)then
            lo=vs;hi=vs
         endif
         ! Branch target, 1-(tmin/t)**(1/3) in cancellation-free form.
         pt=tmin/t;ut=cbrt_c(pt);rt2=(1-pt)/(1+ut+ut*ut);rt=br*sqrt(rt2)
         ! Values at v*: sin, F, y, sqrt(y).
         s=2*hss*hcs;ss3=s**3
         fs=(vs-s*(hcs-hss)*(hcs+hss)+m*pi)/(sqrt2*ss3)
         ys=y0+2*b*hss*hss;if(b<0)ys=yend-2*b*hcs*hcs
         qs=sqrt(ys)
         ! rho at the nodes; on a scaled interval also at its ends (j=ng+1:
         ! lo2, ng+2: hi2), where the level set must change sign.
         nj=ng
         if(.not.whole)then
            nj=ng+2;xn2(ng+1)=lo2;xn2(ng+2)=hi2
            call sincos_half(lo2,sx2(ng+1),cx2(ng+1));call sincos_half(hi2,sx2(ng+2),cx2(ng+2))
            call time_derivs(lo2,1,m,y0,yend,a,b,dd,sx2(ng+1),cx2(ng+1));tv2(ng+1)=dd(0);t12(ng+1)=dd(1)
            call time_derivs(hi2,1,m,y0,yend,a,b,dd,sx2(ng+2),cx2(ng+2));tv2(ng+2)=dd(0);t12(ng+2)=dd(1)
         endif
         xn2(1:ng)=xn;sx2(1:ng)=sxn;cx2(1:ng)=cxn;tv2(1:ng)=tv;t12(1:ng)=t1
         do j=1,nj
            ! Half-angle functions of d=v-v* and of (v+v*)/2 by angle
            ! addition; sin(d/2) by its series where the addition cancels.
            dv=xn2(j)-vs
            if(abs(dv)<.5_dp)then
               sd=dv/2;sn=sd*sd
               sd=sd*(1-sn/6*(1-sn/20*(1-sn/42*(1-sn/72*(1-sn/110*(1-sn/156))))))
            else
               sd=sx2(j)*hcs-cx2(j)*hss
            endif
            cd=cx2(j)*hcs+sx2(j)*hss
            sm=sx2(j)*hcs+cx2(j)*hss;cm=cx2(j)*hcs-sx2(j)*hss
            s=2*sx2(j)*cx2(j)
            dyv=2*b*sm*sd                                  ! y(v)-y(v*)
            dsv=2*cm*sd                                    ! sin v - sin v*
            cs2=1-2*sm*sm                                  ! cos(v+v*)
            dphi=dv-cs2*2*sd*cd                            ! v-sin v cos v, differenced
            ds3=dsv*(s*s+s*(2*hss*hcs)+(2*hss*hcs)**2)
            dfv=(dphi/sqrt2-fs*ds3)/s**3                   ! F(v)-F(v*)
            yv=ys+dyv;qv=sqrt(yv)
            pv=yv*(fs+dfv)+a
            dtv=dyv/(qv+qs)*pv+qs*(dyv*(fs+dfv)+ys*dfv)    ! T(v)-tmin
            u=cbrt_c(tmin/tv2(j))
            rho2=(dtv/tv2(j))/(1+u+u*u)
            rho=sign(sqrt(max(rho2,0._dp)),dv)
            ! 1/(rho-rt) = (rho+rt)/(rho**2-rt**2); the difference is taken
            ! where it is accurate: in u near the poles, in rho**2 near v*.
            den=rho2-rt2
            if(u<.5_dp)den=ut-u
            du=-u*t12(j)/(3*tv2(j))
            ff2(j)=(rho+rt)/den;rr(j)=rho
            if(rho/=0)then
               fp2(j)=(-du/(2*rho))/den+(rho+rt)*du/(den*den)
            else
               fp2(j)=sqrt(max(t2s,0._dp)/(6*tmin))/den     ! d rho/dv at v*
            endif
         enddo
         if(whole)then
            fa=1/(-1-rt);fb=1/(1-rt)
         else
            if(.not.(rr(ng+1)<rt.and.rt<rr(ng+2)))return  ! root outside: starter path
            fa=ff2(ng+1);fb=ff2(ng+2)
         endif
         call global_ratio(ng,ff2,fp2,fa,fb,lo2,hi2,xn,wn,(lo2+hi2)/2,vs,err)
         if(.not.(err<=1e-13_dp*(hi2-lo2).and.(vs-x)*br>0.and.vs>lo2.and.vs<hi2))return
         ! Near either pole recentre the moment there, so the small distance
         ! keeps its relative precision.
         k=1
         if(vs<.25_dp)then
            call global_ratio(ng,ff2,fp2,fa,fb,lo2,hi2,xn,wn,0._dp,vs,err)
         elseif(vs>pi-.25_dp)then
            call global_ratio(ng,ff2,fp2,fa,fb,lo2,hi2,xn,wn,pi,vs,err)
            vs=pi-vs;k=2
         endif
         if(.not.(vs>0))then
            k=1;return
         endif
         x=vs;call sincos_half(x,hs,hc);trig=.true.
         status=0;r=0
      end function

      ! Multi-revolution branch outside the minimum's interval [lo,hi] (lo=hi=v*
      ! after the whole-interval ratio): sigma=(tmin/T)**(1/3) from
      ! the pole to the interval, in zeta about the pole when the sqrt(y)
      ! branch points crowd it (kinds 4, 3), else in v or d (kinds 5, 2).
      logical function branch_fixed() result(ok)
         real(dp)::lo2,hi2
         integer::kz
         s3=cbrt_c(tref/t);lo2=0
         if(br<0)then
            if(b>0.and.lo<hi)then
               kz=4;s0=sqrt(y0/(2*b));hi2=asinh(sin(lo/2)/s0)
            else
               kz=5;hi2=lo
            endif
         else
            if(b<0.and.lo<hi)then
               kz=3;s0=sqrt(yend/(-2*b));hi2=asinh(sin((pi-hi)/2)/s0)
            else
               kz=2;hi2=pi-hi
            endif
         endif
         ok=ratio_root(kz,lo2,hi2)
      end function

      ! One fixed 24-node ratio in coordinate kz on [l0,h0], certified by
      ! the signs at the ends and the second moment; sets k, x, hs, hc.
      ! In zeta the singularities lie at height pi/2 while the interval grows
      ! like log(1/s0), so it is cut into equal pieces of length <= 4 and the
      ! piece holding the sign change is taken.
      logical function ratio_root(kz,l0,h0) result(ok)
         integer,intent(in)::kz
         real(dp),intent(in)::l0,h0
         real(dp)::zn(nz),wn(nz),ff(nz),fp(nz),fa,fb,z,err,tz,dz,cen,hl,lo2,hi2
         integer::j,np
         ok=.false.
         lo2=l0;hi2=h0;np=1
         if(kz==3.or.kz==4)np=ceiling((h0-l0)/4)
         hl=(h0-l0)/np
         do j=1,np-1
            call zeval(kz,l0+j*hl,fa,dz,tz)
            if(fa>0)then
               hi2=l0+j*hl;exit
            endif
            lo2=l0+j*hl
         enddo
         cen=(lo2+hi2)/2;hl=(hi2-lo2)/2
         if(.not.(hl>0))return
         zn=[cen+hl*gzx,cen-hl*gzx];wn=[gzw*hl,gzw*hl]
         do j=1,nz
            call zeval(kz,zn(j),ff(j),fp(j),tz)
         enddo
         call zeval(kz,lo2,fa,dz,tz);call zeval(kz,hi2,fb,dz,tz)
         if(.not.(fa<0.and.fb>0))return
         call global_ratio(nz,ff,fp,fa,fb,lo2,hi2,zn,wn,cen,z,err)
         ! A root near the lower end keeps its relative precision with the
         ! moment taken about that end.
         if(z-lo2<.25_dp*(hi2-lo2))call global_ratio(nz,ff,fp,fa,fb,lo2,hi2,zn,wn,lo2,z,err)
         if(.not.(z>lo2.and.z<hi2.and.err<=1e-10_dp*(hi2-lo2)))return
         if(.not.(z>0))return
         trig=.true.
         select case(kz)
          case(2)
            k=2;x=z;call sincos_half(x,hs,hc)
          case(3,4)
            call shch(z,tz,dz);hs=s0*tz;hc=sqrt((1-hs)*(1+hs))
            k=1;if(kz==3)k=2
            x=2*asin(hs)
          case(5)
            k=1;x=z;call sincos_half(x,hs,hc)
         end select
         ok=.true.
      end function

      ! f, f' and T at z in coordinate kz (multi-revolution branches), for
      ! f = 1/g with g = sigma-s3 increasing, sigma=(tref/T)**(1/3), which is
      ! zero at a pole (z=0 there): 2 d=pi-v, 3 zeta about pi
      ! (sin(d/2)=s0*sinh), 4 zeta about 0 (sin(v/2)=s0*sinh), 5 v.
      subroutine zeval(kz,z,f,fp,tt)
         integer,intent(in)::kz
         real(dp),intent(in)::z
         real(dp),intent(out)::f,fp,tt
         real(dp)::u1,u2,dt,sg,dd(0:3),sh,ch
         if(z==0)then
            tt=huge(z);f=-1/s3;fp=0;return               ! pole: sigma=0
         endif
         select case(kz)
          case(2)
            call sincos_half(z,sh,ch);call time_derivs(z,2,m,y0,yend,a,b,dd,sh,ch)
            dt=dd(1)
          case(3,4)
            call shch(z,sh,ch);u1=s0*sh;u2=sqrt((1-u1)*(1+u1))
            if(kz==3)then
               call time_derivs(2*asin(u1),2,m,y0,yend,a,b,dd,u1,u2)
            else
               call time_derivs(2*asin(u1),1,m,y0,yend,a,b,dd,u1,u2)
            endif
            dt=dd(1)*2*s0*ch/u2
          case default
            call sincos_half(z,sh,ch);call time_derivs(z,1,m,y0,yend,a,b,dd,sh,ch)
            dt=dd(1)
         end select
         tt=dd(0);sg=cbrt_c(tref/tt)
         f=1/(sg-s3);fp=sg*dt/(3*tt)*f*f                ! -dsigma*f**2
      end subroutine

   end subroutine general_root







   ! Batch of normalized solves, same arguments per case as excali_one, in
   ! chunks of nb cases through solve_chunk (the code excali_one runs with
   ! one case), so the results are identical.
   subroutine excali_batch(n,r1,r2,tof,direction,v1,v2,status,half,nrev,branch)
   integer,intent(in)   :: n
   real(dp),intent(in)  :: r1(3,n),r2(3,n),tof(n)
   integer,intent(in)   :: direction(n),nrev(n),branch(n)
   real(dp),intent(out) :: v1(3,n),v2(3,n)
   integer,intent(out)  :: status(n),half(n)
   integer :: i0,nc
   ! Chunks are independent; under OpenMP they run in parallel (dynamic
   ! schedule), and small batches stay serial.
   !$omp parallel do schedule(dynamic) private(nc) if(n>8*nb)
   do i0=1,n,nb
      nc=min(nb,n-i0+1)
      call solve_chunk(nc,r1(:,i0:i0+nc-1),r2(:,i0:i0+nc-1),tof(i0:i0+nc-1),direction(i0:i0+nc-1), &
         v1(:,i0:i0+nc-1),v2(:,i0:i0+nc-1),status(i0:i0+nc-1),half(i0:i0+nc-1),nrev(i0:i0+nc-1),branch(i0:i0+nc-1))
   enddo
   !$omp end parallel do
   end subroutine excali_batch

   ! nc <= nb normalized solves, stage by stage so independent cases overlap:
   ! geometry; table or closed-form starters; one rectangle (3 or 8 nodes) at a
   ! fixed fraction of a lower bound on the distance to every singularity;
   ! general_root for multiple revolutions off the tables; velocities. Every
   ! zero-revolution case takes exactly one path.
   subroutine solve_chunk(nc,r1c,r2c,tofc,dirc,v1c,v2c,statc,halfc,nrev,branch,pnormal)
   integer,intent(in)   :: nc
   real(dp),intent(in)  :: r1c(3,nc),r2c(3,nc),tofc(nc)
   integer,intent(in)   :: dirc(nc),nrev(nc),branch(nc)
   real(dp),intent(out) :: v1c(3,nc),v2c(3,nc)
   integer,intent(out)  :: statc(nc),halfc(nc)
   real(dp),optional,intent(in) :: pnormal(3)
   type(frame) :: gc(nb)
   real(dp) :: xc(nb),cc(nb),hsc(nb),hcc(nb),wbc(nb),lbc(nb),ubc(nb),ddc(nb),elc,vsc,chord,semi
   integer :: kc(nb),kd(nb),krc(nb),accc(nb),j
   logical :: fast(nb),nearc(nb),trig(nb),eightc(nb)
   ! Stage 1: geometry; the multi-revolution tables.
   do j=1,nc
      fast(j)=.false.;nearc(j)=.false.;halfc(j)=0;statc(j)=-1
      if(abs(dirc(j))/=1.or.nrev(j)<0.or.abs(branch(j))/=1.or..not.(tofc(j)>0) &
         .or..not.ieee_is_finite(tofc(j)))cycle
      call geometry(r1c(:,j),r2c(:,j),tofc(j),dirc(j),gc(j),statc(j),halfc(j),pnormal)
      if(statc(j)/=-4)cycle
      lbc(j)=0;ubc(j)=pi;accc(j)=0
      if(nrev(j)>0)then
         ! Every ellipse through both points has a>=s/2, so T > m*2*pi*(s/2)**1.5
         ! (within 31% of the minimum time) rejects most infeasible times.
         chord=sqrt((gc(j)%rn1-gc(j)%rn2)**2+4*gc(j)%rn1*gc(j)%rn2*gc(j)%ss**2);semi=(gc(j)%rn1+gc(j)%rn2+chord)/2
         statc(j)=5
         if(gc(j)%t>=nrev(j)*pi*semi*sqrt(semi)/sqrt2*(1-1e-12_dp)) &
            call mrev_table(gc(j)%a/(sqrt2*semi)/(1+sqrt(chord/semi)),semi,gc(j)%t,nrev(j),branch(j),xc(j),vsc,statc(j))
         if(statc(j)==2)call mrev_fold(gc(j),nrev(j),branch(j),vsc,xc(j),statc(j))
         fast(j)=statc(j)==1;kc(j)=1
         if(statc(j)==3)then
            statc(j)=1;xc(j)=vsc                    ! t = t_min: the root is v* itself
         elseif(statc(j)/=5)then
            statc(j)=-4
         endif
         if(branch(j)<0)then
            ubc(j)=vsc
         else
            lbc(j)=vsc
         endif
         if(xc(j)>.75_dp*pi)then                    ! distance from the upper pole
            elc=lbc(j);lbc(j)=pi-ubc(j);ubc(j)=pi-elc;xc(j)=pi-xc(j);kc(j)=2
         endif
         if(statc(j)==1)call sincos_half(xc(j),hsc(j),hcc(j))
      endif
   enddo
   ! Stage 2: the zero-revolution region selector (its own loop: inside the
   ! geometry loop it costs 5-10 ns per solve in code generation), then the
   ! rectangle's coordinate, centre and half-width.
   do j=1,nc
      if(statc(j)/=-4.or.nrev(j)>0)cycle
      call zero_rev_start(gc(j),kc(j),xc(j),wbc(j),accc(j),nearc(j))
      fast(j)=xc(j)>0.or.kc(j)==0
   enddo
   do j=1,nc
      if(.not.fast(j).or.kc(j)==0)cycle
      call fast_coordinate(kc(j),xc(j),nrev(j),gc(j)%b,wbc(j),lbc(j),ubc(j),accc(j),kd(j),cc(j),elc,eightc(j),hsc(j),hcc(j))
      fast(j)=elc>0;krc(j)=kexp(elc)
   enddo
   ! Stage 3: one rectangle each, then its acceptance.
   do j=1,nc
      ddc(j)=0
      if(fast(j).and.kc(j)/=0.and.accc(j)>=0)call rect_ratio(kd(j),cc(j),hsc(j),hcc(j),krc(j),eightc(j),gc(j)%y0,merge(wbc(j),gc(j)%yend,kd(j)==-2),gc(j)%a,gc(j)%b, &
         gc(j)%t,nrev(j),ddc(j))
   enddo
   do j=1,nc
      if(fast(j).and.kc(j)/=0)call fast_accept(kd(j),cc(j),ddc(j),hsc(j),hcc(j),kc(j),xc(j),fast(j))
   enddo
   ! Stage 4: the multi-revolution region formulas for the rest, then velocities.
   do j=1,nc
      trig(j)=.true.
      if(fast(j))then
         statc(j)=0
         if(kc(j)==0)then
            kc(j)=1;xc(j)=0;hsc(j)=0;hcc(j)=1              ! parabola
         endif
      elseif(statc(j)==-4)then
         call general_root(gc(j),nrev(j),branch(j),kc(j),xc(j),hsc(j),hcc(j),trig(j),statc(j))
         if(.not.trig(j))then
            if(kc(j)>0)then
               call sincos_half(xc(j),hsc(j),hcc(j))
            else
               call shch(xc(j)/2,hsc(j),hcc(j))
            endif
         endif
      endif
      if(statc(j)==0.or.statc(j)==1)then
         call velocities(gc(j),r1c(:,j),r2c(:,j),dirc(j),kc(j),xc(j),hsc(j),hcc(j),nearc(j),v1c(:,j),v2c(:,j))
         if(.not.all(ieee_is_finite(v1c(:,j))).or..not.all(ieee_is_finite(v2c(:,j))))statc(j)=-4
      endif
      if(statc(j)/=0.and.statc(j)/=1)then
         v1c(:,j)=ieee_value(0.0_dp,ieee_quiet_nan);v2c(:,j)=v1c(:,j)
      endif
   enddo

   ! The stages' steps, one call site each (internal, so the compiler inlines them).
   contains

      ! Transfer frame, normalized by the larger radius: unit vectors, the
      ! transverse directions, |cos(theta/2)|, |sin(theta/2)|, s1=sqrt(r1),
      ! s2=sqrt(r2), q=s1-s2 (cancellation-free), a, b, y0=y(v=0), yend and
      ! ypi=y(pi), where y(v)=y0+2b*sin(v/2)**2 (b>=0) or yend-2b*cos(v/2)**2
      ! (b<0), and the normalized time t. status -4 means usable.
      pure subroutine geometry(r1,r2,tof,direction,g,status,half,plane_normal)
         real(dp),intent(in)::r1(3),r2(3),tof
         integer,intent(in)::direction
         type(frame),intent(out)::g
         integer,intent(out)::status,half
         real(dp),optional,intent(in)::plane_normal(3)
         real(dp)::h(3),h2(3),hn,dr
         half=0;status=-1
         g%rn1=sqrt(dot_product(r1,r1));g%rn2=sqrt(dot_product(r2,r2))
         if(.not.(g%rn1>0.and.g%rn2>0).or..not.ieee_is_finite(g%rn1+g%rn2))return
         dr=g%rn1-g%rn2
         if(abs(dr)<.5_dp*max(g%rn1,g%rn2))dr=dot_product(r1-r2,r1+r2)/(g%rn1+g%rn2)  ! without cancellation
         ! Correctly rounded unit vectors; the larger normalized radius is exactly 1.
         g%e1=r1/g%rn1;g%e2=r2/g%rn2
         if(g%rn1>=g%rn2)then
            g%scale=g%rn1;g%rn2=g%rn2/g%rn1;g%rn1=1
         else
            g%scale=g%rn2;g%rn1=g%rn1/g%rn2;g%rn2=1
         endif
         h=cross(g%e1,g%e2);hn=sqrt(dot_product(h,h))
         status=-2
         if(hn<32*epsilon(hn))then
            half=-1
            if(.not.present(plane_normal))return
            if(norm2(g%e1+g%e2)>32*epsilon(hn))then
               status=-3;return
            endif
            if(.not.all(ieee_is_finite(plane_normal)))return
            h=plane_normal;hn=norm2(h)
            if(hn<=0.or..not.ieee_is_finite(hn))return
            if(max(abs(dot_product(h/hn,g%e1)),abs(dot_product(h/hn,g%e2)))>64*epsilon(hn))return
            half=1
         endif
         if(hn<1e-4_dp)half=1
         h=h*(direction/hn)
         h2=g%e1+g%e2;g%cs=sqrt(dot_product(h2,h2))/2;h2=g%e2-g%e1;g%ss=sqrt(dot_product(h2,h2))/2
         g%t1v=cross(h,g%e1);g%t2v=cross(h,g%e2)
         g%s1=sqrt(g%rn1);g%s2=sqrt(g%rn2)
         g%p=g%s1*g%s2;g%a=direction*sqrt2*g%p*g%cs;g%b=sqrt2*g%a
         g%q=(g%rn1-g%rn2)/(g%s1+g%s2)                  ! q=s1-s2
         if(min(g%rn1,g%rn2)>.5_dp)g%q=(dr/g%scale)/(g%s1+g%s2)
         g%yend=g%q*g%q+2*g%p*g%ss*g%ss/(1+g%cs)
         g%y0=g%yend;if(direction<0)g%y0=g%q*g%q+2*g%p*(1+g%cs)
         g%ypi=g%yend;if(g%b>=0)g%ypi=g%y0+2*g%b
         g%t=tof/(g%scale*sqrt(g%scale))
         status=-4
      end subroutine

      ! Multiple revolutions near the fold (|sigma| < 1e-2, or t within the
      ! table's t_min margin), where the table's v* error (6e-5 of
      ! R=min(v*,pi-v*)) would dominate |v-v*|. v* is the only zero of
      ! D = T'*sin(v)**4 near the table's value: one 8-node ratio at R/8. Then
      ! t_min=T(v*) is exact to rounding and decides feasibility (stat 5; stat
      ! 3: t = t_min, the root is v*), and the branch root starts from the fold
      ! series v = v* + k*s*(1 - T3*k*s/(6*T2)), k = sqrt(6 t_min/T2) with T2,
      ! T3 the second and third derivatives at v*, s = branch*sqrt(1-(t_min/t)**(1/3)),
      ! within 0.6% of |v-v*| for |s| < 3e-2 (stat 1, vs = v*).
      pure subroutine mrev_fold(g,m,br,vs,x,stat)
         type(frame),intent(in)::g
         integer,intent(in)::m,br
         real(dp),intent(inout)::vs
         real(dp),intent(out)::x
         integer,intent(out)::stat
         real(dp)::el,xn(n1),ff(n1),fp(n1),dd(0:3),sn,cn,s,co,fa,fb,err,pt,ut,sg,k,z
         integer::i
         el=pow2m(kexp(min(vs,pi-vs)/8));fa=0;fb=0
         do i=1,n1+2
            z=vs+el*u8(i)
            call sincos_half(z,sn,cn)
            call time_derivs(z,1,m,g%y0,g%yend,g%a,g%b,dd,sn,cn)
            s=2*sn*cn;co=(cn-sn)*(cn+sn)
            if(i<=n1)then
               xn(i)=z;ff(i)=1/(dd(1)*s**4);fp(i)=-(dd(2)*s+4*dd(1)*co)*s**3*ff(i)**2
            elseif(i==n1+1)then
               fb=1/(dd(1)*s**4)
            else
               fa=1/(dd(1)*s**4)
            endif
         enddo
         call global_ratio(n1,ff,fp,fa,fb,vs-el,vs+el,xn,w8(1:n1)*el,vs,x,err)
         vs=x
         call sincos_half(vs,sn,cn)
         call time_derivs(vs,1,m,g%y0,g%yend,g%a,g%b,dd,sn,cn)
         stat=5
         if(g%t<dd(0)*(1-16*epsilon(dd(0))))return
         stat=3
         if(abs(g%t-dd(0))<=16*epsilon(dd(0))*dd(0))return
         pt=dd(0)/g%t;ut=cbrt_c(pt);sg=br*sqrt((1-pt)/(1+ut+ut*ut))
         k=sqrt(6*dd(0)/dd(2))
         x=vs+k*sg*(1-dd(3)*k*sg/(6*dd(2)))
         stat=1
      end subroutine

      ! Zero-revolution region selector: the coordinate k, the starter x in it
      ! and the starter's accuracy class acc, which fixes the rectangle before
      ! it is evaluated (fast_coordinate). wb is the short-way branch point;
      ! near marks the near-parabolic window (velocities use Lagrange f,g).
      !
      !   region                         k      starter                     acc
      !   parabola                       0      x = 0                        -
      !   |t-tp| <= 3.1% tp (near)       +-1    series reversion in X        X
      !     2nd-order term < 1e-8 of X          the series itself            -1
      !   ellipse, |tau| < 0.75          1, 2   table_v (tab)                0
      !   ellipse, -0.9 <= tau < -0.75   1, 2   table_v (tab)                1
      !   ellipse, tau >= 0.75           1, 2   table_vp (tabp)              0
      !     tau > 1-1e-6                        tabp clamped                 2
      !   ellipse, tau < -0.9            1, 2   long_d (closed form)         2
      !   hyperbola, t/tp >= 0.1         -1     table_w (tabh)               0
      !   hyperbola, t/tp < 0.1          -1,-2  large-w closed form          1
      !
      ! acc 0: within 0.2% of the singularity bound R, 3 nodes at R/256;
      ! acc 1: within 4% (1.1% for hyperbolas), 8 nodes at R/8; acc 2: within
      ! 6.4%, 8 nodes at R/4; X: 8 nodes at |X|/8 (or less). k=2 is pi-v
      ! (used when v > 3pi/4, or d < pi/4 from long_d); k=-2 is xi=wb-w.
      pure subroutine zero_rev_start(g,k,x,wb,acc,near)
         type(frame),intent(in)::g
         integer,intent(out)::k,acc
         real(dp),intent(out)::x,wb
         logical,intent(out)::near
         real(dp)::tp,u,tau,dc,sl,sw,r0,r1,dt
         tp=sqrt(g%y0)*(g%y0*sqrt2/3+g%a)
         k=0;x=0;wb=0;acc=0;near=.false.
         if(abs(g%t-tp)<=4*epsilon(tp)*tp)return    ! parabola
         near=abs(g%t-tp)<=.031_dp*tp      ! meets the table at phi=0.99
         k=1;u=rcbrt(tp/g%t)
         if(g%t<tp)then
            k=-1;u=g%t/tp
            if(g%b>0)then
               sw=sqrt(g%y0/(2*g%b))
               wb=2*log1p_c(sw+sw*sw/(1+sqrt(1+sw*sw)))
            endif
         endif
         if(near)then
            ! When the second-order term is below 1e-8 of X the series is exact
            ! to rounding (the dropped term is its square), while T-t at a
            ! rectangle's nodes is then only a few hundred ulps and can round
            ! to 0: the series is the root (acc -1).
            call near_parabola_x(g%y0,g%a,g%b,g%t,tp,u,dt)
            if(u*k>0)x=sqrt(abs(u))
            acc=merge(-1,1,abs(dt)<1e-8_dp)
            return
         endif
         if(k<0.and.u<tabh_u0)then
            ! Hyperbola, t/tp<0.1: large-w closed form T = sqrt2*c*sqrt(y)/E,
            ! E=exp(w), c=y0+b, y=c-b*E/2, a quadratic in E (rationalized on the
            ! short way): w=log(E), or on the short way with y=(t*E/(sqrt2*c))**2
            ! < y0/2 (xi < log 2) xi=wb-w (k=-2) from y, exactly and without
            ! cancellation (cosh(w)=(c-y)/b). Within 1.1% of the bound for every tau.
            sw=g%y0+g%b;dc=sqrt((sw*sw*g%b)**2+8*g%t*g%t*sw**3)
            if(g%b>=0)then
               dc=4*sw**3/(sw*sw*g%b+dc)
            else
               dc=(dc-sw*sw*g%b)/(2*g%t*g%t)
            endif
            x=log(dc);acc=1
            if(g%b>0)then
               sl=(g%t*dc/(sqrt2*sw))**2                  ! y
               if(sl<g%y0/2)then
                  r0=sqrt(g%y0*g%ypi);r1=sqrt((g%y0-sl)*(g%ypi-sl))
                  k=-2;x=log1p_c((sl+sl*(2*sw-sl)/(r0+r1))/(sw-sl+r1))
               else
                  x=min(x,.99_dp*wb)
               endif
            endif
            return
         endif
         dc=sqrt(g%q*g%q*(g%s1+g%s2)**2+4*g%p*g%p*g%ss*g%ss);sl=(g%rn1+g%rn2+dc)/2   ! chord, semiperimeter
         tau=g%a/(sqrt2*sl)/(1+sqrt(dc/sl))       ! Lancaster-Blanchard lam/(1+sqrt(1-lam**2))
         if(k<0)then
            x=table_w(tau,u,wb)
            return
         endif
         if(tau>=.75_dp)then
            ! Short-way table tabp, split at Th=T(pi/2): s1=(Th/t)**(1/3) with
            ! scale (K/Th)**(1/3) (side 1), or s2=(t-tp)/(Th-tp) with scale
            ! s2*(s2+2*tp/(Th-tp)) (side 2).
            sw=g%y0+g%b
            sw=sqrt(sw)*(sw*(pi/(2*sqrt2))+g%a)          ! Th
            if(g%t>=sw)then
               u=rcbrt(sw/g%t);x=table_vp(tau,u,rcbrt(pi*g%ypi*sqrt(g%ypi)/(sqrt2*sw)),1)
            else
               u=(g%t-tp)/(sw-tp);x=table_vp(tau,u,u*(u+2*tp/(sw-tp)),2)
            endif
            if(tau>1-1e-6_dp)acc=2
         elseif(tau<-.9_dp)then
            ! Long way toward tau -> -1: th=t/|b|**1.5, eta**2=2*yend/|b|; d=pi-v.
            x=long_d(g%t/((-g%b)*sqrt(-g%b)),sqrt(2*g%yend/(-g%b)));acc=2
            if(x>.25_dp*pi)then
               x=pi-x
            else
               k=2
            endif
            return
         else
            ! Elliptic table: scale (K/tp)**(1/3), K=pi*y(pi)**1.5/sqrt2, for
            ! s<=1/2, else (t-tp)/T0', T0'=dT/d(v**2) at the parabola.
            if(u<=.5_dp)then
               dt=rcbrt(pi*g%ypi*sqrt(g%ypi)/(sqrt2*tp))
            else
               sw=sqrt(g%y0)
               dt=(g%t-tp)/(sw*(g%y0*(rsqrt2/5)+g%b*(rsqrt2/2))+g%a*g%b/(4*sw))
            endif
            x=table_v(tau,u,dt)
            if(abs(tau)>=.75_dp)acc=1
         endif
         if(x>.75_dp*pi)then
            x=pi-x;k=2
         endif
      end subroutine

      ! Coordinate, half-width el and centre half-angle functions of the
      ! rectangle. Every singularity of 1/(T-t) lies at least R from the root:
      ! min(v,pi-v) for ellipses (sqrt(y) branch points, mirror root, pole),
      ! min(w,pi) or min(2w,wb-w,pi) for hyperbolas in w, min(xi,2w,pi) in
      ! xi=wb-w, |X| near the parabola, the bracket [lo,hi] for m>0. el is
      ! R/256 (3 nodes), R/8 or R/4 (8 nodes) by the zero-revolution accuracy
      ! class acc; X 8 nodes at |X|/8 (or less); multiple revolutions 8 nodes
      ! at R/16. kexp rounds el down to a power of two (el<=1 keeps the node
      ! tables in range).
      pure subroutine fast_coordinate(k,x,m,b,wb,lo,hi,acc,kind1,c,el,eight,hs,hc)
         integer,intent(in)::k,m,acc
         real(dp),intent(in)::x,b,wb,lo,hi
         integer,intent(out)::kind1
         real(dp),intent(out)::c,el,hs,hc
         logical,intent(out)::eight
         hs=0;hc=1;kind1=k;c=x;eight=.true.
         if(m>0)then
            el=min(c-lo,hi-c)/8
            call sincos_half(c,hs,hc)
            return
         elseif(abs(k)==1.and.x<.3_dp)then
            kind1=0;c=k*x*x
            el=min(abs(c),.2_dp-abs(c))
            if(k>0)el=min(2*abs(c),.2_dp-abs(c))
            if(k<0.and.b>0)el=min(el,wb*wb-abs(c))
            el=el/16
            return
         elseif(k>0)then
            el=min(c,pi-c)                   ! k=2: x is already pi-v
            call sincos_half(c,hs,hc)
         elseif(k==-2)then
            el=min(c,2*(wb-c),pi)            ! xi: branch point, mirror root, w=+-i*pi
            call shch(c/2,hs,hc)
         else
            ! w: mirror root -w (short way) or w=+-i*pi, the branch point wb, and
            ! at large w the near-copies of the root at w+-2*pi*i (T ~ exp(-w)).
            el=min(c,pi)
            if(b>0)el=min(2*c,wb-c,pi)
            call shch(c/2,hs,hc)
         endif
         select case(acc)
          case(0)
            el=el/256;eight=.false.
          case(1)
            el=el/8
          case default
            el=el/4
         end select
      end subroutine

      ! The rectangle's root is final: move the root and its half-angle
      ! functions to it (ok false only for a non-finite result).
      pure subroutine fast_accept(kind1,c,dd,hs,hc,k,x,ok)
         integer,intent(in)::kind1
         real(dp),intent(in)::c,dd
         real(dp),intent(inout)::hs,hc
         integer,intent(inout)::k
         real(dp),intent(inout)::x
         logical,intent(out)::ok
         real(dp)::cn
         ok=abs(dd)<huge(dd)
         if(.not.ok)return
         cn=c+dd
         if(kind1==0)then
            k=1;if(cn<0)k=-1
            x=sqrt(abs(cn))
            if(k>0)then
               call sincos_half(x,hs,hc)
            else
               call shch(x/2,hs,hc)
            endif
         else
            call rotate(kind1,cn,dd/2,hs,hc)
            k=kind1;x=cn
         endif
         ok=x>0
      end subroutine

      ! Velocities from the root x in coordinate k with half-angle functions
      ! sw, dw: radial and transverse components, cancellation-free, or near the parabola and the
      ! short-way branch point (fine) Lagrange f,g in double-double.
      pure subroutine velocities(g,r1,r2,direction,k,x,sw,dw,fine,v1,v2)
         type(frame),intent(in)::g
         real(dp),intent(in)::r1(3),r2(3),x,sw,dw
         integer,intent(in)::direction,k
         logical,intent(in)::fine
         real(dp),intent(out)::v1(3),v2(3)
         real(dp)::p,c,s,y,dc,sh,ch,vr1,vr2,vt1,vt2
         logical::fg
         p=sign(1.0_dp,g%a);fg=fine;sh=sw;ch=dw
         if(k>0)then
            c=(ch-sh)*(ch+sh);s=(2*sh*ch)**2
         else
            c=1+2*sh*sh;s=-(2*sh*ch)**2
         endif
         if(k==2)then
            c=-c
            y=g%y0+2*g%b*ch*ch
            if(g%b<0)y=g%yend-2*g%b*sh*sh
         elseif(k>0)then
            if(g%b>=0)then
               y=g%y0+2*g%b*sh*sh
            else
               y=g%yend-2*g%b*ch*ch
            endif
         elseif(k==-2)then
            ! From xi = wb-w: sinh(w/2), cosh(w/2) and y in cancellation-free form.
            p=sqrt(g%y0/(2*g%b));c=sqrt(1+p*p)                  ! sinh, cosh of wb/2
            y=2*g%b*(2*p*c*ch-(1+2*p*p)*sh)*sh
            dc=p*ch-c*sh;ch=c*ch-p*sh;sh=dc
            c=1+2*sh*sh;s=-(2*sh*ch)**2;fg=.true.
            p=sign(1.0_dp,g%a)
         else
            y=g%y0-2*g%b*sh*sh
            ! Short way near the branch point, y cancels; T = sqrt(y)*(y*F+a) = t
            ! gives it with that error reduced by y*F/(y*F+a) when a
            ! dominates, and the half-angle functions follow from it in
            ! cancellation-free form.
            if(g%b>0.and.y<g%y0/4)then
               call fhyp_d(x,2*sh*ch,c,dc,s);dc=y*dc
               if(dc<g%a)then
                  y=(g%t/(dc+g%a))**2
                  sh=sqrt((g%y0-y)/(2*g%b));ch=sqrt(1+sh*sh);c=1+2*sh*sh;s=-(2*sh*ch)**2
                  fg=.true.
               endif
            endif
         endif
         dc=p*g%cs-c
         if(p*g%cs*c>0.and.abs(p*g%cs+c)>0.5_dp)dc=(s-g%ss*g%ss)/(p*g%cs+c)
         if(fg.and.g%cs>.3_dp)then
            ! Lagrange f,g from the exact inputs: correctly rounded for this y
            ! (radial/transverse needs the rounded unit vectors).
            call fg_velocities(r1,r2,y,g%scale,direction,v1,v2)
         else
            ! With q=s1-s2: (r2-r1)/(r1*(s2/s1+1)) = -q/s1, and vt uses s2/s1.
            y=sqrt2/sqrt(y*g%scale)                ! common factor sqrt(2/(y*scale))
            dc=y*dc;p=y*p*g%cs*g%q;s=y*g%ss
            sh=1/g%s1;ch=1/g%s2
            vr1=dc-p*sh;vt1=s*g%s2*sh
            vr2=-(dc+p*ch);vt2=s*g%s1*ch
            v1=vr1*g%e1+vt1*g%t1v
            v2=vr2*g%e2+vt2*g%t2v
         endif
      end subroutine
   end subroutine solve_chunk

   ! Zero-height thin rectangle on [c-L,c+L], L=2**(-kr): Gauss-Legendre
   ! 8 or 3 nodes (eight) plus the two endpoints. dd=Q1/Q0 locates the
   ! enclosed root; de=|Q2/Q0-dd**2|/L vanishes for a single pole and
   ! estimates the error; sc reports a sign change of T-t.
   ! Nodes (kind 0: X=k*x**2 series; 1: v; 2: pi-v; -1: w; -2: xi=wb-w) are
   ! evaluated in one fused loop per coordinate; hs,hc are the half-angle
   ! functions of c and the tables those of L*u, combined by angle addition.
   ! Zero-height ratio over [lo,hi] with an n-node rule (nodes xn, weights
   ! wn including the half-width): root = c + Q[((v-c)f)']/Q[f'] for the one
   ! simple pole of f inside, where Q[g'] = g(hi)-g(lo)-sum wn*g'. f, fp: f
   ! and f' at the nodes; fa, fb: f at lo and hi; c: moment centre (any; one
   ! near the root keeps its precision). err = |Q2/Q0-(Q1/Q0)**2|/(hi-lo)
   ! vanishes for a pure pole.
   pure subroutine global_ratio(n,f,fp,fa,fb,lo,hi,xn,wn,c,root,err)
      integer,intent(in)::n
      real(dp),intent(in)::f(n),fp(n),fa,fb,lo,hi,xn(n),wn(n),c
      real(dp),intent(out)::root,err
      real(dp)::q0,q1,q2,e(n),d
      e=xn-c
      q0=fb-fa-sum(wn*fp)
      q1=(hi-c)*fb-(lo-c)*fa-sum(wn*(f+e*fp))
      q2=(hi-c)**2*fb-(lo-c)**2*fa-sum(wn*e*(2*f+e*fp))
      d=q1/q0;root=c+d
      err=abs(q2/q0-d*d)/(hi-lo)
   end subroutine

   pure subroutine rect_ratio(kind1,c,hs,hc,kr,eight,y0,yend,a,b,t,m,dd)
      integer,intent(in)::kind1,kr,m
      real(dp),intent(in)::c,hs,hc,y0,yend,a,b,t
      logical,intent(in)::eight
      real(dp),intent(out)::dd
      real(dp),dimension(n1+2)::uu,wq,ct,st,f,f1
      real(dp)::el,q0,q1,r0,r1,fa,fb,x,sh,ch,sn,cn,rs,yv,y1,ff,ff1,q,g,g1,ser,ser1,pk,y00,ab,mt,sb,cb
      integer::np,j,n
      logical::lw,us
      el=pow2m(kr)
      if(kr>kmax)then
         ! Below 2**-62 (xi ~ t**2 at the shortest times): cos=cosh=1 and
         ! sin=sinh=L*u/2 to rounding.
         if(eight)then
            np=n1+2;uu=u8;wq=w8
         else
            np=5;uu(1:5)=u3;wq(1:5)=w3
         endif
         ct=1;st=el*uu/2
      elseif(eight)then
         np=n1+2;uu=u8;wq=w8
         if(kind1<0)then
            ct=ch8f(:,kr);st=sh8f(:,kr)
         elseif(kind1>0)then
            ct=c8f(:,kr);st=s8f(:,kr)
         endif
      else
         np=5;uu(1:5)=u3;wq(1:5)=w3
         if(kind1<0)then
            ct(1:5)=ch3f(:,kr);st(1:5)=sh3f(:,kr)
         elseif(kind1>0)then
            ct(1:5)=c3f(:,kr);st(1:5)=s3f(:,kr)
         endif
      endif
      select case(kind1)
      case(0)
         do j=1,np
            x=c+el*uu(j)
            ser=cser(1);ser1=0
            do n=2,9
               ser1=ser1*x+ser;ser=ser*x+cser(n)
            enddo
            yv=y0+b*x*ser;y1=b*(ser+x*ser1)
            ff=pser2(1);ff1=0
            do n=2,10
               ff1=ff1*x+ff;ff=ff*x+pser2(n)
            enddo
            q=sqrt(yv);g=yv*ff+a
            f(j)=1/(q*g-t)
            f1(j)=-(.5_dp*y1/q*g+q*(y1*ff+yv*ff1))*f(j)**2
         enddo
      case(1,2)
         ! y = y0+2b*sin(v/2)**2 or yend-2b*cos(v/2)**2; v=pi-x swaps them.
         us=(kind1==1).eqv.(b>=0)
         y00=merge(y0,yend,b>=0);ab=2*abs(b)
         pk=1;if(kind1==2)pk=-1
         mt=merge(m*pi,(m+1)*pi,kind1==1)
         do j=1,np
            x=c+el*uu(j)
            sh=hs*ct(j)+hc*st(j);ch=hc*ct(j)-hs*st(j)
            sn=2*sh*ch;cn=(ch-sh)*(ch+sh);rs=1/sn
            yv=y00+ab*merge(sh,ch,us)**2
            y1=pk*b*sn
            ff=(pk*(x-sn*cn)+mt)*rsqrt2*rs**3
            ff1=(pk*sqrt2-3*cn*ff)*rs
            q=sqrt(yv);g=yv*ff+a
            f(j)=1/(q*g-t)
            f1(j)=-(.5_dp*y1/q*g+q*(y1*ff+yv*ff1))*f(j)**2
         enddo
      case(-2)
         ! xi = wb-w, the distance from the short-way branch point (yend
         ! carries wb): y = 2b*sinh(wb-xi/2)*sinh(xi/2) keeps full relative
         ! precision as y -> 0; sinh(w), cosh(w) by angle addition.
         call shch(yend,sb,cb)
         do j=1,np
            x=c+el*uu(j)
            sh=hs*ct(j)+hc*st(j);ch=hc*ct(j)+hs*st(j)
            sn=sb*(1+2*sh*sh)-cb*2*sh*ch;cn=cb*(1+2*sh*sh)-sb*2*sh*ch   ! sinh(w), cosh(w)
            yv=2*b*(sb*ch-cb*sh)*sh;y1=b*sn
            call fhyp_d(yend-x,sn,cn,ff,ff1)
            q=sqrt(yv);g=yv*ff+a
            f(j)=1/(q*g-t)
            f1(j)=-(.5_dp*y1/q*g+q*(y1*ff-yv*ff1))*f(j)**2
         enddo
      case default
         lw=b<0.and.c>=0.3_dp
         do j=1,np
            x=c+el*uu(j)
            sh=hs*ct(j)+hc*st(j);ch=hc*ct(j)+hs*st(j)
            sn=2*sh*ch;cn=1+2*sh**2;rs=1/sn
            yv=y0-2*b*sh**2;y1=-b*sn
            q=sqrt(yv)
            if(lw)then
               ! Cancellation-free long-way hyperbolic form and its derivative.
               g=((y0+b)*cn-b)*rs**2-yv*x*rs**3
               g1=(y0+b)*rs-2*((y0+b)*cn-b)*cn*rs**3-(y1*x+yv)*rs**3+3*yv*x*cn*rs**4
               f(j)=1/(q*g*rsqrt2-t)
               f1(j)=-(.5_dp*y1/q*g+q*g1)*rsqrt2*f(j)**2
            else
               ff=(sn*cn-x)*rsqrt2*rs**3
               ff1=(sqrt2-3*cn*ff)*rs
               g=yv*ff+a
               f(j)=1/(q*g-t)
               f1(j)=-(.5_dp*y1/q*g+q*(y1*ff+yv*ff1))*f(j)**2
            endif
         enddo
      end select
      ! Q[g'] = g(c+L)-g(c-L) - L*sum w g'(x) for g = f, (x-c)f; dd=Q1/Q0.
      ! Two interleaved partial sums shorten the dependent chain.
      fb=f(np-1);fa=f(np)
      q0=0;q1=0;r0=0;r1=0
      do j=1,np-3,2
         x=el*uu(j)
         q0=q0+wq(j)*f1(j)
         q1=q1+wq(j)*(x*f1(j)+f(j))
         x=el*uu(j+1)
         r0=r0+wq(j+1)*f1(j+1)
         r1=r1+wq(j+1)*(x*f1(j+1)+f(j+1))
      enddo
      if(mod(np,2)==1)then
         j=np-2;x=el*uu(j)
         q0=q0+wq(j)*f1(j)
         q1=q1+wq(j)*(x*f1(j)+f(j))
      endif
      q0=fb-fa-el*(q0+r0)
      q1=el*(fb+fa)-el*(q1+r1)
      dd=q1/q0
      if(.not.(dd==dd))dd=huge(dd)
   end subroutine



   ! Half-angle functions hs,hc (sin/cos, or sinh/cosh for kind1<0) carried
   ! to the new centre x=c+2h by angle addition; the Taylor series of the
   ! shift is exact to rounding for |h|<1/4. A result that cancels is
   ! recomputed at x.
   pure subroutine rotate(kind1,x,h,hs,hc)
      integer,intent(in)::kind1
      real(dp),intent(in)::x,h
      real(dp),intent(inout)::hs,hc
      real(dp)::s,c,q,r
      if(kind1==0)return
      if(abs(h)<.25_dp)then
         q=h*h;if(kind1>0)q=-q
         s=h*(1+q*(rf(3)+q*(rf(5)+q*(rf(7)+q*(rf(9)+q*(rf(11)+q*rf(13)))))))
         c=1+q*(rf(2)+q*(rf(4)+q*(rf(6)+q*(rf(8)+q*(rf(10)+q*(rf(12)+q*rf(14)))))))
      elseif(kind1>0)then
         s=sin(h);c=cos(h)
      else
         call shch(h,s,c)
      endif
      q=hs;r=hc
      hs=q*c+hc*s
      if(kind1>0)then
         hc=hc*c-q*s
      else
         hc=hc*c+q*s
      endif
      if(abs(hs)<abs(q)/4.or.abs(hc)<abs(r)/4)then
         if(kind1>0)then
            hs=sin(x/2);hc=cos(x/2)
         else
            call shch(x/2,hs,hc)
         endif
      endif
   end subroutine

   !> F and its first three derivatives in x: k=+1, x=v (elliptic);
   !! k=-1, x=w (hyperbolic, sinh/cosh).
   pure subroutine f_derivs(x,k,m,f,sx,cx)
   real(dp),intent(in)  :: x
   integer,intent(in)   :: k,m
   real(dp),intent(out) :: f(0:3)
   real(dp),optional,intent(in) :: sx,cx   ! sin/cos (sinh/cosh) of x if known
   real(dp) :: s,c,r,g(0:3),xx
   integer  :: j
   if(k==2)then
      if(present(sx))then
         s=sx;c=cx
      else
         s=sin(x);c=cos(x)
      endif
      r=1/s
      f(0)=((m+1)*pi-x+s*c)/sqrt2*r**3
      f(1)=(-sqrt2-3*c*f(0))*r
      f(2)=(3*s*f(0)-4*c*f(1))*r
      f(3)=(3*c*f(0)+7*s*f(1)-5*c*f(2))*r
   elseif(m==0.and.x<0.3_dp)then
      ! Parabolic series G(X), X=k*x**2, with Taylor coefficients g(j)=G^(j)/j!.
      xx=k*x*x;g=0;g(0)=pser(1)
      do j=2,10
         g(3)=g(3)*xx+g(2);g(2)=g(2)*xx+g(1);g(1)=g(1)*xx+g(0);g(0)=g(0)*xx+pser(j)
      enddo
      g=g/sqrt2
      f(0)=g(0);f(1)=2*k*x*g(1);f(2)=8*x*x*g(2)+2*k*g(1);f(3)=48*k*x**3*g(3)+24*x*g(2)
   else
      ! Closed form and the ODE  s*F' + 3*c*F = sqrt2  (s'=c, c'=-k*s).
      if(present(sx))then
         s=sx;c=cx
      elseif(k>0)then
         s=sin(x);c=cos(x)
      else
         call shch(x,s,c)
      endif
      r=1/s
      f(0)=k*(x+m*pi-s*c)/sqrt2*r**3
      f(1)=(sqrt2-3*c*f(0))*r
      f(2)=(3*k*s*f(0)-4*c*f(1))*r
      f(3)=(3*k*c*f(0)+7*k*s*f(1)-5*c*f(2))*r
   endif
   end subroutine f_derivs

   !> T and its first three derivatives in the anomaly coordinate x (k as above).
   pure subroutine time_derivs(x,k,m,y0,yend,a,b,d,hsx,hcx)
   real(dp),intent(in)  :: x,y0,yend,a,b
   integer,intent(in)   :: k,m
   real(dp),intent(out) :: d(0:3)
   real(dp),optional,intent(in) :: hsx,hcx   ! sin/cos(x/2) if known (k>0)
   real(dp) :: sh,ch,s,c,y,y1,y2,y3,f(0:3),q,q1,q2,q3,p,p1,p2,p3
   if(k>0)then
      if(present(hsx))then
         sh=hsx;ch=hcx
      else
         sh=sin(x/2);ch=cos(x/2)
      endif
      c=(ch-sh)*(ch+sh)
   endif
   if(k==2)then
      y=y0+2*b*ch*ch
      if(b<0)y=yend-2*b*sh*sh
   elseif(k>0)then
      if(b>=0)then
         y=y0+2*b*sh*sh
      else
         y=yend-2*b*ch*ch
      endif
   else
      call shch(x/2,sh,ch);c=1+2*sh*sh
      y=y0-2*b*sh*sh
   endif
   s=2*sh*ch
   y1=k*b*s;y2=k*b*c;y3=-b*s
   if(k==2)then
      y1=-b*s;y2=-b*c;y3=b*s
   endif
   call f_derivs(x,k,m,f,s,c)
   q=sqrt(y);q1=y1/(2*q);q2=(y2/2-q1*q1)/q;q3=(y3/2-3*q1*q2)/q
   p=y*f(0)+a;p1=y1*f(0)+y*f(1);p2=y2*f(0)+2*y1*f(1)+y*f(2)
   p3=y3*f(0)+3*y2*f(1)+3*y1*f(2)+y*f(3)
   d(0)=q*p
   d(1)=q1*p+q*p1
   if(k<0.and.b<0.and.x>=0.3_dp)then
      p=((y0+b)*c-b)/s**2-y*x/s**3
      p1=(y0+b)/s-2*((y0+b)*c-b)*c/s**3-(y1*x+y)/s**3+3*y*x*c/s**4
      d(0)=q*p/sqrt2;d(1)=(q1*p+q*p1)/sqrt2
   endif
   d(2)=q2*p+2*q1*p1+q*p2
   d(3)=q3*p+3*q2*p1+3*q1*p2+q*p3
   end subroutine time_derivs

   ! sinh and cosh together: Taylor series below 1/2 (exact to rounding),
   ! otherwise one exp. A separate sinh/cosh pair costs about twice as much
   ! in some libms (glibc).
   elemental subroutine shch(x,s,c)
      real(dp),intent(in)::x
      real(dp),intent(out)::s,c
      real(dp)::q,e
      if(abs(x)<.5_dp)then
         q=x*x
         s=x*(1+q*(rf(3)+q*(rf(5)+q*(rf(7)+q*(rf(9)+q*(rf(11)+q*rf(13)))))))
         c=1+q*(rf(2)+q*(rf(4)+q*(rf(6)+q*(rf(8)+q*(rf(10)+q*(rf(12)+q*rf(14)))))))
      elseif(abs(x)<700)then
         e=exp(x);q=1/e
         s=(e-q)/2;c=(e+q)/2
      else
         s=sinh(x);c=cosh(x)
      endif
   end subroutine

   ! sin(c/2) and cos(c/2) for 0<=c<=pi from polynomials: reduction
   ! to |y|<=pi/4 (pi/2-c/2 is exact there, Sterbenz) and Taylor polynomials
   ! (truncation below 5e-17). Other c use the intrinsics.
   elemental subroutine sincos_half(c,s,co)
      real(dp),intent(in)::c
      real(dp),intent(out)::s,co
      real(dp),parameter::pio2_hi=1.5707963267948966_dp,pio2_lo=6.123233995736766e-17_dp
      real(dp)::y,z,ps,pc
      if(.not.(c>=0.and.c<=pi))then
         s=sin(c/2);co=cos(c/2);return
      endif
      y=c/2
      if(y>pio2_hi/2)y=(pio2_hi-y)+pio2_lo
      z=y*y
      ps=y+y*z*(-rf(3)+z*(rf(5)+z*(-rf(7)+z*(rf(9)+z*(-rf(11)+z*(rf(13)-z*rf(15)))))))
      pc=1+z*(-rf(2)+z*(rf(4)+z*(-rf(6)+z*(rf(8)+z*(-rf(10)+z*(rf(12)+z*(-rf(14)+z*rf(16))))))))
      if(c/2>pio2_hi/2)then
         s=pc;co=ps
      else
         s=ps;co=pc
      endif
   end subroutine

   ! 2**(-k) (|k|<1000) and 1-exponent(x) (x>0) from the bit pattern;
   ! gfortran implements the scale/exponent intrinsics as library calls.
   elemental real(dp) function pow2m(k) result(y)
      integer,intent(in)::k
      if(ieee64)then
         y=transfer(shiftl(int(1023-k,int64),52),y)
      else
         y=scale(1._dp,-k)
      endif
   end function

   elemental integer function kexp(x) result(k)
      real(dp),intent(in)::x
      if(ieee64.and.x>=tiny(x))then
         k=1023-int(shiftr(transfer(x,0_int64),52))
      elseif(x>0)then
         k=1-exponent(x)
      else
         k=huge(k)
      endif
   end function

   ! Natural logarithm to ~1e-9 relative: exponent bits, mantissa reduced to
   ! [sqrt(1/2),sqrt(2)], and the atanh series in z=(m-1)/(m+1), |z|<0.18.
   elemental real(dp) function rlog(x) result(y)
      real(dp),intent(in)::x
      real(dp),parameter::ln2=log(2._dp)
      integer(int64)::i,e
      real(dp)::m,z,q
      if(.not.(ieee64.and.x>=tiny(x).and.x<=huge(x)))then
         y=log(x);return
      endif
      i=transfer(x,i);e=shiftr(i,52)-1023
      m=transfer(ior(iand(i,4503599627370495_int64),4607182418800017408_int64),m)
      if(m>sqrt2)then
         m=m/2;e=e+1
      endif
      z=(m-1)/(m+1);q=z*z
      y=e*ln2+2*z*(1+q*(1._dp/3+q*(.2_dp+q*(1._dp/7+q*(1._dp/9)))))
   end function

   ! Cube root to ~2e-5 for starters and table coordinates (the contour sets
   ! the accuracy): exponent/3 bit estimate (3%) and one Halley correction; a
   ! second correction (fine) gives ~1e-15.
   elemental real(dp) function rcbrt(x,fine) result(y)
      real(dp),intent(in)::x
      logical,optional,intent(in)::fine
      real(dp)::w,r
      if(.not.(ieee64.and.x>=tiny(x).and.x<=huge(x)))then
         y=cbrt_c(x);return
      endif
      y=transfer(transfer(x,0_int64)/3+3071306045661749710_int64,y)
      w=1/y;r=(x*w)*(w*w)
      y=y*(1+2*r)/(2+r)
      if(present(fine))then
         if(fine)then
            w=1/y;r=(x*w)*(w*w)
            y=y*(1+2*r)/(2+r)
         endif
      endif
   end function

   ! Zero-revolution elliptic root v from the tables at tau=lam/(1+sqrt(1-lam**2))
   ! and phi=(tp/t)**(1/3), by local 4x4 cubic Lagrange interpolation:
   ! A=(pi-v)/(phi*sc) for phi<=1/2, sc=(K/tp)**(1/3) from the pole asymptote
   ! T ~ K/(pi-v)**3, and B=v**2/sc above, sc=(t-tp)/T0', T0'=dT/d(v**2) at
   ! the parabola. Both tend to 1 at their ends. 0 outside the tables.
   pure real(dp) function table_v(tau,phi,sc) result(v)
      real(dp),intent(in)::tau,phi,sc
      real(dp)::u,wt(4),ws(4),r(4)
      integer::i,j,h
      v=0
      if(.not.(abs(tau)<=tab_tmax.and.phi>0.and.phi<=.99_dp))return
      u=(tau+tab_tmax)*((tab_n-1)/(2*tab_tmax))
      i=min(max(int(u)-1,0),tab_n-4);call lagrange4(u-i,wt)
      h=1;u=phi*(2*(tab_n-1))
      if(phi>.5_dp)then
         h=2;u=u-(tab_n-1)
      endif
      j=min(max(int(u)-1,0),tab_n-4);call lagrange4(u-j,ws)
      ! Four independent column sums, then one combination (short chains).
      r=wt(1)*tab(i+1,j+1:j+4,h)+wt(2)*tab(i+2,j+1:j+4,h)+wt(3)*tab(i+3,j+1:j+4,h)+wt(4)*tab(i+4,j+1:j+4,h)
      u=(ws(1)*r(1)+ws(2)*r(2))+(ws(3)*r(3)+ws(4)*r(4))
      if(h==1)then
         v=pi-phi*sc*u
      else
         v=sqrt(max(u*sc,0._dp))
      endif
      if(.not.(v>0.and.v<pi))v=0
   end function

   ! Short-way elliptic root v from tabp at tau in [0.75,1-1e-6] (grid uniform in
   ! 1-sqrt(1-tau)) and s in [0,1]: side 1 (t>=T(pi/2)) v=pi-s*sc*A, side 2
   ! v=sqrt(B*sc) (sc from zero_rev_start). 0 outside the table.
   pure real(dp) function table_vp(tau,s,sc,side) result(v)
      real(dp),intent(in)::tau,s,sc
      integer,intent(in)::side
      real(dp)::u,wt(4),ws(4),r(4)
      integer::i,j
      v=0
      u=min(1-sqrt(1-tau),tabp_z1)              ! clamped beyond tau=1-1e-6
      if(.not.(s>=0.and.s<=1))return
      u=(u-tabp_z0)*((tab_n-1)/(tabp_z1-tabp_z0))
      i=min(max(int(u)-1,0),tab_n-4);call lagrange4(u-i,wt)
      u=s*(tab_n-1)
      j=min(max(int(u)-1,0),tab_n-4);call lagrange4(u-j,ws)
      r=wt(1)*tabp(i+1,j+1:j+4,side)+wt(2)*tabp(i+2,j+1:j+4,side)+wt(3)*tabp(i+3,j+1:j+4,side)+wt(4)*tabp(i+4,j+1:j+4,side)
      u=(ws(1)*r(1)+ws(2)*r(2))+(ws(3)*r(3)+ws(4)*r(4))
      if(side==1)then
         v=pi-s*sc*u
      else
         v=sqrt(max(u*sc,0._dp))
      endif
      if(.not.(v>0.and.v<pi))v=0
   end function

   ! Long-way elliptic root as d=pi-v toward tau -> -1, th=t/|b|**1.5,
   ! eta**2=2*yend/|b|, within 6.4% of min(v,d). With dl=th-pi/4 (the plateau
   ! pi/4): the pole layer (pi/4)((1+eta**2/d**2)**1.5-1)=dl when that gives
   ! d<2*eta; else the transition d**3/2+dl*d**2=3*pi*eta**2/8 (linear body,
   ! linearized layer) when that gives d<0.1; else the body limit from tabb.
   pure real(dp) function long_d(th,eta) result(d)
      real(dp),intent(in)::th,eta
      real(dp)::dl,x,e23,chi,u,wt(4),h
      integer::j
      dl=th-pi/4
      if(dl>0)then
         x=4*dl/pi
         if(x<1e-3_dp)then
            x=x*(2._dp/3-x*(1._dp/9-x*(4._dp/81)))     ! (1+x)**(2/3)-1
         else
            x=exp(log1p_c(x)*(2._dp/3))-1
         endif
         d=eta/sqrt(x)
         if(d<2*eta)return
      endif
      e23=cbrt_c(eta*eta);chi=dl/e23
      if(chi>0)then
         ! 1/x is the largest root of w**3-(8chi/(3pi))w-4/(3pi)=0
         x=1/depressed_max(-8*chi/(3*pi),-4/(3*pi))
      elseif(chi<0)then
         x=depressed_max(-4*chi*chi/3,16*chi**3/27-3*pi/4)-2*chi/3
      else
         x=cbrt_c(3*pi/4)
      endif
      d=e23*x
      if(d<.1_dp)return
      u=min(max((th-1._dp/3)/(pi/4-1._dp/3),0._dp),1._dp)
      h=merge(u/.5_dp,(u-.5_dp)/.5_dp,u<=.5_dp)*(tabb_n-1)
      j=min(max(int(h)-1,0),tabb_n-4);call lagrange4(h-j,wt)
      if(u<=.5_dp)then
         h=wt(1)*tabb(j+1,1)+wt(2)*tabb(j+2,1)+wt(3)*tabb(j+3,1)+wt(4)*tabb(j+4,1)
         d=pi-sqrt(max(h*u,0._dp))
      else
         h=wt(1)*tabb(j+1,2)+wt(2)*tabb(j+2,2)+wt(3)*tabb(j+3,2)+wt(4)*tabb(j+4,2)
         d=h*(1-u)
      endif
   end function

   ! Largest real root of y**3+p*y+q=0 for p<0, in the trigonometric or
   ! hyperbolic form (cancellation-free).
   pure real(dp) function depressed_max(p,q) result(y)
      real(dp),intent(in)::p,q
      real(dp)::m,c
      m=sqrt(-p/3);c=-q/(2*m**3)
      if(abs(c)<=1)then
         y=2*m*cos(acos(c)/3)
      else
         y=sign(2*m*cosh(log(abs(c)+sqrt(c*c-1))/3),c)
      endif
   end function

   ! Hyperbolic table starter for u=t/tp in [0.1,0.99] on a grid in sqrt|tau|
   ! (clamped at 1-1e-6, where the limit w/wb -> sqrt(1-u**2) is regular):
   ! w**2 = tabh*(1-u)*wn**2 with L=1-log(u), wn=L on the long way and
   ! wn**2 = L**2*(1-exp(-wb**2/L**2)) on the short way (wb the branch point).
   pure real(dp) function table_w(tau,u,wb) result(w)
      real(dp),intent(in)::tau,u,wb
      real(dp)::r,l2,z,s,c,g,wt(4),ws(4),cl(4)
      integer::i,j,h
      w=0
      if(.not.(u>=tabh_u0.and.u<=.99_dp))return
      h=merge(1,2,tau>0)
      g=min(sqrt(abs(tau))*((tab_n-1)/tabh_s1),tab_n-1._dp)
      i=min(max(int(g)-1,0),tab_n-4);call lagrange4(g-i,wt)
      g=(u-tabh_u0)*((tab_n-1)/(1-tabh_u0))
      j=min(max(int(g)-1,0),tab_n-4);call lagrange4(g-j,ws)
      cl=wt(1)*tabh(i+1,j+1:j+4,h)+wt(2)*tabh(i+2,j+1:j+4,h)+wt(3)*tabh(i+3,j+1:j+4,h)+wt(4)*tabh(i+4,j+1:j+4,h)
      r=(ws(1)*cl(1)+ws(2)*cl(2))+(ws(3)*cl(3)+ws(4)*cl(4))
      l2=(1-rlog(u))**2
      if(wb>0)then
         z=wb*wb/l2
         if(z<1)then
            call shch(z/2,s,c);l2=l2*2*s*(c-s)          ! 1-exp(-z) without cancellation
         else
            l2=l2*(1-exp(-z))
         endif
      endif
      w=sqrt(max(r*(1-u)*l2,0._dp))
      if(wb>0.and..not.(w<wb))w=0
      if(.not.(w>0))w=0
   end function

   ! Near the parabola, X = k*x**2 (v**2 for an ellipse, -w**2 for a
   ! hyperbola) from the series T = T0 + T1*X + T2*X**2 + ..., with
   ! y = y0 + b*X*(1/2 - X/24 + ...) and y*F + a from the parabolic series of
   ! F, reverted to second order: X = d/T1 - T2*d**2/T1**3, d = t - T0.
   pure subroutine near_parabola_x(y0,a,b,t,tp,xx,rc)
      real(dp),intent(in)::y0,a,b,t,tp
      real(dp),intent(out)::xx,rc       ! X, and the relative size of its second-order term
      real(dp)::q0,q1,q2,g0,g1,g2,p0,p1,p2,t1,t2,d
      g0=pser2(10);g1=pser2(9);g2=pser2(8)
      q0=sqrt(y0);q1=b/(4*q0);q2=(-b/24-q1*q1)/(2*q0)
      p0=y0*g0+a;p1=y0*g1+b/2*g0;p2=y0*g2+b/2*g1-b/24*g0
      t1=q0*p1+q1*p0;t2=q0*p2+q1*p1+q2*p0
      d=t-tp
      rc=t2*d/(t1*t1)
      xx=d/t1*(1-rc)
   end subroutine

   ! Multi-revolution table starter: tau and the semiperimeter sl, time t,
   ! revolutions m and branch br. stat=5: t is below the tabulated minimum
   ! time by more than the table error (infeasible); stat=1: v is the branch
   ! root, within 2.5% of its distance to v* (or to the pole), from the
   ! normalized table v = v* + sigma*(P-v*)*q; stat=2: near the fold
   ! (|sigma| < 1e-2, or within the table's t_min margin), vs is v* within
   ! 6e-5 of min(v*,pi-v*) for mrev_fold; stat=0: outside the tables. Slices
   ! in mu=1/m by 4-point Lagrange.
   pure subroutine mrev_table(tau,sl,t,m,br,v,vs,stat)
      real(dp),intent(in)::tau,sl,t
      integer,intent(in)::m,br
      real(dp),intent(out)::v,vs
      integer,intent(out)::stat
      real(dp),parameter::margin=1e-6_dp
      real(dp)::u,wt(4),wc(4),ws(4),wm(4),mu,cm,r,sig
      integer::i,ic,j,k,km(4),nk,ia,jb,kk,h
      stat=0;v=0;vs=0
      if(.not.(abs(tau)<=tab_tmax))return
      mu=1._dp/m;nk=0
      do k=1,mr_nm
         if(mr_m(k)==m)then
            nk=1;km(1)=k;wm(1)=1
         endif
      enddo
      if(nk==0)then
         kk=mr_nm-1
         do k=1,mr_nm-1
            if(mu>=mr_mu(k+1))then
               kk=k;exit
            endif
         enddo
         kk=min(max(kk-1,1),mr_nm-3);nk=4
         do ia=1,4
            km(ia)=kk+ia-1;wm(ia)=1
            do jb=1,4
               if(jb/=ia)wm(ia)=wm(ia)*(mu-mr_mu(kk+jb-1))/(mr_mu(kk+ia-1)-mr_mu(kk+jb-1))
            enddo
         enddo
      endif
      u=(tau+tab_tmax)*((mr_nc-1)/(2*tab_tmax))
      ic=min(max(int(u)-1,0),mr_nc-4);call lagrange4(u-ic,wc)
      cm=0
      do k=1,nk
         cm=cm+wm(k)*(wc(1)*mr_c(ic+1,km(k))+wc(2)*mr_c(ic+2,km(k))+wc(3)*mr_c(ic+3,km(k))+wc(4)*mr_c(ic+4,km(k)))
         vs=vs+wm(k)*(wc(1)*mr_s(ic+1,km(k))+wc(2)*mr_s(ic+2,km(k))+wc(3)*mr_s(ic+3,km(k))+wc(4)*mr_s(ic+4,km(k)))
      enddo
      r=cm*sqrt(sl)*rcbrt(m/t,.true.)         ! (t_min/t)**(1/3)
      if(r>1+margin)then
         stat=5;return
      endif
      stat=2
      if(r>=1-margin)return
      sig=sqrt(1-r)
      if(sig<1e-2_dp)return
      stat=0;h=merge(2,1,br>0)
      u=(tau+tab_tmax)*((mr_nt-1)/(2*tab_tmax))
      i=min(max(int(u)-1,0),mr_nt-4);call lagrange4(u-i,wt)
      u=sig*(mr_ns-1)
      j=min(max(int(u)-1,0),mr_ns-4);call lagrange4(u-j,ws)
      do k=1,nk
         do jb=1,4
            v=v+wm(k)*ws(jb)*(wt(1)*mr_q(i+1,j+jb,h,km(k))+wt(2)*mr_q(i+2,j+jb,h,km(k)) &
               +wt(3)*mr_q(i+3,j+jb,h,km(k))+wt(4)*mr_q(i+4,j+jb,h,km(k)))
         enddo
      enddo
      v=vs+sig*(merge(pi,0._dp,br>0)-vs)*br*v
      if(v>0.and.v<pi.and.vs>0.and.vs<pi.and.(v-vs)*br>0)stat=1
   end subroutine

   ! Hyperbolic F(w) = (sinh(w)*cosh(w)-w)/(sqrt2*sinh(w)**3) and dF/dw from
   ! sinh(w) and cosh(w); the parabolic series below w = 0.3, where the closed
   ! form cancels.
   pure subroutine fhyp_d(w,sh,ch,f,f1)
      real(dp),intent(in)::w,sh,ch
      real(dp),intent(out)::f,f1
      real(dp)::xx
      integer::j
      if(w<.3_dp)then
         xx=-w*w;f=pser2(1);f1=0
         do j=2,10
            f1=f1*xx+f;f=f*xx+pser2(j)
         enddo
         f1=-2*w*f1
      else
         f=(sh*ch-w)*rsqrt2/sh**3;f1=(sqrt2-3*ch*f)/sh
      endif
   end subroutine

   ! Lagrange velocities v1 = (r2-f*r1)/g, v2 = (gd*r2-r1)/g with f = 1-y/|r1|,
   ! gd = 1-y/|r2|, g = A*sqrt(y), A = d*sqrt(|r1||r2|+r1.r2) (y normalized by
   ! scale), from the exact inputs with every scalar carried as a head and a
   ! tail (exact products by FMA, compensated sums), so each component has
   ! one final rounding.
   pure subroutine fg_velocities(r1,r2,yn,scale,d,v1,v2)
      real(dp),intent(in)::r1(3),r2(3),yn,scale
      integer,intent(in)::d
      real(dp),intent(out)::v1(3),v2(3)
      real(dp)::yh,yl,s1h,s1l,s2h,s2l,dh,dl,n1h,n1l,n2h,n2l,ah,al,qh,ql,gh,gl,ih,il, &
                fh,fl,uh,ul,p,e,s,es,t
      integer::k
      yh=yn*scale;yl=ieee_fma(yn,scale,-yh)
      call dot2(r1,r1,s1h,s1l);call dot2(r2,r2,s2h,s2l);call dot2(r1,r2,dh,dl)
      call sqrt2s(s1h,s1l,n1h,n1l);call sqrt2s(s2h,s2l,n2h,n2l)
      ! A**2 = |r1||r2| + r1.r2
      p=n1h*n2h;e=ieee_fma(n1h,n2h,-p)+(n1h*n2l+n1l*n2h)
      call fsum(p,dh,s,es);es=es+e+dl
      call sqrt2s(s,es,ah,al);ah=d*ah;al=d*al
      ! g = A*sqrt(y), then 1/g
      call sqrt2s(yh,yl,qh,ql)
      gh=ah*qh;gl=ieee_fma(ah,qh,-gh)+(ah*ql+al*qh)
      ih=1/gh;il=(ieee_fma(-gh,ih,1._dp)-gl*ih)*ih
      ! f = 1-y/|r1|, gd = 1-y/|r2|
      call quot(yh,yl,n1h,n1l,qh,ql);call fsum(1._dp,-qh,fh,fl);fl=fl-ql
      call quot(yh,yl,n2h,n2l,qh,ql);call fsum(1._dp,-qh,uh,ul);ul=ul-ql
      do k=1,3
         ! v1 = (r2 - f*r1)*(1/g)
         p=fh*r1(k);e=ieee_fma(fh,r1(k),-p)+fl*r1(k)
         call fsum(r2(k),-p,s,es);es=es-e
         t=s*ih;v1(k)=t+(ieee_fma(s,ih,-t)+(s*il+es*ih))
         ! v2 = (gd*r2 - r1)*(1/g)
         p=uh*r2(k);e=ieee_fma(uh,r2(k),-p)+ul*r2(k)
         call fsum(p,-r1(k),s,es);es=es+e
         t=s*ih;v2(k)=t+(ieee_fma(s,ih,-t)+(s*il+es*ih))
      enddo
   contains
      pure subroutine fsum(a,b,sm,er)            ! sm+er = a+b exactly
         real(dp),intent(in)::a,b
         real(dp),intent(out)::sm,er
         real(dp)::v
         sm=a+b;v=sm-a;er=(a-(sm-v))+(b-v)
      end subroutine
      pure subroutine dot2(u,v,h,l)            ! u.v as h+l (compensated)
         real(dp),intent(in)::u(3),v(3)
         real(dp),intent(out)::h,l
         real(dp)::pr,q,er,w
         integer::i
         h=u(1)*v(1);l=ieee_fma(u(1),v(1),-h)
         do i=2,3
            pr=u(i)*v(i);er=ieee_fma(u(i),v(i),-pr)
            call fsum(h,pr,w,q);h=w;l=l+(q+er)
         enddo
      end subroutine
      pure subroutine sqrt2s(xh,xl,h,l)        ! sqrt(xh+xl) as h+l
         real(dp),intent(in)::xh,xl
         real(dp),intent(out)::h,l
         h=sqrt(xh);l=(ieee_fma(-h,h,xh)+xl)/(2*h)
      end subroutine
      pure subroutine quot(xh,xl,zh,zl,h,l)    ! (xh+xl)/(zh+zl) as h+l
         real(dp),intent(in)::xh,xl,zh,zl
         real(dp),intent(out)::h,l
         h=xh/zh;l=(ieee_fma(-zh,h,xh)+xl-zl*h)/zh
      end subroutine
   end subroutine

   ! Cubic Lagrange weights on the nodes 0,1,2,3 at p.
   pure subroutine lagrange4(p,w)
      real(dp),intent(in)::p
      real(dp),intent(out)::w(4)
      real(dp)::p0,p1,p2,p3
      p0=p;p1=p-1;p2=p-2;p3=p-3
      w(1)=-(1._dp/6)*(p1*p2)*p3;w(2)=.5_dp*(p0*p2)*p3
      w(3)=-.5_dp*(p0*p1)*p3;w(4)=(1._dp/6)*(p0*p1)*p2
   end subroutine

   pure function cross(u,v) result(w)
   real(dp),intent(in) :: u(3),v(3)
   real(dp) :: w(3)
   w=[u(2)*v(3)-u(3)*v(2),u(3)*v(1)-u(1)*v(3),u(1)*v(2)-u(2)*v(1)]
   end function cross

end module excali_solver
