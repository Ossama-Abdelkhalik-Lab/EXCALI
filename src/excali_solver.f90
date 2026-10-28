! EXCALI (EXact Cauchy-Analytic Lambert Integration): finite rectangular-contour Lambert root extraction across conics.
! Extended from Exact v2.2's thin-rectangle ratio. The root comes from the
! zero-height limit of that rectangle (real Gauss-Legendre quadrature of
! 1/(T-t) and its derivative), applied in one or two fixed stages without a
! Newton polish. Geometry, half-anomaly time functions, derivatives and
! closed-form starters are adapted from ConfluentContour's former contour_simple. No
! confluent/Householder root steps and no calls to another Lambert solver.
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
   integer,parameter  :: max_steps=120
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
   ! The tables are initialized read-only module data, not parameters:
   ! some gfortran versions copy a whole parameter array to the stack when
   ! it is indexed at run time.
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
   ! Reciprocal factorials 1/j! (constant multiplies, not divisions).
   real(dp),parameter :: rf(16)=[(1/gamma(j_+1._dp),j_=1,16)]
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

   ! General path (no tables): closed-form starters, the minimum-time
   ! search for multiple revolutions, zero-height rectangles in one or two
   ! stages and the bracketed contour for rare starter misses. Returns the
   ! root in coordinate k (or Y=sqrt(y) in yy for the shortest hyperbolas),
   ! with its half-angle functions when trig. status 0/1 on success.
   subroutine general_root(g,m,br,k,x,hs,hc,trig,yy,status)
   type(frame),intent(in) :: g
   integer,intent(in)     :: m,br
   integer,intent(out)    :: k,status
   real(dp),intent(out)   :: x,hs,hc,yy
   logical,intent(out)    :: trig
   real(dp) :: rn1,rn2,ss,a,b,y0,yend,ypi,t,tp,tv2,lo,hi,phi,sig,tmin,d(0:3),f(0:3), &
               k0,kpi,c,dc,p,wb,psi,sw,dt
   logical :: merged
   rn1=g%rn1;rn2=g%rn2;ss=g%ss;a=g%a;b=g%b;y0=g%y0;yend=g%yend;ypi=g%ypi;t=g%t
   merged=.false.;trig=.false.;status=-4;x=0;hs=0;hc=1

   !======= Starters and admissible interval =======
   k=1;lo=0;hi=pi;yy=-1
   if(m==0)then
      tp=sqrt(y0)*(y0*sqrt2/3+a)                                 ! parabolic time
      tv2=y0*sqrt(y0)*(rsqrt2/5)+b*(sqrt(y0)*rsqrt2+a/(2*sqrt(y0)))/2 ! dT/d(v**2) at 0
      if(abs(t-tp)<=4*epsilon(t)*tp)then
         x=0                                                     ! parabola
      elseif(t>tp)then
         phi=rcbrt(tp/t)
         ! Ellipse: cubic Hermite for v**2 in phi=(tp/t)**(1/3), between
         ! v**2=(t-tp)/tv2 at phi=1 and pi-v=phi*(K/tp)**(1/3) at phi=0,
         ! with the exact upper pole T ~ pi*y(pi)**1.5/(sqrt2*(pi-v)**3).
         x=hermite(phi,pi*pi,-2*pi*sqrt(ypi)*rcbrt((pi*rsqrt2)/tp),0._dp,-3*tp/tv2)
         if(.not.(x>0))x=(t-tp)/tv2         ! model undershoot: parabolic limit
         x=min(sqrt(x),pi*(1-1e-3_dp))
         if(x>.75_dp*pi)then
            ! Near the pole use its asymptote in d=pi-v with y~y(pi)-b*d**2/2:
            ! T ~ pi*y**1.5/(sqrt2*d**3)+a*sqrt(y), so y/d**2 = Q, solved for d.
            p=t-a*sqrt(ypi)
            if(p>0)then
               p=rcbrt(sqrt2*p/pi)**2+b/2
               if(p>0)x=pi-min(sqrt(ypi/p),pi/4)
            endif
            ! Long way with y_end<<|b|: T plateaus at ~pi*|b|**1.5/4 before
            ! the pole. Below it, cubic Hermite for v**2 in T between the
            ! parabola (slope 1/tv2) and v=pi (dT/dv=|b|**1.5/2).
            if(x>.75_dp*pi.and.b<0)then
               p=pi*(-b)*sqrt(-b)/4
               if(t<p.and.p>tp)then
                  x=hermite((t-tp)/(p-tp),0._dp,(p-tp)/tv2,pi*pi,(p-tp)*4*pi/((-b)*sqrt(-b)))
                  x=min(sqrt(max(x,0._dp)),pi*(1-1e-3_dp))
               endif
            endif
         endif
      else
         ! Hyperbola in w=sqrt(-z)/2: w**2=(tp-t)/tv2 near the parabola.
         k=-1;hi=huge(hi)
         if(b>0)then
            ! Short way: T -> 0 at the branch point y(wb)=0, wb=2*asinh(s).
            sw=sqrt(y0/(2*b));dt=sqrt(1+sw*sw)
            wb=2*log1p_c(sw+sw*sw/(1+dt));hi=wb
            if(t<0.1_dp*tp)then
               ! Solve in Y=sqrt(y): T=Y*(Y**2*F+a) with the exact
               ! w(Y)=2*asinh(sqrt((y0-Y**2)/(2b))). Either term of
               ! T = a*Y + F*Y**3 alone overestimates Y; start from the smaller.
               call f_derivs(wb,-1,0,f)
               yy=min(t/a,rcbrt(t/f(0)),0.9_dp*sqrt(y0))
               ! Zero-height rectangles in Y avoid the branch-point
               ! cancellation of hyperbolic anomaly space.
               lo=0;hi=sqrt(y0)
               if(.not.converge_y())return
               sw=(y0-yy*yy)/(2*b);x=2*asinh(sqrt(sw))
            else
               ! Cubic Hermite for w**2 in (t/tp)**2 between the parabola and
               ! wb, where T**2 ~ a**2*y has slope -a**2*b*sinh(wb).
               x=hermite((t/tp)**2,wb*wb,-wb*tp*tp/(a*a*b*sw*dt),0._dp,-tp/(2*tv2))
               if(.not.(x>0))x=(tp-t)/tv2
               x=min(sqrt(x),wb*(1-1e-6_dp))
            endif
         else
            ! Long way: T ~ sqrt(-b)*(y0+b)*exp(-w/2) as w -> infinity.
            psi=rlog(tp/t)
            x=sqrt(tp*psi/tv2)
            if(psi>1)x=max(x,2*rlog(sqrt(-b)*(y0+b)/t))
         endif
      endif
   else
      ! Multiple revolutions. Every ellipse through both points has a>=s/2,
      ! so T > m*2*pi*(s/2)**1.5 (within 31% of the minimum time) rejects
      ! most infeasible times without evaluating T.
      status=5
      dc=sqrt((rn1-rn2)**2+4*rn1*rn2*ss*ss);c=(rn1+rn2+dc)/2     ! chord, semiperimeter
      if(t<m*pi*c*sqrt(c)/sqrt2*(1-1e-12_dp))return
      if(.not.mrev_start())return
   endif

   ! Store distance from the upper pole directly on near-full-period arcs.
   ! Subtracting a converged v from pi loses the small arc information.
   call pole_coordinate()

   !======= Root: zero-height thin rectangles =======
   if(x>0.and.yy<0.and..not.merged)then
      ! Rare starter misses continue as a bracketed contour.
      if(.not.contour_root())then
         if(.not.converge())return
      endif
   endif
   status=0
   if(merged)status=1

   contains

      ! Multi-revolution starter without the tables: minimum time by Halley,
      ! then the branch starter. False (status 5) when t is below t_min.
      logical function mrev_start() result(ok)
         integer::it
         real(dp)::p,dt,s
         ok=.false.
         ! Minimum time by Halley on T'(v) = 0, started at the large-m minimum
         ! of y**1.5/sin(v)**3, cos v = b/(y(0)+b+sqrt(y(0)*y(pi))), which is
         ! off by ~0.15/m. Halley leaves at most ~12*step**3 in v and ~150*
         ! step**4 relative in tmin, so it stops once tmin decides feasibility
         ! and locates the minimum well inside the distance to the root.
         ! Half-angle functions follow from cos v, then by angle addition.
         c=b/(y0+b+sqrt(y0*ypi))
         x=acos(c);hs=sqrt((1-c)/2);hc=sqrt((1+c)/2)
         do it=1,20
            call time_derivs(x,1,m,y0,yend,a,b,d,hs,hc)
            p=-2*d(1)*d(2)/(2*d(2)**2-d(1)*d(3))
            dt=min(max(x+p,x/2),(x+pi)/2)-x
            ! T and T'' at the step from the Taylor data (error O(step**4)).
            tmin=d(0)+dt*(d(1)+dt*(d(2)/2+dt*d(3)/6))
            x=x+dt;call rotate(1,x,dt/2,hs,hc)
            ! A step below 1e-3 leaves ~1e-8 in v, and tmin is quadratic in it.
            if(abs(dt)<=1e-3_dp)exit
            if(abs(dt)<.3_dp.and.tmin>0)then
               p=t-tmin;s=300*tmin*dt**4
               if(p<-s)exit
               if(p>s.and.12*abs(dt)**3<sqrt(2*p/max(d(2),tiny(p)))/8)exit
            endif
         enddo
         ! T'' at the minimum; after a larger early-exit step near the pole the
         ! linear update is unreliable and the starter needs only its scale.
         p=d(2)+dt*d(3)
         if(p>d(2)/2.or.abs(dt)<=1e-3_dp)d(2)=p
         if(t<tmin*(1-16*epsilon(t)).or..not.(d(2)>0))return
         merged=abs(t-tmin)<=16*epsilon(t)*tmin
         if(merged)then
            status=1
            lo=x;hi=x
         endif
         status=-4
         ! Each branch: cubic Hermite in sig=sqrt(1-(tmin/t)**(1/3)) between the
         ! minimum (v*-v ~ sqrt(6*tmin/T'')*sig) and the pole T ~ K/v**3 at sig=1.
         p=tmin/t;sig=rcbrt(p)
         sig=sqrt(max((1-p)/(1+sig+sig*sig),0._dp))   ! 1-p**(1/3) without cancellation
         k0=y0*sqrt(y0)*m*pi/sqrt2;kpi=ypi*sqrt(ypi)*(m+1)*pi/sqrt2
         ! Near a pole the asymptote with y quadratic in the pole distance d,
         ! T ~ n*pi*y**1.5/(sqrt2*d**3)+a*sqrt(y), y ~ y_p+b*d**2/2 (v=0) or
         ! y_p-b*d**2/2 (v=pi), gives y/d**2 = Q in closed form; it also covers
         ! y_p<<|b|, where T plateaus before the pole; used when the Hermite
         ! starter collapses onto the pole and the model root is well inside.
         if(br<0)then
            hi=x
            x=hermite(sig,x,-sqrt(6*tmin/d(2)),0._dp,-2*rcbrt(k0/tmin))
            p=t-a*sqrt(y0)
            if(p>0.and.x<.02_dp*hi)then
               p=rcbrt(sqrt2*p/(m*pi))**2-b/2
               if(p>0)then
                  if(sqrt(y0/p)<hi/2)x=sqrt(y0/p)
               endif
            endif
         else
            lo=x
            x=pi-hermite(sig,pi-x,-sqrt(6*tmin/d(2)),0._dp,-2*rcbrt(kpi/tmin))
            p=t-a*sqrt(ypi)
            if(p>0.and.pi-x<.02_dp*(pi-lo))then
               p=rcbrt(sqrt2*p/((m+1)*pi))**2+b/2
               if(p>0)then
                  if(sqrt(ypi/p)<(pi-lo)/2)x=pi-sqrt(ypi/p)
               endif
            endif
         endif
         if(t-tmin<1e-2_dp*tmin)then
            ! Close to the minimum, invert t-tmin = T''/2*e**2 + T'''/6*e**3 by
            ! series reversion; the root may lie far inside the bracket margin.
            p=sqrt(2*(t-tmin)/d(2))
            x=merge(hi,lo,br<0)+merge(-p,p,br<0)-d(3)/(6*d(2))*p*p
            x=min(max(x,lo+p/2),hi-p/2)
         else
            x=min(max(x,lo+(hi-lo)*1e-3_dp),hi-(hi-lo)*1e-3_dp)
         endif
         ok=.true.
      end function

      subroutine pole_coordinate()
         real(dp)::p
         if(k==1.and.x>0.75_dp*pi.and..not.merged)then
            p=lo;lo=pi-hi;hi=pi-p;x=pi-x;k=2
         endif
      end subroutine

      ! Root from zero-height thin rectangles; no derivative steps. Stage 1
      ! (8 nodes, half-width L~R/4, R the distance to the nearest singularity
      ! of 1/(T-t)) needs only to contain the starter's root, and then usually
      ! reaches rounding level. The ratio stays exact for a pole outside the
      ! interval, so a consistent estimate outside it re-centres. The
      ! disagreement de of two independent ratios measures the error; when
      ! above rounding, stage 2 (3 nodes, half-width ~16*de, far inside R)
      ! converges like (L/2R)**6. Every accepted interval must show a sign
      ! change of T-t. Singularities: branch points of sqrt(y), the mirror
      ! root of the even m=0 time function, the other multirevolution branch,
      ! and the complex roots of T=t near a pole of T. Near the parabola the
      ! coordinate is X=k*x**2, where T is a series in X. On failure the
      ! state is restored for the safeguarded fallback.
      logical function contour_root() result(ok)
         integer::k0
         real(dp)::lo0,hi0
         k0=k;lo0=lo;hi0=hi
         ok=rect_root()
         if(.not.ok)then
            k=k0;lo=lo0;hi=hi0;trig=.false.
         endif
      end function

      logical function rect_root() result(ok)
         real(dp)::c,rr,el,dd,de,beta,eta,grow,clamped,stp,bim
         integer::kind1,kr,pass
         logical::sc,outside,nox,near
         ok=.false.
         ! Branch points of sqrt(y), computed once: for b>0 at v=+-i*beta
         ! (w=beta); for b<0, with s**2=y0/(-2b)=1+(y0+2b)/(-2b), at
         ! v=pi+-i*eta, w=eta+-i*pi (s>1), or w=+-i*bim (s<=1, hyperbolic).
         ! Only distances are needed, so the inverse functions are approximate.
         ! The short-way hyperbola keeps its exact w bound hi, which the root
         ! may approach arbitrarily closely.
         beta=0;eta=0;bim=-1
         if(b>0)then
            beta=2*rasinh(sqrt(y0/(2*b)))
            if(k<0)beta=hi
         endif
         if(b<0)then
            if(y0+2*b>=0)then
               eta=2*rasinh(sqrt((y0+2*b)/(-2*b)))
            else
               bim=2*rasin(sqrt(y0/(-2*b)))
            endif
         endif
         c=x;kind1=k;nox=.false.;de=0
         ! Stage 1 starts at L~R/4, where 8 nodes usually reach rounding
         ! level. A self-consistent estimate outside the interval re-centres
         ! (going straight to stage 2 if already precise enough); an
         ! inconsistent one (root far away) doubles L, up to ~R.
         grow=4
         outside=.false.
         do pass=1,5
            call canonical(kind1,c,de,nox)
            if(.not.trig.and.kind1>0)then
               call sincos_half(c,hs,hc);trig=.true.
            elseif(.not.trig.and.kind1<0)then
               call shch(c/2,hs,hc);trig=.true.
            endif
            rr=reach(kind1,c,beta,eta,bim)
            el=min(1._dp,rr/grow);clamped=el
            if(kind1==0)el=min(el,.2_dp-abs(c))
            ! Never reach a pole of T: the sign test would see its jump.
            if(kind1==2.or.(kind1==1.and.m>0))el=min(el,.9_dp*c)
            if(kind1==1.or.(kind1==2.and.m>0))el=min(el,.9_dp*(pi-c))
            clamped=el-clamped
            if(.not.(el>0))return
            kr=kexp(el)
            if(kr>kmax)return
            call rect_ratio(kind1,c,hs,hc,kr,.true.,y0,yend,a,b,t,m,dd,de,sc)
            el=pow2m(kr)
            if(abs(dd)<.8_dp*el.and.de<.1_dp*el.and.sc)then
               c=c+dd;call rotate(kind1,c,dd/2,hs,hc);exit
            endif
            if(pass==5)return
            if(abs(dd)<8*el.and.de<.5_dp*abs(dd))then
               near=abs(dd)<1.5_dp*el
               c=c+dd;call rotate(kind1,c,dd/2,hs,hc)
               if(kind1/=0.and..not.(c>lo.and.c<hi))return
               call canonical(kind1,c,de,nox)
               ! Just outside and precise: stage 2 certifies its own sign
               ! change. Farther out de is no error bound; repeat stage 1.
               if(near.and.de<reach(kind1,c,beta,eta,bim)/256)then
                  outside=.true.;exit
               endif
            elseif(clamped<0.and.abs(dd)<huge(dd).and.kind1/=0)then
               ! L is capped by a pole: step towards the estimate instead.
               stp=sign(min(abs(dd),el),dd)
               c=c+stp;call rotate(kind1,c,stp/2,hs,hc)
               if(.not.(c>lo.and.c<hi))return
               call canonical(kind1,c,de,nox)
            elseif(kind1==0)then
               ! Root beyond the series region: widen in the anomaly coordinate.
               kind1=1;if(c<0)kind1=-1
               c=sqrt(abs(c));grow=2;nox=.true.;trig=.false.
            elseif(grow>1.1_dp)then
               grow=grow/2
            else
               return
            endif
         enddo
         call canonical(kind1,c,de,nox)
         ! Rounding-level error needs no second stage. Velocities are smooth in
         ! X near the parabola, so there the floor is absolute.
         rr=max(abs(c),el*1e-3_dp)
         if(kind1==0)rr=max(abs(c),1e-2_dp)
         if(outside.or.de>4*epsilon(c)*rr)then
            rr=reach(kind1,c,beta,eta,bim)
            el=min(rr/16,max(16*de,1e-7_dp*rr))
            if(outside)el=min(rr/16,max(16*de,1e-4_dp*rr))
            kr=kexp(el)
            if(kr>kmax)return
            ! 3 nodes suffice far inside R; otherwise use the 8-node rule.
            if(.not.trig.and.kind1>0)then
               call sincos_half(c,hs,hc);trig=.true.
            elseif(.not.trig.and.kind1<0)then
               call shch(c/2,hs,hc);trig=.true.
            endif
            call rect_ratio(kind1,c,hs,hc,kr,el>rr/1024,y0,yend,a,b,t,m,dd,de,sc)
            el=pow2m(kr)
            if(.not.(abs(dd)<.8_dp*el.and.de<1e-3_dp*el.and.sc))return
            c=c+dd;call rotate(kind1,c,dd/2,hs,hc)
         endif
         if(kind1==0)then
            if(c==0)return
            k=1;if(c<0)k=-1
            x=sqrt(abs(c));trig=.false.
         else
            k=kind1;x=c
         endif
         ok=x>0
      end function

      ! Coordinate for a zero-revolution estimate: X=k*x**2 (series) when
      ! |x|<0.3 unless excluded (nox), the anomaly beyond |X|=0.15. The error
      ! estimate de is carried into the new coordinate.
      subroutine canonical(kind1,c,de,nox)
         integer,intent(inout)::kind1
         real(dp),intent(inout)::c,de
         logical,intent(in)::nox
         real(dp)::q
         ! Elliptic v versus distance from the upper pole, with the bracket.
         if((kind1==1.and.c>.8_dp*pi).or.(kind1==2.and.c>.3_dp*pi))then
            q=lo;lo=pi-hi;hi=pi-q;c=pi-c;kind1=3-kind1
            q=hs;hs=hc;hc=q
         endif
         if(m/=0)return
         if(kind1==0.and.abs(c)>.15_dp)then
            kind1=1;if(c<0)kind1=-1
            c=sqrt(abs(c));de=de/(2*c);trig=.false.
         elseif(abs(kind1)==1.and.c<.3_dp.and..not.nox)then
            c=max(c,0._dp);de=(2*c+de)*de;c=kind1*c*c;kind1=0;trig=.false.
         endif
      end subroutine

      ! Distance from c to the nearest singularity of 1/(T-t) in the local
      ! coordinate (kind 0: X; 1: v; 2: pi-v; -1: w).
      real(dp) function reach(kind1,c,beta,eta,bim) result(rr)
         integer,intent(in)::kind1
         real(dp),intent(in)::c,beta,eta,bim
         rr=huge(rr)
         select case(kind1)
         case(0)
            if(b>0)then
               rr=abs(c+beta*beta)
            elseif(bim>=0)then
               rr=abs(c-bim*bim)
            elseif(b<0)then
               rr=sqrt((c-pi*pi+eta*eta)**2+(2*pi*eta)**2)
            endif
            rr=min(rr,.4_dp)
         ! Near a pole of T (v=pi; also v=0 for m>0) T~K/d**3, so T=t has two
         ! complex roots (poles of 1/(T-t)) at sqrt(3) times the distance d.
         case(1)
            if(b>0)rr=sqrt(c*c+beta*beta)
            if(b<0)rr=sqrt((pi-c)**2+eta*eta)
            rr=min(rr,1.6_dp*(pi-c))
            if(m>0)rr=min(rr,1.6_dp*c)
         case(2)
            if(b>0)rr=sqrt((pi-c)**2+beta*beta)
            if(b<0)rr=sqrt(c*c+eta*eta)
            rr=min(rr,1.6_dp*c)
            if(m>0)rr=min(rr,1.6_dp*(pi-c))
         case default
            if(b>0)then
               rr=beta-c
            elseif(bim>=0)then
               rr=sqrt(c*c+bim*bim)
            elseif(b<0)then
               rr=sqrt((c-eta)**2+pi*pi)
            endif
         end select
         if(kind1/=0)then
            if(m==0.and.kind1/=2)then
               rr=min(rr,2*c)
            elseif(m>0)then
               ! Only the minimum end of the bracket matters: beyond it lies the
               ! other branch's root. The pole of T at the other end is a zero
               ! of 1/(T-t), not a singularity.
               if((kind1==1).eqv.(br<0))then
                  rr=min(rr,hi-c)
               else
                  rr=min(rr,c-lo)
               endif
            endif
         endif
      end function



      ! Very short hyperbolas: the bracketed contour in Y=sqrt(y), where
      ! T(Y)=Y*(Y**2*F(w)+a) is increasing and w(Y)=2*asinh(sqrt((y0-Y**2)/(2b)))
      ! avoids the branch-point cancellation of w. Intervals stay inside
      ! (0,sqrt(y0)); w=0 at Y=sqrt(y0) is the only nearby singularity.
      logical function converge_y() result(ok)
         real(dp)::res,trial,el,dd,de,tt,dtt,ft(n1+2),f1(n1+2),q0,q1,q2,xj
         integer::j,i
         ok=.false.;el=(hi-lo)/4
         do j=1,max_steps
            call ytime(yy,tt,dtt)
            res=tt-t
            if(abs(res)<=16*epsilon(t)*t)then
               ok=.true.;return
            endif
            if(res<0)then
               lo=yy
            else
               hi=yy
            endif
            if(.not.(hi-lo>8*epsilon(yy)*yy))then
               ok=abs(res)<=1e-12_dp*t;return
            endif
            el=min(el,hi-lo,.9_dp*yy,.9_dp*(sqrt(y0)-yy))
            trial=(lo+hi)/2
            if(el>0)then
               ! The zero-height rectangle of rect_ratio, in Y.
               do i=1,n1+2
                  call ytime(yy+el*u8(i),tt,dtt)
                  ft(i)=1/(tt-t);f1(i)=-dtt*ft(i)**2
               enddo
               q0=0;q1=0;q2=0
               do i=1,n1
                  xj=el*u8(i)
                  q0=q0+w8(i)*f1(i);q1=q1+w8(i)*(xj*f1(i)+ft(i));q2=q2+w8(i)*xj*(xj*f1(i)+2*ft(i))
               enddo
               q0=ft(n1+1)-ft(n1+2)-el*q0
               q1=el*(ft(n1+1)+ft(n1+2))-el*q1
               q2=el*el*(ft(n1+1)-ft(n1+2))-el*q2
               dd=q1/q0;de=abs(q2/q0-dd*dd)/el
               if(ft(n1+1)*ft(n1+2)<0)then
                  lo=max(lo,yy-el);hi=min(hi,yy+el)
                  if(de<=4*epsilon(yy)*yy.and.abs(dd)<el)then
                     yy=yy+dd;ok=yy>0;return
                  endif
                  el=max(16*de,1e-7_dp*yy)
               elseif(.not.(abs(dd)<huge(dd)))then
                  el=el/2
               elseif(abs(dd)>=el)then
                  if(res<0)then
                     lo=max(lo,yy+el)
                  else
                     hi=min(hi,yy-el)
                  endif
                  el=2*el
               endif
               if(abs(dd)<huge(dd).and.yy+dd>lo.and.yy+dd<hi)trial=yy+dd
            endif
            if(abs(trial-yy)<=8*epsilon(t)*max(abs(yy),tiny(t)))then
               ok=abs(res)<=1e-12_dp*t;return
            endif
            yy=trial
         enddo
      end function

      ! T(Y) and dT/dY on the short-way hyperbola, with sinh(w/2)=s:
      ! dw/dY = -2Y/(b*sinh(w)).
      subroutine ytime(yv,tt,dtt)
         real(dp),intent(in)::yv
         real(dp),intent(out)::tt,dtt
         real(dp)::s,shw,ff(0:3)
         s=sqrt((y0-yv*yv)/(2*b))
         shw=2*s*sqrt(1+s*s)
         call f_derivs(2*asinh(s),-1,0,ff,shw,1+2*s*s)
         tt=yv*(yv*yv*ff(0)+a)
         dtt=3*yv*yv*ff(0)+a-2*yv**4*ff(1)/(b*shw)
      end subroutine

      ! Bracketed contour for rare starter misses: every estimate is a
      ! zero-height rectangle ratio, with no derivative step. T is monotone
      ! on the branch, so the sign of T-t at the centre and the rectangle's
      ! own sign test shrink [lo,hi]; an estimate outside it is replaced by
      ! the midpoint. A rectangle that leaves the domain (dd not finite)
      ! carries no sign information and only halves L.
      logical function converge() result(ok)
         real(dp)::res,trial,sgn,el,dd,de,hs2,hc2,cr,elr,xl,xh,xe,wide
         integer::j,kr
         logical::sc,series
         ok=.false.
         ! Monotone direction of the selected branch.
         sgn=1
         if(k<0.or.(m>0.and.br<0))sgn=-1
         if(k==2)sgn=-sgn
         if(k<0.and.b<0)then
            hi=max(2*x,1._dp)
            do j=1,100
               call time_derivs(hi,k,m,y0,yend,a,b,d)
               if(d(0)<t)exit
               hi=hi*1.5_dp
               if(hi>300)return
            enddo
         endif
         trig=.false.;el=(hi-lo)/4
         do j=1,max_steps
            call time_derivs(x,k,m,y0,yend,a,b,d)
            res=d(0)-t
            if(abs(res)<=16*epsilon(t)*t)then
               ok=.true.;return
            endif
            ! Rectangles stay inside the bracket as it was before this step,
            ! so they never reach the other branch beyond t_min or a pole.
            wide=min(x-lo,hi-x)
            if(res*sgn<0)then
               lo=x
            else
               hi=x
            endif
            if(.not.(hi-lo>8*epsilon(x)*x))then
               ok=abs(res)<=1e-12_dp*t;return
            endif
            ! Never reach a pole of T (v=pi; also v=0 for m>0) or w<0.
            el=min(el,1._dp,.9_dp*x,.9_dp*wide)
            if(k==1.or.(k==2.and.m>0))el=min(el,.9_dp*(pi-x))
            trial=(lo+hi)/2
            ! Near the parabola the rectangle works in X=k*x**2, where T is a
            ! series; in x the closed form cancels (x - sin x cos x ~ x**3).
            series=m==0.and.k/=2.and.x<.3_dp
            cr=x;elr=el
            if(series)then
               cr=k*x*x;elr=min(2*x*el,.2_dp-abs(cr))
            endif
            kr=kexp(elr)
            if(elr>0.and.kr<=kmax)then
               if(series)then
                  call rect_ratio(0,cr,0._dp,1._dp,kr,.true.,y0,yend,a,b,t,m,dd,de,sc)
               else
                  if(k>0)then
                     hs2=sin(x/2);hc2=cos(x/2)
                  else
                     call shch(x/2,hs2,hc2)
                  endif
                  call rect_ratio(k,x,hs2,hc2,kr,.true.,y0,yend,a,b,t,m,dd,de,sc)
               endif
               elr=pow2m(kr)
               ! The rectangle's ends and its estimate, mapped back to x.
               xl=xof(cr-elr,series,k);xh=xof(cr+elr,series,k)
               if(xl>xh)then
                  xe=xl;xl=xh;xh=xe
               endif
               xe=-1
               if(abs(dd)<huge(dd))xe=xof(cr+dd,series,k)
               if(sc)then
                  lo=max(lo,xl);hi=min(hi,xh)
                  if(de<=4*epsilon(x)*max(abs(cr),merge(1e-2_dp,elr*1e-3_dp,series)).and.abs(dd)<elr.and.xe>0)then
                     x=xe;ok=.true.;return
                  endif
                  el=max(16*de,1e-7_dp*abs(cr))
                  if(series)el=el/(2*x)
               elseif(.not.(abs(dd)<huge(dd)))then
                  el=el/2
               else
                  ! No sign change. When the ratio also puts the root outside,
                  ! the bracket moves past the rectangle; near the root both can
                  ! be rounding noise, and then neither moves it.
                  if(abs(dd)>=elr)then
                     if(res*sgn<0)then
                        lo=max(lo,xh)
                     else
                        hi=min(hi,xl)
                     endif
                  endif
                  el=hi-lo
               endif
               if(xe>lo.and.xe<hi)trial=xe
            endif
            if(abs(trial-x)<=8*epsilon(t)*max(abs(x),1e-100_dp))then
               ok=abs(res)<=1e-12_dp*t;return
            endif
            x=trial
         enddo
      end function

   end subroutine general_root







   ! Batch of normalized solves, same arguments per case as excali_one. The
   ! common case (a starter from the tables or the near-parabolic series, or
   ! a parabola, with the transfer angle not within 1e-4 of 0 or pi) runs in stages over
   ! blocks of cases, so independent cases overlap: geometry and table
   ! starter, one 8-node rectangle, velocities. Results are those of
   ! excali_one, which solves every other case and any rectangle a stage does
   ! not accept.
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

   ! nc <= nb normalized solves, stage by stage so that independent cases
   ! overlap: geometry; table starters (zero revolutions: tables and the
   ! near-parabolic series; multiple revolutions: feasibility, minimum-time
   ! location and branch root); one 8-node rectangle whose width is at most
   ! 1/16 of a lower bound on the distance to every singularity, so that it
   ! is final; the general path for anything it does not certify; velocities.
   subroutine solve_chunk(nc,r1,r2,tof,direction,v1,v2,status,half,nrev,branch,plane_normal)
   integer,intent(in)   :: nc
   real(dp),intent(in)  :: r1(3,nc),r2(3,nc),tof(nc)
   integer,intent(in)   :: direction(nc),nrev(nc),branch(nc)
   real(dp),intent(out) :: v1(3,nc),v2(3,nc)
   integer,intent(out)  :: status(nc),half(nc)
   real(dp),optional,intent(in) :: plane_normal(3)
   type(frame) :: g(nb)
   real(dp) :: x(nb),c(nb),hs(nb),hc(nb),wb(nb),lo(nb),hi(nb),dd(nb),de(nb),yy(nb),el,tp,dc,sl
   integer :: k(nb),kd(nb),kr(nb),j
   logical :: fast(nb),sc(nb),near(nb),trig(nb),eight(nb)
   ! Stage 1: geometry and the table coordinates (multiple revolutions: the
   ! tables themselves).
   do j=1,nc
      fast(j)=.false.;near(j)=.false.;half(j)=0;status(j)=-1
      if(abs(direction(j))/=1.or.nrev(j)<0.or.abs(branch(j))/=1.or..not.(tof(j)>0) &
         .or..not.ieee_is_finite(tof(j)))cycle
      call geometry(r1(:,j),r2(:,j),tof(j),direction(j),g(j),status(j),half(j),plane_normal)
      if(status(j)/=-4)cycle
      lo(j)=0;hi(j)=pi;fast(j)=.true.
      if(nrev(j)==0)then
         call table_coordinates(g(j),k(j),c(j),dd(j),wb(j),near(j))
      else
         ! Every ellipse through both points has a>=s/2, so T > m*2*pi*(s/2)**1.5
         ! (within 31% of the minimum time) rejects most infeasible times.
         dc=sqrt((g(j)%rn1-g(j)%rn2)**2+4*g(j)%rn1*g(j)%rn2*g(j)%ss**2);sl=(g(j)%rn1+g(j)%rn2+dc)/2
         status(j)=5
         if(g(j)%t>=nrev(j)*pi*sl*sqrt(sl)/sqrt2*(1-1e-12_dp)) &
            call mrev_table(g(j)%a/(sqrt2*sl)/(1+sqrt(dc/sl)),sl,g(j)%t,nrev(j),branch(j),x(j),tp,status(j))
         fast(j)=status(j)==1;k(j)=1
         if(status(j)/=5)status(j)=-4
         if(branch(j)<0)then
            hi(j)=tp
         else
            lo(j)=tp
         endif
         if(x(j)>.75_dp*pi)then                    ! distance from the upper pole
            el=lo(j);lo(j)=pi-hi(j);hi(j)=pi-el;x(j)=pi-x(j);k(j)=2
         endif
      endif
   enddo
   ! Stage 2: zero-revolution table starters, then the rectangle's
   ! coordinate and width.
   do j=1,nc
      if(status(j)/=-4.or.nrev(j)>0)cycle
      call table_root(k(j),c(j),dd(j),wb(j),near(j),x(j))
      fast(j)=x(j)>0.or.k(j)==0
   enddo
   do j=1,nc
      if(.not.fast(j).or.k(j)==0)cycle
      call fast_coordinate(k(j),x(j),nrev(j),g(j)%b,wb(j),lo(j),hi(j),kd(j),c(j),el,eight(j),hs(j),hc(j))
      fast(j)=el>0;kr(j)=kexp(el)
   enddo
   ! Stage 3: one 8-node rectangle each, then its certificate.
   do j=1,nc
      if(fast(j).and.k(j)/=0)call rect_ratio(kd(j),c(j),hs(j),hc(j),kr(j),eight(j),g(j)%y0,g(j)%yend,g(j)%a,g(j)%b, &
         g(j)%t,nrev(j),dd(j),de(j),sc(j))
   enddo
   do j=1,nc
      if(fast(j).and.k(j)/=0)call fast_accept(kd(j),c(j),kr(j),dd(j),de(j),sc(j),hs(j),hc(j),k(j),x(j),fast(j))
   enddo
   ! Stage 4: the general path for the rest, then velocities.
   do j=1,nc
      yy(j)=-1;trig(j)=.true.
      if(fast(j))then
         status(j)=0
         if(k(j)==0)then
            k(j)=1;x(j)=0;hs(j)=0;hc(j)=1              ! parabola
         endif
      elseif(status(j)==-4)then
         call general_root(g(j),nrev(j),branch(j),k(j),x(j),hs(j),hc(j),trig(j),yy(j),status(j))
         if(yy(j)<0.and..not.trig(j))then
            if(k(j)>0)then
               call sincos_half(x(j),hs(j),hc(j))
            else
               call shch(x(j)/2,hs(j),hc(j))
            endif
         endif
      endif
      if(status(j)==0.or.status(j)==1)then
         call velocities(g(j),r1(:,j),r2(:,j),direction(j),k(j),x(j),hs(j),hc(j),near(j),yy(j),v1(:,j),v2(:,j))
         if(.not.all(ieee_is_finite(v1(:,j))).or..not.all(ieee_is_finite(v2(:,j))))status(j)=-4
      endif
      if(status(j)/=0.and.status(j)/=1)then
         v1(:,j)=ieee_value(0.0_dp,ieee_quiet_nan);v2(:,j)=v1(:,j)
      endif
   enddo

   ! The stages' steps, one call site each (internal, so the compiler inlines them).
   contains

      ! Transfer frame, normalized by the larger radius: unit vectors, the
      ! transverse directions, |cos(theta/2)|, |sin(theta/2)|, s1=sqrt(r1),
      ! s2=sqrt(r2), q=s1-s2 (without cancellation), a, b, y0=y(v=0), yend and
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
         dr=dot_product(r1-r2,r1+r2)/(g%rn1+g%rn2)     ! |r1|-|r2| without cancellation
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
         ! q=s1-s2; from |r1|-|r2| without cancellation when the radii are close.
         g%q=(g%rn1-g%rn2)/(g%s1+g%s2)
         if(min(g%rn1,g%rn2)>.5_dp)g%q=(dr/g%scale)/(g%s1+g%s2)
         g%yend=g%q*g%q+2*g%p*g%ss*g%ss/(1+g%cs)
         g%y0=g%yend;if(direction<0)g%y0=g%q*g%q+2*g%p*(1+g%cs)
         g%ypi=g%yend;if(g%b>=0)g%ypi=g%y0+2*g%b
         g%t=tof/(g%scale*sqrt(g%scale))
         status=-4
      end subroutine

      ! Zero-revolution starter in two steps. table_coordinates: k=0 for a
      ! parabola; k=1 (ellipse) or -1 (hyperbola, with the short-way branch
      ! point wb), the Lancaster-Blanchard tau and the table coordinate u
      ! ((tp/t)**(1/3) or t/tp), or within 3% of the parabolic time (near) the
      ! second-order series reversion X in u. table_root: the root x in
      ! coordinate k (2: pi-v near the upper pole); 0 when no starter applies.
      pure subroutine table_coordinates(g,k,tau,u,wb,near)
         type(frame),intent(in)::g
         integer,intent(out)::k
         real(dp),intent(out)::tau,u,wb
         logical,intent(out)::near
         real(dp)::tp,dc,sl,sw
         tp=sqrt(g%y0)*(g%y0*sqrt2/3+g%a)
         near=abs(g%t-tp)<=.03_dp*tp
         tau=0;u=0;wb=0;k=0
         if(abs(g%t-tp)<=4*epsilon(tp)*tp)return
         k=1;u=rcbrt(tp/g%t)
         if(g%t<tp)then
            k=-1;u=g%t/tp
            if(g%b>0)then
               sw=sqrt(g%y0/(2*g%b))
               wb=2*log1p_c(sw+sw*sw/(1+sqrt(1+sw*sw)))
            endif
         endif
         if(near)then
            u=near_parabola_x(g%y0,g%a,g%b,g%t,tp)
         else
            dc=sqrt(g%q*g%q*(g%s1+g%s2)**2+4*g%p*g%p*g%ss*g%ss);sl=(g%rn1+g%rn2+dc)/2   ! chord, semiperimeter
            tau=g%a/(sqrt2*sl)/(1+sqrt(dc/sl))
         endif
      end subroutine

      pure subroutine table_root(k,tau,u,wb,near,x)
         integer,intent(inout)::k
         real(dp),intent(in)::tau,u,wb
         logical,intent(in)::near
         real(dp),intent(out)::x
         x=0
         if(k==0)return
         if(near)then
            if(u*k>0)x=sqrt(abs(u))
         elseif(k>0)then
            x=table_v(tau,u)
            if(x>.75_dp*pi)then
               x=pi-x;k=2
            endif
         else
            x=table_w(tau,u,wb)
         endif
      end subroutine

      ! Coordinate, half-width el and half-angle functions of the centre for
      ! the one-rectangle path. The bound on the distance to every
      ! singularity of 1/(T-t) is min(v,pi-v), min(w,pi) or min(2w,wb-w), |X|
      ! near the parabola (X=+-x**2, a series for |X|<=0.2), or the bracket
      ! [lo,hi] ending at the minimum-time location for m>0. el is 1/16 of it
      ! for 8 nodes (error ~(1/32)**16); the zero-revolution elliptic table
      ! is within 0.4% of it, so there el is 1/256 of it for 3 nodes
      ! (error ~(1/512)**6).
      pure subroutine fast_coordinate(k,x,m,b,wb,lo,hi,kind1,c,el,eight,hs,hc)
         integer,intent(in)::k,m
         real(dp),intent(in)::x,b,wb,lo,hi
         integer,intent(out)::kind1
         real(dp),intent(out)::c,el,hs,hc
         logical,intent(out)::eight
         hs=0;hc=1;kind1=k;c=x;eight=.true.
         if(m>0)then
            el=min(c-lo,hi-c)
            call sincos_half(c,hs,hc)
         elseif(abs(k)==1.and.x<.3_dp)then
            kind1=0;c=k*x*x
            el=min(abs(c),.2_dp-abs(c))
            if(k<0.and.b>0)el=min(el,wb*wb-abs(c))
         elseif(k>0)then
            el=min(c,pi-c)/16;eight=.false.  ! k=2: x is already pi-v
            call sincos_half(c,hs,hc)
         else
            el=min(c,pi)
            if(b>0)el=min(2*c,wb-c)
            call shch(c/2,hs,hc)
         endif
         el=el/16
      end subroutine

      ! Certify the rectangle (sign change, estimate inside, error estimate at
      ! rounding) and move the root and its half-angle functions to it.
      pure subroutine fast_accept(kind1,c,kr,dd,de,sc,hs,hc,k,x,ok)
         integer,intent(in)::kind1,kr
         real(dp),intent(in)::c,dd,de
         logical,intent(in)::sc
         real(dp),intent(inout)::hs,hc
         integer,intent(inout)::k
         real(dp),intent(inout)::x
         logical,intent(out)::ok
         real(dp)::cn
         ok=sc.and.abs(dd)<.8_dp*pow2m(kr).and.de<=4*epsilon(c)*max(abs(c),merge(1e-2_dp,0._dp,kind1==0))
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
      ! sw, dw (or from Y=sqrt(y)=yy>=0 for the shortest hyperbolas): radial and
      ! transverse components, cancellation-free, or near the parabola and the
      ! short-way branch point (fine) Lagrange f,g in double-double.
      pure subroutine velocities(g,r1,r2,direction,k,x,sw,dw,fine,yy,v1,v2)
         type(frame),intent(in)::g
         real(dp),intent(in)::r1(3),r2(3),x,sw,dw
         integer,intent(in)::direction,k
         logical,intent(in)::fine
         real(dp),intent(out)::v1(3),v2(3)
         real(dp),intent(in)::yy
         real(dp)::p,c,s,y,dc,sh,ch,vr1,vr2,vt1,vt2
         logical::fg
         p=sign(1.0_dp,g%a);fg=fine;sh=sw;ch=dw
         if(yy>=0)then
            y=yy*yy;dc=-2*(g%y0-y)/(2*g%b)-g%ss*g%ss/(1+g%cs)   ! cs-cosh(w), from y
         else
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
            else
               y=g%y0-2*g%b*sh*sh
               ! Short way near the branch point, y cancels; T = sqrt(y)*(y*F+a) = t
               ! gives it with that error reduced by y*F/(y*F+a) when a
               ! dominates, and the half-angle functions follow from it
               ! without cancellation.
               if(g%b>0.and.y<g%y0/4)then
                  dc=y*fhyp(x,2*sh*ch,c)
                  if(dc<g%a)then
                     y=(g%t/(dc+g%a))**2
                     sh=sqrt((g%y0-y)/(2*g%b));ch=sqrt(1+sh*sh);c=1+2*sh*sh;s=-(2*sh*ch)**2
                     fg=.true.
                  endif
               endif
            endif
            dc=p*g%cs-c
            if(p*g%cs*c>0.and.abs(p*g%cs+c)>0.5_dp)dc=(s-g%ss*g%ss)/(p*g%cs+c)
         endif
         if(fg.and.g%cs>.3_dp)then
            ! Lagrange f,g from the exact inputs: correctly rounded for this y
            ! (radial/transverse needs the rounded unit vectors).
            call fg_velocities(r1,r2,y,g%scale,direction,v1,v2)
         else
            ! With q=s1-s2: (r2-r1)/(r1*(s2/s1+1)) = -q/s1, and vt uses s2/s1.
            y=sqrt2/sqrt(y*g%scale)                ! common factor sqrt(2/(y*scale))
            dc=y*dc;p=y*p*g%cs*g%q;s=y*g%ss
            vr1=dc-p/g%s1;vt1=s*g%s2/g%s1
            vr2=-(dc+p/g%s2);vt2=s*g%s1/g%s2
            v1=vr1*g%e1+vt1*g%t1v
            v2=vr2*g%e2+vt2*g%t2v
         endif
      end subroutine
   end subroutine solve_chunk

   ! Zero-height thin rectangle on [c-L,c+L], L=2**(-kr): Gauss-Legendre
   ! 8 (stage 1) or 3 (stage 2) nodes plus the two endpoints. dd=Q1/Q0
   ! locates the enclosed root; de=|Q2/Q0-dd**2|/L vanishes for a single
   ! pole and measures the error; sc certifies a sign change of T-t.
   ! Nodes (kind 0: X=k*x**2 series; 1: v; 2: pi-v; -1: w) are evaluated
   ! in one fused loop per coordinate; hs,hc are the half-angle functions
   ! of c and the tables those of L*u, combined by angle addition.
   pure subroutine rect_ratio(kind1,c,hs,hc,kr,eight,y0,yend,a,b,t,m,dd,de,sc)
      integer,intent(in)::kind1,kr,m
      real(dp),intent(in)::c,hs,hc,y0,yend,a,b,t
      logical,intent(in)::eight
      real(dp),intent(out)::dd,de
      logical,intent(out)::sc
      real(dp),dimension(n1+2)::uu,wq,ct,st,f,f1
      real(dp)::el,q0,q1,q2,r0,r1,r2,fa,fb,x,sh,ch,sn,cn,rs,yv,y1,ff,ff1,q,g,g1,ser,ser1,pk,y00,ab,mt
      integer::np,j,n
      logical::lw,us
      el=pow2m(kr)
      if(eight)then
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
      ! Q[g'] = g(c+L)-g(c-L) - L*sum w g'(x) for g = f, (x-c)f, (x-c)**2 f.
      ! Two interleaved partial sums shorten the dependent chain.
      fb=f(np-1);fa=f(np)
      q0=0;q1=0;q2=0;r0=0;r1=0;r2=0
      do j=1,np-3,2
         x=el*uu(j)
         q0=q0+wq(j)*f1(j)
         q1=q1+wq(j)*(x*f1(j)+f(j))
         q2=q2+wq(j)*x*(x*f1(j)+2*f(j))
         x=el*uu(j+1)
         r0=r0+wq(j+1)*f1(j+1)
         r1=r1+wq(j+1)*(x*f1(j+1)+f(j+1))
         r2=r2+wq(j+1)*x*(x*f1(j+1)+2*f(j+1))
      enddo
      if(mod(np,2)==1)then
         j=np-2;x=el*uu(j)
         q0=q0+wq(j)*f1(j)
         q1=q1+wq(j)*(x*f1(j)+f(j))
         q2=q2+wq(j)*x*(x*f1(j)+2*f(j))
      endif
      q0=fb-fa-el*(q0+r0)
      q1=el*(fb+fa)-el*(q1+r1)
      q2=el*el*(fb-fa)-el*(q2+r2)
      q=1/q0
      dd=q1*q
      de=abs(q2*q-dd*dd)*pow2m(-kr)
      sc=fa*fb<0                     ! T-t changes sign: a root is inside
      if(.not.(de==de).or..not.(dd==dd))then
         dd=huge(dd);de=huge(de)
      endif
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

   ! sin(c/2) and cos(c/2) for 0<=c<=pi without a library call: reduction
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

   ! asinh(s), s>=0, to ~1e-4 for singularity distances (only the rectangle
   ! width depends on them): series below 1/2, else a polynomial logarithm.
   elemental real(dp) function rasinh(s) result(y)
      real(dp),intent(in)::s
      real(dp)::q
      if(s<.5_dp)then
         q=s*s
         y=s*(1+q*(-1._dp/6+q*(3._dp/40+q*(-5._dp/112+q*(35._dp/1152)))))
      else
         y=qlog(s+sqrt(1+s*s))
      endif
   end function

   ! Natural logarithm to ~1e-5 absolute: exponent bits and a quintic in the
   ! mantissa (Abramowitz and Stegun 4.1.44); no division.
   elemental real(dp) function qlog(x) result(y)
      real(dp),intent(in)::x
      real(dp),parameter::ln2=log(2._dp)
      integer(int64)::i
      real(dp)::u
      if(.not.(ieee64.and.x>=tiny(x).and.x<=huge(x)))then
         y=log(x);return
      endif
      i=transfer(x,i)
      u=transfer(ior(iand(i,4503599627370495_int64),4607182418800017408_int64),u)-1
      y=(shiftr(i,52)-1023)*ln2+u*(.99949556_dp+u*(-.49190896_dp+u*(.28947478_dp+u*(-.13606275_dp+u*.03215845_dp))))
   end function

   ! asin(x), 0<=x<=1, to ~5e-5 (Abramowitz and Stegun 4.4.45).
   elemental real(dp) function rasin(x) result(y)
      real(dp),intent(in)::x
      y=pi/2-sqrt(1-x)*(1.5707288_dp+x*(-.2121144_dp+x*(.0742610_dp-x*.0187293_dp)))
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

   ! Cube root to ~2e-5 for the closed-form starters (the contour sets the
   ! accuracy): exponent/3 bit estimate (3%) and one Halley correction; a
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

   ! x from a rectangle coordinate: X=k*x**2 near the parabola (series), else x
   ! itself; -1 when X lies on the other side of the parabola.
   elemental real(dp) function xof(c,series,k) result(x)
      real(dp),intent(in)::c
      logical,intent(in)::series
      integer,intent(in)::k
      x=c
      if(series)then
         x=-1
         if(k*c>=0)x=sqrt(k*c)
      endif
   end function

   ! Zero-revolution root v from the tables at tau=lam/(1+sqrt(1-lam**2)) and
   ! phi=(tp/t)**(1/3): A=(pi-v)/phi for phi<=1/2, B=v**2/(1-phi) above, by
   ! local 4x4 cubic Lagrange interpolation. 0 outside the tables.
   pure real(dp) function table_v(tau,phi) result(v)
      real(dp),intent(in)::tau,phi
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
         v=pi-phi*u
      else
         v=sqrt(max(u*(1-phi),0._dp))
      endif
      if(.not.(v>0.and.v<pi))v=0
   end function

   ! Zero-revolution hyperbolic root w from the tables at tau and u=t/tp
   ! (see excali_tables.inc); wb is the short-way branch point. 0 outside.
   pure real(dp) function table_w(tau,u,wb) result(w)
      real(dp),intent(in)::tau,u,wb
      real(dp)::g,p,wt(4),ws(4),r(4)
      integer::i,j,h
      w=0
      if(.not.(abs(tau)>=tabh_tmin.and.abs(tau)<=tabh_tmax.and.u<=.99_dp))return
      if(tau>0)then
         h=1;p=u
         if(.not.(p>=tabh_ps))return
         g=(p-tabh_ps)*((tab_n-1)/(1-tabh_ps))
      else
         h=2;p=rcbrt(u)
         if(.not.(p>=tabh_pl))return
         g=(p-tabh_pl)*((tab_n-1)/(1-tabh_pl))
      endif
      j=min(max(int(g)-1,0),tab_n-4);call lagrange4(g-j,ws)
      g=(abs(tau)-tabh_tmin)*((tab_n-1)/(tabh_tmax-tabh_tmin))
      if(h==2)g=(tab_n-1)-g
      i=min(max(int(g)-1,0),tab_n-4);call lagrange4(g-i,wt)
      r=wt(1)*tabh(i+1,j+1:j+4,h)+wt(2)*tabh(i+2,j+1:j+4,h)+wt(3)*tabh(i+3,j+1:j+4,h)+wt(4)*tabh(i+4,j+1:j+4,h)
      g=(ws(1)*r(1)+ws(2)*r(2))+(ws(3)*r(3)+ws(4)*r(4))
      if(h==1)then
         w=sqrt(max(g*(1-u),0._dp))*wb
         if(.not.(w<wb))w=0
      else
         w=sqrt(max(g*(1-u),0._dp))*(1-rlog(u))
      endif
      if(.not.(w>0))w=0
   end function

   ! Near the parabola, X = k*x**2 (v**2 for an ellipse, -w**2 for a
   ! hyperbola) from the series T = T0 + T1*X + T2*X**2 + ..., with
   ! y = y0 + b*X*(1/2 - X/24 + ...) and y*F + a from the parabolic series of
   ! F, reverted to second order: X = d/T1 - T2*d**2/T1**3, d = t - T0.
   pure real(dp) function near_parabola_x(y0,a,b,t,tp) result(xx)
      real(dp),intent(in)::y0,a,b,t,tp
      real(dp)::q0,q1,q2,g0,g1,g2,p0,p1,p2,t1,t2,d
      g0=pser2(10);g1=pser2(9);g2=pser2(8)
      q0=sqrt(y0);q1=b/(4*q0);q2=(-b/24-q1*q1)/(2*q0)
      p0=y0*g0+a;p1=y0*g1+b/2*g0;p2=y0*g2+b/2*g1-b/24*g0
      t1=q0*p1+q1*p0;t2=q0*p2+q1*p1+q2*p0
      d=t-tp
      xx=d/t1-t2*d*d/(t1*t1*t1)
   end function

   ! Multi-revolution table starter: tau and the semiperimeter sl, time t,
   ! revolutions m and branch br. stat=5: t is below the tabulated minimum
   ! time by more than the table error (infeasible); stat=1: v is the
   ! branch root estimate and vs the minimum-time location; stat=0: outside
   ! the tables or too close to t_min. Slices in mu=1/m by 4-point Lagrange.
   pure subroutine mrev_table(tau,sl,t,m,br,v,vs,stat)
      real(dp),intent(in)::tau,sl,t
      integer,intent(in)::m,br
      real(dp),intent(out)::v,vs
      integer,intent(out)::stat
      real(dp),parameter::margin=1e-6_dp
      real(dp)::u,wt(4),wc(4),ws(4),wm(4),mu,cm,r,sig
      integer::i,ic,j,k,km(4),nk,ia,jb,kk
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
      if(r>=1-margin)return
      sig=br*sqrt(1-r)
      u=(tau+tab_tmax)*((mr_nt-1)/(2*tab_tmax))
      i=min(max(int(u)-1,0),mr_nt-4);call lagrange4(u-i,wt)
      u=(sig+1)*((mr_ns-1)/2._dp)
      j=min(max(int(u)-1,0),mr_ns-4);call lagrange4(u-j,ws)
      do k=1,nk
         do jb=1,4
            v=v+wm(k)*ws(jb)*(wt(1)*mr_v(i+1,j+jb,km(k))+wt(2)*mr_v(i+2,j+jb,km(k)) &
               +wt(3)*mr_v(i+3,j+jb,km(k))+wt(4)*mr_v(i+4,j+jb,km(k)))
         enddo
      enddo
      if(v>0.and.v<pi.and.vs>0.and.vs<pi.and.(v-vs)*br>0)stat=1
   end subroutine

   ! Hyperbolic F(w) = (sinh(w)*cosh(w)-w)/(sqrt2*sinh(w)**3) from sinh(w) and
   ! cosh(w); the parabolic series below w = 0.3, where the closed form cancels.
   pure real(dp) function fhyp(w,sh,ch) result(f)
      real(dp),intent(in)::w,sh,ch
      real(dp)::xx
      integer::j
      if(w<.3_dp)then
         xx=-w*w;f=pser2(1)
         do j=2,10
            f=f*xx+pser2(j)
         enddo
      else
         f=(sh*ch-w)*rsqrt2/sh**3
      endif
   end function

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
      pure subroutine fsum(a,b,s,e)            ! s+e = a+b exactly
         real(dp),intent(in)::a,b
         real(dp),intent(out)::s,e
         real(dp)::v
         s=a+b;v=s-a;e=(a-(s-v))+(b-v)
      end subroutine
      pure subroutine dot2(u,v,h,l)            ! u.v as h+l (compensated)
         real(dp),intent(in)::u(3),v(3)
         real(dp),intent(out)::h,l
         real(dp)::p,q,e,w
         integer::i
         h=u(1)*v(1);l=ieee_fma(u(1),v(1),-h)
         do i=2,3
            p=u(i)*v(i);e=ieee_fma(u(i),v(i),-p)
            call fsum(h,p,w,q);h=w;l=l+(q+e)
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

   !> Cubic Hermite on [0,1]: values p0,p1 and slopes m0,m1 at 0 and 1.
   pure real(dp) function hermite(u,p0,m0,p1,m1)
   real(dp),intent(in) :: u,p0,m0,p1,m1
   hermite=p0*(2*u**3-3*u*u+1)+m0*(u**3-2*u*u+u)+p1*(3*u*u-2*u**3)+m1*(u**3-u*u)
   end function hermite

   pure function cross(u,v) result(w)
   real(dp),intent(in) :: u(3),v(3)
   real(dp) :: w(3)
   w=[u(2)*v(3)-u(3)*v(2),u(3)*v(1)-u(1)*v(3),u(1)*v(2)-u(2)*v(1)]
   end function cross

end module excali_solver
