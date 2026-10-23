! eLam: finite rectangular-contour Lambert root extraction across conics.
! Extended from Exact v2.2's thin-rectangle ratio. The root comes from the
! zero-height limit of that rectangle (real Gauss-Legendre quadrature of
! 1/(T-t) and its derivative), applied in one or two fixed stages without a
! Newton polish. Geometry, half-anomaly time functions, derivatives and
! closed-form starters are adapted from this project's contour_simple. No
! confluent/Householder root steps and no calls to another Lambert solver.
module elam_solver
   use iso_fortran_env, only: real64, int64
   use ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use iso_c_binding, only: c_double
   implicit none
   private
   public :: elam_one, elam_solve, rect_ratio
   integer,parameter :: dp=real64
   real(dp),parameter :: pi=acos(-1.0_dp),sqrt2=sqrt(2.0_dp),rsqrt2=1/sqrt2
   integer,parameter  :: max_steps=120
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
   real(dp),parameter :: u8(n1+2)=[g8x(1:h1),-g8x(1:h1),1._dp,-1._dp],w8(n1+2)=[g8w,g8w,0._dp,0._dp]
   real(dp),parameter :: u3(5)=[0._dp,g3x(2),-g3x(2),1._dp,-1._dp],w3(5)=[g3w(1),g3w(2),g3w(2),0._dp,0._dp]
   real(dp),parameter :: c8f(n1+2,0:kmax)=reshape([((cos(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   real(dp),parameter :: s8f(n1+2,0:kmax)=reshape([((sin(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   real(dp),parameter :: ch8f(n1+2,0:kmax)=reshape([((cosh(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   real(dp),parameter :: sh8f(n1+2,0:kmax)=reshape([((sinh(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   real(dp),parameter :: c3f(5,0:kmax)=reshape([((cos(.5_dp**(k_+1)*u3(j_)),j_=1,5),k_=0,kmax)],[5,kmax+1])
   real(dp),parameter :: s3f(5,0:kmax)=reshape([((sin(.5_dp**(k_+1)*u3(j_)),j_=1,5),k_=0,kmax)],[5,kmax+1])
   real(dp),parameter :: ch3f(5,0:kmax)=reshape([((cosh(.5_dp**(k_+1)*u3(j_)),j_=1,5),k_=0,kmax)],[5,kmax+1])
   real(dp),parameter :: sh3f(5,0:kmax)=reshape([((sinh(.5_dp**(k_+1)*u3(j_)),j_=1,5),k_=0,kmax)],[5,kmax+1])
   ! Reciprocal factorials 1/j! (constant multiplies, not divisions).
   real(dp),parameter :: rf(14)=[(1/gamma(j_+1._dp),j_=1,14)]
   real(dp),parameter :: c8(h1+1,0:kmax)=reshape([((cos(.5_dp**(k_+1)*g8x(j_)),j_=1,h1+1),k_=0,kmax)],[h1+1,kmax+1])
   real(dp),parameter :: s8(h1+1,0:kmax)=reshape([((sin(.5_dp**(k_+1)*g8x(j_)),j_=1,h1+1),k_=0,kmax)],[h1+1,kmax+1])
   real(dp),parameter :: ch8(h1+1,0:kmax)=reshape([((cosh(.5_dp**(k_+1)*g8x(j_)),j_=1,h1+1),k_=0,kmax)],[h1+1,kmax+1])
   real(dp),parameter :: sh8(h1+1,0:kmax)=reshape([((sinh(.5_dp**(k_+1)*g8x(j_)),j_=1,h1+1),k_=0,kmax)],[h1+1,kmax+1])
   real(dp),parameter :: c3(3,0:kmax)=reshape([((cos(.5_dp**(k_+1)*g3x(j_)),j_=1,3),k_=0,kmax)],[3,kmax+1])
   real(dp),parameter :: s3(3,0:kmax)=reshape([((sin(.5_dp**(k_+1)*g3x(j_)),j_=1,3),k_=0,kmax)],[3,kmax+1])
   real(dp),parameter :: ch3(3,0:kmax)=reshape([((cosh(.5_dp**(k_+1)*g3x(j_)),j_=1,3),k_=0,kmax)],[3,kmax+1])
   real(dp),parameter :: sh3(3,0:kmax)=reshape([((sinh(.5_dp**(k_+1)*g3x(j_)),j_=1,3),k_=0,kmax)],[3,kmax+1])
   ! (1-cos(v))/X in X=v**2 (highest first), for the near-parabolic coordinate.
   real(dp),parameter :: cser(9)=[1._dp/6402373705728000._dp,-1._dp/20922789888000._dp,1._dp/87178291200._dp, &
      -1._dp/479001600._dp,1._dp/3628800._dp,-1._dp/40320._dp,1._dp/720._dp,-1._dp/24._dp,.5_dp]
   ! Parabolic series of (v-sin(v)cos(v))/sin(v)**3 in v**2 (highest first),
   ! used for |v|<0.3 on zero-revolution arcs where the closed form cancels.
   real(dp),parameter :: pser(10)=[115665628927._dp/2601853536522240000._dp, &
      3821869001._dp/10751460894720000._dp,133541._dp/48117888000._dp,763967._dp/36324288000._dp, &
      1393481._dp/9081072000._dp,1181._dp/1108800,29._dp/4200,17._dp/420,1._dp/5,2._dp/3]
   real(dp),parameter :: pser2(10)=pser/sqrt(2._dp)

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

   subroutine elam_solve(r1,r2,tof,mu,direction,nrev,branch,v1,v2,status,plane_normal)
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
      call elam_one(r1/sc,r2/sc,tt,direction,v1,v2,status,half,nrev,branch,plane_normal)
      if(status==0.or.status==1)then
         v1=v1*vscale;v2=v2*vscale
      endif
   end subroutine

   subroutine elam_one(r1,r2,tof,direction,v1,v2,status,half,nrev,branch,plane_normal)
   real(dp),intent(in)          :: r1(3),r2(3),tof
   integer,intent(in)           :: direction
   real(dp),intent(out)         :: v1(3),v2(3)
   integer,intent(out)          :: status,half
   integer,optional,intent(in)  :: nrev,branch
   real(dp),optional,intent(in) :: plane_normal(3)
   real(dp) :: scale,rn1,rn2,e1(3),e2(3),h(3),hn,t1v(3),t2v(3),cs,ss,p,q,a,b,y0,yend,ypi, &
               t,tp,tv2,x,lo,hi,phi,sig,tmin,d(0:3),f(0:3),k0,kpi,c,s,y,dc,ratio, &
               vr1,vr2,vt1,vt2,wb,psi,yy,sw,dt,hs,hc
   integer  :: m,br,k,it
   logical :: merged,trig

   v1=ieee_value(0.0_dp,ieee_quiet_nan);v2=v1;half=0;status=-1
   merged=.false.;trig=.false.
   m=0;br=-1
   if(present(nrev))m=nrev
   if(present(branch))br=branch
   if(abs(direction)/=1.or.m<0.or.abs(br)/=1.or..not.(tof>0).or..not.ieee_is_finite(tof))return

   !======= Geometry, normalized by the larger radius =======
   ! Normalized entry (elam_solve pre-scales): plain dot products suffice.
   rn1=sqrt(dot_product(r1,r1));rn2=sqrt(dot_product(r2,r2))
   if(.not.(rn1>0.and.rn2>0).or..not.ieee_is_finite(rn1+rn2))return
   scale=max(rn1,rn2);e1=r1/rn1;e2=r2/rn2;rn1=rn1/scale;rn2=rn2/scale
   h=cross(e1,e2);hn=sqrt(dot_product(h,h))
   status=-2
   if(hn<32*epsilon(hn))then
      half=-1
      if(.not.present(plane_normal))return
      if(norm2(e1+e2)>32*epsilon(hn))then
         status=-3;return
      endif
      if(.not.all(ieee_is_finite(plane_normal)))return
      h=plane_normal;hn=norm2(h)
      if(hn<=0.or..not.ieee_is_finite(hn))return
      if(max(abs(dot_product(h/hn,e1)),abs(dot_product(h/hn,e2)))>64*epsilon(hn))return
      half=1
   endif
   if(hn<1e-4_dp)half=1
   h=h*(direction/hn)
   t1v=e1+e2;t2v=e2-e1
   cs=sqrt(dot_product(t1v,t1v))/2;ss=sqrt(dot_product(t2v,t2v))/2 ! |cos(theta/2)|, |sin(theta/2)|
   t1v=cross(h,e1);t2v=cross(h,e2)
   p=sqrt(rn1*rn2);a=direction*sqrt2*p*cs;b=sqrt2*a
   q=(rn1-rn2)/(sqrt(rn1)+sqrt(rn2))
   ! y(v)=y0+2b*sin(v/2)**2 (b>=0) or yend-2b*cos(v/2)**2 (b<0); ypi=y(pi).
   yend=q*q+2*p*ss*ss/(1+cs)
   y0=yend;if(direction<0)y0=q*q+2*p*(1+cs)   ! y at v=0
   ypi=yend;if(b>=0)ypi=y0+2*b
   t=tof/(scale*sqrt(scale))
   status=-4

   !======= Starters and admissible interval =======
   k=1;lo=0;hi=pi;yy=-1
   if(m==0)then
      tp=sqrt(y0)*(y0*sqrt2/3+a)                                 ! parabolic time
      tv2=y0*sqrt(y0)/(5*sqrt2)+b*(sqrt(y0)/sqrt2+a/(2*sqrt(y0)))/2 ! dT/d(v**2) at 0
      if(abs(t-tp)<=4*epsilon(t)*tp)then
         x=0                                                     ! parabola
      elseif(t>tp)then
         ! Ellipse: cubic Hermite for v**2 in phi=(tp/t)**(1/3), between
         ! v**2=(t-tp)/tv2 at phi=1 and pi-v=phi*(K/tp)**(1/3) at phi=0,
         ! with the exact upper pole T ~ pi*y(pi)**1.5/(sqrt2*(pi-v)**3).
         phi=rcbrt(tp/t)
         x=hermite(phi,pi*pi,-2*pi*sqrt(ypi)*rcbrt(pi/sqrt2/tp),0._dp,-3*tp/tv2)
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
               ! Newton in Y=sqrt(y): T=Y*(Y**2*F+a) with the exact
               ! w(Y)=2*asinh(sqrt((y0-Y**2)/(2b))). Either term of
               ! T = a*Y + F*Y**3 alone overestimates Y; start from the smaller.
               call f_derivs(wb,-1,0,f)
               yy=min(t/a,rcbrt(t/f(0)),0.9_dp*sqrt(y0))
               ! A separate thin rectangle in Y=sqrt(y) avoids the
               ! branch-point cancellation of hyperbolic anomaly space.
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
      c=(rn1+rn2+sqrt((rn1-rn2)**2+4*rn1*rn2*ss*ss))/2
      if(t<m*pi*c*sqrt(c)/sqrt2*(1-1e-12_dp))return
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
   endif

   ! Store distance from the upper pole directly on near-full-period arcs.
   ! Subtracting a converged v from pi loses the small arc information.
   if(k==1.and.x>0.75_dp*pi.and..not.merged)then
      p=lo;lo=pi-hi;hi=pi-p;x=pi-x;k=2
   endif

   !======= Root: zero-height thin rectangles =======
   if(x>0.and.yy<0.and..not.merged)then
      ! Rare starter failures fall back to the safeguarded bracket iteration.
      if(.not.contour_root())then
         if(.not.converge())return
      endif
   endif

   !======= Velocities (radial/transverse, cancellation-free) =======
   p=sign(1.0_dp,a)
   if(yy>=0)then
      y=yy*yy;dc=-2*sw-ss*ss/(1+cs)             ! cs-cosh(w), from y
   else
      ! One pair of half-angle functions gives c, s=sin(x)**2 and y.
      ! The contour carries them along from its starter (trig).
      if(trig)then
         sw=hs;dt=hc
      elseif(k>0)then
         sw=sin(x/2);dt=cos(x/2)
      else
         sw=sinh(x/2);dt=cosh(x/2)
      endif
      if(k>0)then
         c=(dt-sw)*(dt+sw);s=(2*sw*dt)**2
      else
         c=1+2*sw*sw;s=-(2*sw*dt)**2
      endif
      if(k==2)then
         c=-c
         y=y0+2*b*dt*dt
         if(b<0)y=yend-2*b*sw*sw
      elseif(k>0)then
         if(b>=0)then
            y=y0+2*b*sw*sw
         else
            y=yend-2*b*dt*dt
         endif
      else
         y=y0-2*b*sw*sw
      endif
      dc=p*cs-c
      if(p*cs*c>0.and.abs(p*cs+c)>0.5_dp)dc=(s-ss*ss)/(p*cs+c)
   endif
   ratio=sqrt(rn2/rn1)
   y=sqrt2/sqrt(y*scale)                     ! common factor sqrt(2/(y*scale))
   vr1=y*(dc+p*cs*(rn2-rn1)/(rn1*(ratio+1)))
   vr2=-y*(dc+p*cs*(rn1-rn2)/(rn2*(1/ratio+1)))
   vt1=y*ratio*ss;vt2=y*ss/ratio
   v1=vr1*e1+vt1*t1v
   v2=vr2*e2+vt2*t2v
   status=0
   if(merged)status=1
   if(.not.all(ieee_is_finite(v1)).or..not.all(ieee_is_finite(v2)))status=-4

   contains

      ! Full finite upper-half rectangular contour, 8-point Gauss-Legendre
      ! on the top and one midpoint on each vertical side, as in Exact v2.2.
      ! The ratio uses the imaginary integrals; lower-half symmetry cancels
      ! real parts. Centering the numerator avoids loss of root increments.
      real(dp) function rectangle(center,radius,in_y) result(root)
         real(dp),intent(in)::center,radius
         logical,intent(in)::in_y
         real(dp),parameter::nodes(8)=[-.9602898564975363_dp,-.7966664774136267_dp, &
            -.5255324099163290_dp,-.1834346424956498_dp,.1834346424956498_dp, &
             .5255324099163290_dp,.7966664774136267_dp,.9602898564975363_dp]
         real(dp),parameter::weights(8)=[.1012285362903763_dp,.2223810344533745_dp, &
            .3137066458778873_dp,.3626837833783620_dp,.3626837833783620_dp, &
            .3137066458778873_dp,.2223810344533745_dp,.1012285362903763_dp]
         complex(dp)::z,w,den,num,fn
         real(dp)::height
         integer::j
         den=0;num=0
         height=min(1e-7_dp,0.01_dp*radius)
         do j=1,10
            if(j<=8)then
               z=cmplx(center-radius*nodes(j),height,dp)
               w=cmplx(-radius*weights(j),0._dp,dp)
            elseif(j==9)then
               z=cmplx(center+radius,height/2,dp);w=cmplx(0._dp,height,dp)
            else
               z=cmplx(center-radius,height/2,dp);w=cmplx(0._dp,-height,dp)
            endif
            fn=complex_time(z,in_y)-t
            if(abs(fn)==0)then
               root=center;return
            endif
            w=w/fn;den=den+w;num=num+(z-center)*w
         enddo
         root=center
         if(abs(aimag(den))>tiny(t))root=center+aimag(num)/aimag(den)
      end function

      complex(dp) function complex_time(z,in_y) result(tt)
         complex(dp),intent(in)::z
         logical,intent(in)::in_y
         complex(dp)::zz,ss1,cc1,yy1,ff,xx1
         integer::j,kind1
         zz=z;kind1=k
         if(in_y)then
            ! Analytic asinh using its logarithmic identity.
            xx1=sqrt((y0-z*z)/(2*b))
            zz=2*log(xx1+sqrt(1+xx1*xx1));kind1=-1
         endif
         if(kind1==2)then
            ss1=sin(zz);cc1=cos(zz)
            yy1=y0+2*b*cos(zz/2)**2
            if(b<0)yy1=yend-2*b*sin(zz/2)**2
         elseif(kind1>0)then
            ss1=sin(zz);cc1=cos(zz)
            yy1=y0+2*b*sin(zz/2)**2
            if(b<0)yy1=yend-2*b*cos(zz/2)**2
         else
            ss1=sinh(zz);cc1=cosh(zz);yy1=y0-2*b*sinh(zz/2)**2
         endif
         if(kind1==2)then
            ff=((m+1)*pi-zz+ss1*cc1)/(sqrt2*ss1**3)
         elseif(m==0.and.abs(zz)<0.3_dp)then
            xx1=kind1*zz*zz;ff=pser(1)
            do j=2,10
               ff=ff*xx1+pser(j)
            enddo
            ff=ff/sqrt2
         else
            ff=kind1*(zz+m*pi-ss1*cc1)/(sqrt2*ss1**3)
         endif
         if(in_y)then
            tt=z*(z*z*ff+a)
         elseif(kind1<0.and.b<0.and.abs(zz)>=0.3_dp)then
            ! Cancel the large asymptotic terms symbolically, not numerically.
            tt=sqrt(yy1)/sqrt2*(((y0+b)*cc1-b)/ss1**2-yy1*zz/ss1**3)
         else
            tt=sqrt(yy1)*(yy1*ff+a)
         endif
      end function

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
               bim=2*asin(sqrt(y0/(-2*b)))
            endif
         endif
         c=x;kind1=k;nox=.false.;de=0
         ! Stage 1 starts at L~R/4, where 8 nodes usually reach rounding
         ! level. A self-consistent estimate outside the interval re-centres
         ! (going straight to stage 2 if already precise enough); an
         ! inconsistent one (root far away) doubles L, up to ~R.
         grow=4;outside=.false.
         do pass=1,5
            call canonical(kind1,c,de,nox)
            if(.not.trig.and.kind1>0)then
               hs=sin(c/2);hc=cos(c/2);trig=.true.
            elseif(.not.trig.and.kind1<0)then
               hs=sinh(c/2);hc=cosh(c/2);trig=.true.
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
               hs=sin(c/2);hc=cos(c/2);trig=.true.
            elseif(.not.trig.and.kind1<0)then
               hs=sinh(c/2);hc=cosh(c/2);trig=.true.
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



      logical function converge_y() result(ok)
         real(dp)::res,der,trial,rad,old
         integer::j
         ok=.false.
         do j=1,max_steps
            res=real(complex_time(cmplx(yy,0._dp,dp),.true.),dp)-t
            if(abs(res)<=16*epsilon(t)*t)then
               ok=.true.;return
            endif
            if(res<0)then
               lo=yy
            else
               hi=yy
            endif
            der=aimag(complex_time(cmplx(yy,1e-24_dp*max(yy,1e-100_dp),dp),.true.))/(1e-24_dp*max(yy,1e-100_dp))
            trial=yy-res/der
            rad=min(0.45_dp*(hi-lo),max(2*abs(trial-yy),1e-5_dp*abs(yy)))
            if(j==1.and.rad>0.and.yy-rad>0.and.yy+rad<sqrt(y0))then
               old=rectangle(yy,rad,.true.)
               if(ieee_is_finite(old).and.old>lo.and.old<hi)trial=old
            endif
            if(.not.ieee_is_finite(trial).or.trial<=lo.or.trial>=hi)trial=(lo+hi)/2
            if(abs(trial-yy)<=8*epsilon(t)*max(abs(yy),tiny(t)))then
               ok=abs(res)<=1e-12_dp*t;return
            endif
            yy=trial
         enddo
      end function

      logical function converge() result(ok)
         real(dp)::res,trial,rad,old,sgn,left,right
         integer::j
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
         do j=1,max_steps
            call time_derivs(x,k,m,y0,yend,a,b,d)
            res=d(0)-t
            if(abs(res)<=16*epsilon(t)*t)then
               ok=.true.;return
            endif
            left=lo;right=hi
            if(res*sgn<0)then
               lo=x
            else
               hi=x
            endif
            trial=x-res/d(1)
            ! The first finite rectangle encloses the Newton prediction;
            ! all subsequent corrections remain bracketed on this branch.
            rad=min(.45_dp*(right-left),max(2*abs(trial-x),1e-5_dp*abs(x)))
            rad=min(rad,.9_dp*min(x-left,right-x))
            if(j==1.and.rad>0)then
               old=rectangle(x,rad,.false.)
               if(ieee_is_finite(old).and.old>lo.and.old<hi)trial=old
            endif
            if(.not.ieee_is_finite(trial).or.trial<=lo.or.trial>=hi)trial=(lo+hi)/2
            if(abs(trial-x)<=8*epsilon(t)*max(abs(x),1e-100_dp))then
               ok=abs(res)<=1e-12_dp*t;return
            endif
            x=trial
         enddo
      end function

   end subroutine elam_one

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
      real(dp)::el,q0,q1,q2,r0,r1,r2,fa,fb,x,sh,ch,sn,cn,rs,yv,y1,ff,ff1,q,g,g1,ser,ser1,pk,y00,ab,mt,sg
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
         sg=0;if(lw)sg=1
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
         s=sinh(h);c=cosh(h)
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
            hs=sinh(x/2);hc=cosh(x/2)
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
         s=sinh(x);c=cosh(x)
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
      sh=sinh(x/2);ch=cosh(x/2);c=1+2*sh*sh
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

   ! 2**(-k) and 1-exponent(x) (x>0) from the bit pattern, without calls.
   elemental real(dp) function pow2m(k) result(y)
      integer,intent(in)::k
      y=transfer(shiftl(int(1023-k,int64),52),y)
   end function

   elemental integer function kexp(x) result(k)
      real(dp),intent(in)::x
      k=1023-int(shiftr(transfer(x,0_int64),52))
   end function

   ! asinh(s), s>=0, to ~2e-5 for singularity distances: series below 1/2,
   ! else log(s+sqrt(1+s**2)) with a fast logarithm.
   elemental real(dp) function rasinh(s) result(y)
      real(dp),intent(in)::s
      real(dp)::q
      if(s<.5_dp)then
         q=s*s
         y=s*(1+q*(-1._dp/6+q*(3._dp/40+q*(-5._dp/112+q*(35._dp/1152)))))
      else
         y=rlog(s+sqrt(1+s*s))
      endif
   end function

   ! Natural logarithm to ~1e-9 relative: exponent bits, mantissa reduced to
   ! [sqrt(1/2),sqrt(2)], and the atanh series in z=(m-1)/(m+1), |z|<0.18.
   elemental real(dp) function rlog(x) result(y)
      real(dp),intent(in)::x
      real(dp),parameter::ln2=log(2._dp)
      integer(int64)::i,e
      real(dp)::m,z,q
      if(.not.(x>=tiny(x).and.x<=huge(x)))then
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
   ! accuracy): exponent/3 bit estimate (3%) and one Halley correction.
   elemental real(dp) function rcbrt(x) result(y)
      real(dp),intent(in)::x
      real(dp)::w,r
      if(.not.(x>=tiny(x).and.x<=huge(x)))then
         y=cbrt_c(x);return
      endif
      y=transfer(transfer(x,0_int64)/3+3071306045661749710_int64,y)
      w=1/y;r=(x*w)*(w*w)
      y=y*(1+2*r)/(2+r)
   end function

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

end module elam_solver
