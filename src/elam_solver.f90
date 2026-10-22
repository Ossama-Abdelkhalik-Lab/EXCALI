! eLam: finite rectangular-contour Lambert root extraction across conics.
! Extended from Exact v2.2's rectangle-ratio/Newton strategy. Geometry,
! half-anomaly time functions, derivatives and starters are adapted from
! this project's contour_simple. No confluent/Householder root steps and
! no calls to another Lambert solver are used.
module elam_solver
   use iso_fortran_env, only: real64
   use ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: elam_one, elam_solve
   integer,parameter :: dp=real64
   real(dp),parameter :: pi=acos(-1.0_dp),sqrt2=sqrt(2.0_dp)
   integer,parameter  :: max_steps=120
   ! Parabolic series of (v-sin(v)cos(v))/sin(v)**3 in v**2 (highest first),
   ! used for |v|<0.3 on zero-revolution arcs where the closed form cancels.
   real(dp),parameter :: pser(10)=[115665628927._dp/2601853536522240000._dp, &
      3821869001._dp/10751460894720000._dp,133541._dp/48117888000._dp,763967._dp/36324288000._dp, &
      1393481._dp/9081072000._dp,1181._dp/1108800,29._dp/4200,17._dp/420,1._dp/5,2._dp/3]

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
      tt=exp(log(tof)+.5_dp*log(mu)-1.5_dp*log(sc))
      vscale=exp(.5_dp*(log(mu)-log(sc)))
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
   real(dp) :: scale,rn1,rn2,e1(3),e2(3),h(3),hn,t1v(3),t2v(3),cs,ss,p,q,a,b,y0,yend, &
               t,tp,tv2,x,lo,hi,phi,sig,tmin,d(0:3),f(0:3),k0,kpi,c,s,y,dc,ratio, &
               vr1,vr2,vt1,vt2,wb,psi,yy,sw,dt
   integer  :: m,br,k,it
   logical :: merged

   v1=ieee_value(0.0_dp,ieee_quiet_nan);v2=v1;half=0;status=-1
   merged=.false.
   m=0;br=-1
   if(present(nrev))m=nrev
   if(present(branch))br=branch
   if(abs(direction)/=1.or.m<0.or.abs(br)/=1.or..not.(tof>0).or..not.ieee_is_finite(tof))return

   !======= Geometry, normalized by the larger radius =======
   rn1=norm2(r1);rn2=norm2(r2)
   if(.not.(rn1>0.and.rn2>0).or..not.ieee_is_finite(rn1+rn2))return
   scale=max(rn1,rn2);e1=r1/rn1;e2=r2/rn2;rn1=rn1/scale;rn2=rn2/scale
   h=cross(e1,e2);hn=norm2(h)
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
   h=h*(direction/hn);t1v=cross(h,e1);t2v=cross(h,e2)
   cs=norm2(e1+e2)/2;ss=norm2(e2-e1)/2        ! |cos(theta/2)|, |sin(theta/2)|
   p=sqrt(rn1*rn2);a=direction*sqrt2*p*cs;b=sqrt2*a
   q=(rn1-rn2)/(sqrt(rn1)+sqrt(rn2))
   yend=q*q+2*p*ss*ss/(1+cs)                  ! y at v=pi (elliptic, any direction)
   y0=yend;if(direction<0)y0=q*q+2*p*(1+cs)   ! y at v=0
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
         ! with the exact upper pole T ~ pi*yend**1.5/(sqrt2*(pi-v)**3).
         phi=(tp/t)**(1._dp/3)
         x=hermite(phi,pi*pi,-2*pi*sqrt(yend)*(pi/sqrt2/tp)**(1._dp/3),0._dp,-3*tp/tv2)
         if(.not.(x>0))x=(t-tp)/tv2         ! model undershoot: parabolic limit
         x=min(sqrt(x),pi*(1-1e-3_dp))
      else
         ! Hyperbola in w=sqrt(-z)/2: w**2=(tp-t)/tv2 near the parabola.
         k=-1;hi=huge(hi)
         if(b>0)then
            ! Short way: T -> 0 at the branch point y(wb)=0.
            wb=2*asinh(sqrt(y0/(2*b)));hi=wb
            if(t<0.1_dp*tp)then
               ! Newton in Y=sqrt(y): T=Y*(Y**2*F+a) with the exact
               ! w(Y)=2*asinh(sqrt((y0-Y**2)/(2b))). Either term of
               ! T = a*Y + F*Y**3 alone overestimates Y; start from the smaller.
               call f_derivs(wb,-1,0,f)
               yy=min(t/a,(t/f(0))**(1._dp/3),0.9_dp*sqrt(y0))
               ! A separate thin rectangle in Y=sqrt(y) avoids the
               ! branch-point cancellation of hyperbolic anomaly space.
               lo=0;hi=sqrt(y0)
               if(.not.converge_y())return
               sw=(y0-yy*yy)/(2*b);x=2*asinh(sqrt(sw))
            else
               ! Cubic Hermite for w**2 in (t/tp)**2 between the parabola and
               ! wb, where T**2 ~ a**2*y has slope -a**2*b*sinh(wb).
               x=hermite((t/tp)**2,wb*wb,-2*wb*tp*tp/(a*a*b*sinh(wb)),0._dp,-tp/(2*tv2))
               if(.not.(x>0))x=(tp-t)/tv2
               x=min(sqrt(x),wb*(1-1e-6_dp))
            endif
         else
            ! Long way: T ~ sqrt(-b)*(y0+b)*exp(-w/2) as w -> infinity.
            psi=log(tp/t)
            x=sqrt(tp*psi/tv2)
            if(psi>1)x=max(x,2*log(sqrt(-b)*(y0+b)/t))
         endif
      endif
   else
      ! Multiple revolutions: minimum time by Halley on T'(v) = 0, started
      ! at the large-m minimum of y**1.5/sin(v)**3, cos v = b/(y0+b+sqrt(y0*yend)).
      x=acos(b/(y0+b+sqrt(y0*yend)))
      do it=1,20
         call time_derivs(x,1,m,y0,yend,a,b,d)
         p=-2*d(1)*d(2)/(2*d(2)**2-d(1)*d(3))
         p=min(max(x+p,x/2),(x+pi)/2)
         if(abs(p-x)<=1e-12_dp)exit
         x=p
      enddo
      x=p;call time_derivs(x,1,m,y0,yend,a,b,d)
      tmin=d(0)
      status=5
      if(t<tmin*(1-16*epsilon(t)).or..not.(d(2)>0))return
      merged=abs(t-tmin)<=16*epsilon(t)*tmin
      if(merged)then
         status=1
         lo=x;hi=x
      endif
      status=-4
      ! Each branch: cubic Hermite in sig=sqrt(1-(tmin/t)**(1/3)) between the
      ! minimum (v*-v ~ sqrt(6*tmin/T'')*sig) and the pole T ~ K/v**3 at sig=1.
      sig=sqrt(max(1-(tmin/t)**(1._dp/3),0._dp))
      k0=y0*sqrt(y0)*m*pi/sqrt2;kpi=yend*sqrt(yend)*(m+1)*pi/sqrt2
      if(br<0)then
         hi=x
         x=hermite(sig,x,-sqrt(6*tmin/d(2)),0._dp,-2*(k0/tmin)**(1._dp/3))
      else
         lo=x
         x=pi-hermite(sig,pi-x,-sqrt(6*tmin/d(2)),0._dp,-2*(kpi/tmin)**(1._dp/3))
      endif
      x=min(max(x,lo+(hi-lo)*1e-3_dp),hi-(hi-lo)*1e-3_dp)
   endif

   ! Store distance from the upper pole directly on near-full-period arcs.
   ! Subtracting a converged v from pi loses the small arc information.
   if(k==1.and.x>0.75_dp*pi.and..not.merged)then
      p=lo;lo=pi-hi;hi=pi-p;x=pi-x;k=2
   endif

   !======= Finite thin-rectangle ratio and safeguarded Newton =======
   if(x>0.and.yy<0.and..not.merged)then
      if(.not.converge())return
   endif

   !======= Velocities (radial/transverse, cancellation-free) =======
   p=sign(1.0_dp,a)
   if(yy>=0)then
      y=yy*yy;dc=-2*sw-ss*ss/(1+cs)             ! cs-cosh(w), from y
   else
      if(k==2)then
         c=-cos(x);s=sin(x)**2
         y=y0+2*b*cos(x/2)**2
         if(b<0)y=yend-2*b*sin(x/2)**2
      elseif(k>0)then
         c=cos(x);s=sin(x)**2
         if(b>=0)then
            y=y0+2*b*sin(x/2)**2
         else
            y=yend-2*b*cos(x/2)**2
         endif
      else
         c=cosh(x);s=-sinh(x)**2;y=y0-2*b*sinh(x/2)**2
      endif
      dc=p*cs-c
      if(p*cs*c>0.and.abs(p*cs+c)>0.5_dp)dc=(s-ss*ss)/(p*cs+c)
   endif
   ratio=sqrt(rn2/rn1)
   vr1=sqrt2*(dc+p*cs*(rn2-rn1)/(rn1*(ratio+1)))/sqrt(y)
   vr2=-sqrt2*(dc+p*cs*(rn1-rn2)/(rn2*(1/ratio+1)))/sqrt(y)
   vt1=sqrt2*ratio*ss/sqrt(y);vt2=sqrt2*ss/(ratio*sqrt(y))
   v1=(vr1*e1+vt1*t1v)/sqrt(scale)
   v2=(vr2*e2+vt2*t2v)/sqrt(scale)
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

   !> F and its first three derivatives in x: k=+1, x=v (elliptic);
   !! k=-1, x=w (hyperbolic, sinh/cosh).
   pure subroutine f_derivs(x,k,m,f)
   real(dp),intent(in)  :: x
   integer,intent(in)   :: k,m
   real(dp),intent(out) :: f(0:3)
   real(dp) :: s,c,r,g(0:3),xx
   integer  :: j
   if(k==2)then
      s=sin(x);c=cos(x);r=1/s
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
      if(k>0)then
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
   pure subroutine time_derivs(x,k,m,y0,yend,a,b,d)
   real(dp),intent(in)  :: x,y0,yend,a,b
   integer,intent(in)   :: k,m
   real(dp),intent(out) :: d(0:3)
   real(dp) :: sh,ch,s,c,y,y1,y2,y3,f(0:3),q,q1,q2,q3,p,p1,p2,p3
   if(k==2)then
      sh=sin(x/2);ch=cos(x/2);c=(ch-sh)*(ch+sh)
      y=y0+2*b*ch*ch
      if(b<0)y=yend-2*b*sh*sh
   elseif(k>0)then
      sh=sin(x/2);ch=cos(x/2);c=(ch-sh)*(ch+sh)
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
   call f_derivs(x,k,m,f)
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
