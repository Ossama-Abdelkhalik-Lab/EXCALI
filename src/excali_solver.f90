! EXCALI (EXact Cauchy-Analytic Lambert Integration): Lambert's problem for
! every conic and revolution count from the zero-height limit of Exact v2.2's
! thin-rectangle ratio, root = c + Q[((x-c)f)']/Q[f'], f = 1/(T-t), by real
! Gauss-Legendre quadrature. Every case takes one path: a table or closed-form
! starter and one fixed rectangle, its result final, with no search,
! iteration or fallback. Multiple revolutions off the tables (|tau| > 0.9)
! first locate the minimum time by two fixed ratios (mrev_edge).
! Geometry and time functions adapted from ConfluentContour's contour_simple.
module excali_solver
   use iso_fortran_env, only: real32, real64, int64
   use ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
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
   ! Full node lists (nodes, then the endpoints +1,-1) with their weights,
   ! and cos/sin (cosh/sinh) of the half offsets L*u/2 for each L=2**(-k).
   ! They are protected module data, initialized at compile time: indexing a
   ! parameter array at run time makes some gfortran versions copy it to the
   ! stack.
   real(dp),parameter :: u8(n1+2)=[g8x(1:h1),-g8x(1:h1),1._dp,-1._dp],w8(n1+2)=[g8w,g8w,0._dp,0._dp]
   real(dp),protected :: c8f(n1+2,0:kmax)=reshape([((cos(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   real(dp),protected :: s8f(n1+2,0:kmax)=reshape([((sin(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   real(dp),protected :: ch8f(n1+2,0:kmax)=reshape([((cosh(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   real(dp),protected :: sh8f(n1+2,0:kmax)=reshape([((sinh(.5_dp**(k_+1)*u8(j_)),j_=1,n1+2),k_=0,kmax)],[n1+2,kmax+1])
   ! Reciprocal factorials 1/j!, so the Taylor series multiply by constants.
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
      if(status>=0.and.status<=2)then
         v1=v1*vscale;v2=v2*vscale
      endif
   end subroutine

   ! One case, straight-line: the same stages as solve_chunk (shared through
   ! excali_stages.inc) without the block's arrays and stage loops, which only
   ! pay when several cases overlap. Results are identical to excali_batch.
   subroutine excali_one(r1,r2,tof,direction,v1,v2,status,half,nrev,branch,plane_normal)
   real(dp),intent(in)          :: r1(3),r2(3),tof
   integer,intent(in)           :: direction
   real(dp),intent(out)         :: v1(3),v2(3)
   integer,intent(out)          :: status,half
   integer,optional,intent(in)  :: nrev,branch
   real(dp),optional,intent(in) :: plane_normal(3)
   type(frame) :: g
   real(dp) :: x,c,hs,hc,wb,lb,ub,dd,de,rb,el,vs,chord,semi
   integer :: m,br,k,kd,kr
   logical :: fast,near,fin
   m=0;br=-1
   if(present(nrev))m=nrev
   if(present(branch))br=branch
   v1=ieee_value(0.0_dp,ieee_quiet_nan);v2=v1;half=0;status=-1
   if(abs(direction)/=1.or.m<0.or.abs(br)/=1.or..not.(tof>0).or..not.ieee_is_finite(tof))return
   call geometry(r1,r2,tof,direction,g,status,half,plane_normal)
   if(status/=-4)return
   fast=.false.;near=.false.;fin=.false.;rb=0;kr=0;lb=0;ub=pi;hs=0;hc=1;dd=0;de=0;wb=0
   if(m>0)then
      ! Multiple revolutions: the cheap bound, the tables, the fold.
      chord=sqrt((g%rn1-g%rn2)**2+4*g%rn1*g%rn2*g%ss**2);semi=(g%rn1+g%rn2+chord)/2
      status=5
      if(g%t>=m*pi*semi*sqrt(semi)/sqrt2*(1-1e-12_dp)) &
         call mrev_table(g%a/(sqrt2*semi)/(1+sqrt(chord/semi)),semi,g%t,m,br,x,k,vs,status)
      if(status==2)then
         call mrev_fold(g,m,br,vs,x,status);k=1
      endif
      fast=status==1
      if(status==3)then
         status=1;x=vs                             ! t = t_min: the root is v* itself
      elseif(status/=5)then
         status=-4
      endif
      if(status==5)return
      if(br<0)then
         ub=vs
      else
         lb=vs
      endif
      if(k==1.and.x>.75_dp*pi)then
         x=pi-x;k=2
      endif
      if(k==2)then                                 ! distance from the upper pole
         el=lb;lb=pi-ub;ub=pi-el
      endif
      if(status==1)call sincos_half(x,hs,hc)
   else
      call zero_rev_start(g,k,x,wb,fin,near)
      fast=x>0.or.k==0
   endif
   if(fast.and.k/=0)then
      call fast_coordinate(k,x,m,g%b,wb,lb,ub,kd,c,el,rb,hs,hc)
      fast=el>0;kr=kexp(el)
      if(fast.and..not.fin)call rect_ratio(kd,c,hs,hc,kr,g%y0,merge(wb,g%yend,kd==-2),g%a,g%b,g%t,m,dd,de)
      if(fast)call fast_accept(kd,c,dd,hs,hc,k,x,fast)
   endif
   if(fast)then
      status=merge(2,0,de>1e-8_dp*rb*pow2m(kr))     ! status 2: low confidence
      if(k==0)then
         k=1;x=0;hs=0;hc=1                           ! parabola
      endif
   elseif(status==-4.and.m>0)then
      call mrev_edge(g,m,br,k,x,hs,hc,status)
   endif
   if(status>=0.and.status<=2)then
      call velocities(g,r1,r2,direction,k,x,hs,hc,near,v1,v2)
      if(.not.all(ieee_is_finite(v1)).or..not.all(ieee_is_finite(v2)))status=-4
   endif
   if(status<0.or.status>2)then
      v1=ieee_value(0.0_dp,ieee_quiet_nan);v2=v1
   endif

   contains
   include 'excali_stages.inc'
   end subroutine excali_one

   ! Multiple revolutions off the tables (|tau| > 0.9): the same three steps as
   ! the table path, in the coordinate zeta about the crowded pole,
   ! sin(x/2) = s0*sinh(zeta), x = v (short way) or pi-v (long way), where
   ! y = yc*cosh(zeta)**2 exactly and
   !   T = A*coth(zeta)**3/cos(x/2)**3 + sqrt(yc)*cosh(zeta)*(+-yc*cosh(zeta)**2*G(x) + a),
   ! A = n*pi*|b|**1.5/4 (n = m short way, m+1 long way), G the m=0 F. Every
   ! term is evaluated directly, so T and its zeta-derivatives keep full
   ! precision however small v* is.
   ! 1. v*: two 8-node ratios of 1/T' about the large-m estimate (zeta, or x
   !    when that estimate is above 0.2). t_min = T(v*) decides feasibility.
   ! 2. Starter, by the fold distance Delta = (t-t_min)/|t_min-A|: the fold
   !    series; on the crowded side the pole layer A*coth**3 (closed-form
   !    cubics); on the far side the body limit 4*theta*cos(x/2)**3 = n*pi +- (x+sin x),
   !    theta = T/|b|**1.5, inverted from the bt tables (body_x), and the fold
   !    model's quartic by Ferrari (g_root). No step is iterated. Within 15% of R.
   ! 3. One 8-node rectangle at L = R/4, R = the distance to the nearest
   !    singularity (zeta: min(zeta, |zeta-zeta*|, pi/2, zeta_max-zeta); x:
   !    min(|x-x*|, pi-x)), in zeta when the starter has x < 1, in pi-x above 3pi/4.
   ! Returns the root in coordinate k (1: v, 2: pi-v) with its half-angle
   ! functions; status 0, 1 (t = t_min), 2 (low confidence), 5 (t < t_min).
   subroutine mrev_edge(g,m,br,k,x,hs,hc,status)
   type(frame),intent(in) :: g
   integer,intent(in)     :: m,br
   integer,intent(out)    :: k,status
   real(dp),intent(out)   :: x,hs,hc
   real(dp) :: sd,yc,sy,s0,am,b15,t,ve,q,ze,lo,hi,zs,xs,tmin,t2z,t3z,pt,ut,sg,kk,dl,zx,th,d(0:3), &
               err,z0,x0,e2,af,al,be,mm,cc,t2d,mc,dh,kap,db,em,w,x1,rr,el, &
               fr(n1+2),fq(n1+2)
   logical :: crowd
   integer :: i,j,kz                               ! coordinate: 0 zeta, 1 x, 2 pi-x
   t=g%t;sd=sign(1._dp,g%b);k=merge(1,2,g%b>0)
   yc=merge(g%y0,g%yend,g%b>0);sy=sqrt(yc);s0=sqrt(yc/(2*abs(g%b)))
   b15=abs(g%b)*sqrt(abs(g%b));mm=merge(m,m+1,g%b>0)*pi;am=mm*b15/4
   x=0;hs=0;hc=1;status=-4
   ! 1. v* (in x = v or pi-v): the large-m estimate ve, cos ve = b/(y0+b+sqrt(y0*ypi)),
   ! limited by its small-m forms ((3m*pi*y0/b)**(1/3) short way, 2/((m+1)*pi) long way).
   q=g%y0+g%b+sqrt(g%y0*g%ypi)
   if(g%b>0)then
      ve=min(2*asin(sqrt((g%y0+sqrt(g%y0*g%ypi))/(2*q))),cbrt_c(3*m*pi*g%y0/g%b))
   else
      ve=max(2*asin(sqrt((g%ypi+sqrt(g%y0*g%ypi))/(2*q))),2/((m+1)*pi))
   endif
   ze=asinh(sin(ve/2)/s0)
   kz=merge(0,1,ve<.2_dp.and.s0*sinh(ze+2)<.9_dp.and.s0*sinh(2-ze)<.14_dp)
   if(kz==0)then
      lo=ze-2;hi=ze+2
   else
      lo=ve/2;hi=min(2*ve,(ve+pi)/2)
   endif
   ! An 8-node ratio over [lo,hi] (within 4e-5 of R), then one at R/16 about
   ! its result (exact to rounding); R = pi/2 in zeta, min(x,pi-x) in x.
   ze=(lo+hi)/2;el=(hi-lo)/2
   do i=1,2
      do j=1,n1+2
         call tprime(ze+el*u8(j),fr(j),fq(j))
      enddo
      call global_ratio(n1,fr(1:n1),fq(1:n1),fr(n1+2),fr(n1+1),ze-el,ze+el,ze+el*u8(1:n1),w8(1:n1)*el,ze,zs,err)
      ze=zs;el=merge(pi/2,min(zs,pi-zs),kz==0)/16
   enddo
   if(kz==0)then
      xs=2*asin(s0*sinh(zs))
   else
      xs=zs;zs=asinh(sin(xs/2)/s0)
   endif
   if(.not.(zs>0.and.xs>0.and.xs<pi))return
   call tjet(zs,0,d,zx,th)
   tmin=d(0);t2z=d(2);t3z=d(3)
   status=5;pt=tmin/t                              ! one quantity for both tests
   if(pt>1+16*epsilon(pt))return
   if(pt>=1-16*epsilon(pt))then
      status=1;x=xs;call sincos_half(x,hs,hc);return
   endif
   ! 2. Starter: z0 in zeta, or x0 in x (z0 = 0 marks the x form).
   ut=cbrt_c(pt);sg=sd*br*sqrt((1-pt)/(1+ut+ut*ut))  ! sg > 0: zeta > zeta*
   crowd=sg<0;z0=0;x0=0
   dl=(t-tmin)/abs(tmin-am)
   if(dl<merge(.003_dp,.03_dp,g%b>0))then
      kk=sqrt(6*tmin/t2z);z0=zs+kk*sg*(1-t3z*kk*sg/(6*t2z))
   elseif(g%b>0)then
      if(crowd)then
         ! Pole layer coth**3 ~ 1+6e**(-2zeta) fitted at zeta*: 2w**3-D*w**2+1 = 0,
         ! w = e**(zeta-zeta*); below zeta = 2.5 the exact layer a*coth**3 + b*cosh,
         ! fitted at zeta*, with cosh at the layer-only root.
         e2=exp(-2*zs);af=tmin/(1+18*e2)
         z0=zs-log(depressed_max(-(t-af)/(6*af*e2),2._dp))
         if(z0<2.5_dp)then
            q=sinh(zs);rr=tanh(zs)
            al=tmin*rr**3/(1+3/q**2);be=3*al/(rr*rr*q**3)
            z0=atanh(cbrt_c(al/t))
            z0=atanh(cbrt_c(al/(t-be*cosh(z0))))
         endif
      else
         ! The body limit's root vb, final within 8*v* of the far pole;
         ! otherwise its secant slope from v* sets the fold cubic
         ! x + 1/(2x**2) = 1.5 + (t-t_min)/(al*v*), x = v/v*.
         x0=body_x((t-am)/am,1/mm)
         if(pi-x0>=8*xs)then
            q=cos(xs/2);rr=mm+xs+sin(xs)
            if(x0>1.001_dp*xs)then
               al=(t-b15*rr/(4*q**3))/(x0-xs)
            else
               al=b15*(4*q**3+3*rr*sin(xs/2))/(8*q**4)    ! T_body'(v*)
            endif
            dh=1.5_dp+(t-tmin)/(al*xs)
            x0=xs*(depressed_max(-dh*dh/3,.5_dp-2*dh**3/27)+dh/3)
         endif
      endif
   else
      cc=6*am*s0*s0;t2d=t2z*zx*zx;mc=min(6*cc/(xs**4*t2d),1._dp)
      dh=(t-tmin)/(t2d*xs*xs/2)
      if(crowd)then
         if(t>am)z0=atanh(cbrt_c(am/t))
         if(z0>0.and.z0<2.5_dp)then
            continue                                  ! pole layer, exact
         elseif(mc>=.3_dp.or.(mc>=.05_dp.and.dl<1))then
            z0=0;x0=xs*g_root(mc,dh,.false.)
         elseif(t>am)then
            z0=0;x0=cubic_pos(2*(t-am)/b15,2*cc/b15)
         else
            ! Body anchored at the minimum time, linear to the plateau: the
            ! body root x1 by the quadratic through (0,am) and its minimum
            ! (x_b,am*(1+em)), w = (am-t+tmin-tb)/(am-tb), tb = am*(1+em),
            ! x1 = x_b*w/(1+sqrt(1-w)); slope (am-t+tmin-tb)/x1.
            call body_min(1/mm,db,em)
            w=min(max((t-tmin)/(am*em)+1,0._dp),1._dp)
            al=-am*em*(1+sqrt(1-w))/db
            x0=cubic_pos((t-am)/al,cc/al)
         endif
      else
         if(mc>=.3_dp)x0=xs*g_root(mc,dh,.true.)
         if(.not.(x0>0.and.x0<1.5_dp*xs))then
            ! Body plus the layer term cc/d**2, offset to t_min at d*: the body
            ! roots without the layer (x1) and with its value at d* (x0)
            ! bracket the root; one linear interpolation of the residual.
            kap=tmin-b15*(mm-xs-sin(xs))/(4*cos(xs/2)**3)-cc/xs**2
            x1=body_x((t-kap-am)/am,-1/mm);x0=body_x((t-kap-am-cc/xs**2)/am,-1/mm)
            rr=cc/x1**2;el=cc/x0**2-cc/xs**2
            if(rr>el)x0=x0-(x1-x0)*el/(rr-el)
         endif
      endif
   endif
   if(z0>0)x0=2*asin(min(s0*sinh(z0),1._dp))
   kz=1
   if(x0<1)then
      kz=0;if(.not.(z0>0))z0=asinh(sin(x0/2)/s0)
   elseif(x0>.75_dp*pi)then
      kz=2
   endif
   ! 3. One rectangle at L = R/4.
   select case(kz)
    case(0)
      ze=z0;rr=min(z0,abs(z0-zs),pi/2,asinh(1/s0)-z0)
    case(1)
      ze=x0;rr=min(abs(x0-xs),pi-x0)
    case default
      ze=pi-x0;rr=min(abs(x0-xs),ze)
   end select
   el=rr/4;status=-4
   if(.not.(el>0))return
   do j=1,n1+2
      call tval(ze+el*u8(j),fr(j),fq(j))
   enddo
   call global_ratio(n1,fr(1:n1),fq(1:n1),fr(n1+2),fr(n1+1),ze-el,ze+el,ze+el*u8(1:n1),w8(1:n1)*el,ze,x,err)
   status=merge(2,0,err>5e-9_dp*rr)
   if(kz==0)then
      hs=s0*sinh(x);hc=sqrt((1-hs)*(1+hs));x=2*asin(hs)
   else
      call sincos_half(x,hs,hc)
      if(kz==2)k=3-k
   endif
   if(.not.(x>0.and.x<pi))status=-4

   contains

      ! T and its first three zeta-derivatives at z in coordinate kc (0 zeta,
      ! 1 x, 2 pi-x); zx = dzeta/dx, th = tanh(zeta).
      subroutine tjet(z,kc,d,zx,th)
         real(dp),intent(in)::z
         integer,intent(in)::kc
         real(dp),intent(out)::d(0:3),zx,th
         real(dp)::sh,ch,u,cu,xv,r,w,cth(0:3),p(0:3),uj(0:3),gx(0:3),cj(0:3)
         select case(kc)
          case(0)
            call shch(z,sh,ch);u=s0*sh;cu=sqrt((1-u)*(1+u));xv=2*asin(u)
          case(1)
            xv=z;call sincos_half(xv,u,cu)
          case default
            xv=pi-z;call sincos_half(z,cu,u)
         end select
         if(kc>0)then
            sh=u/s0;ch=sqrt(1+sh*sh)
         endif
         r=1/(sh*sh);w=ch/sh
         cth=[w,-r,2*w*r,-2*r*(r+2*w*w)]               ! coth
         uj=s0*[sh,ch,sh,ch]                            ! sin(x/2)
         w=cu*cu
         ! A*coth**3*(1-u**2)**(-3/2)
         p=jmul(jmul(jmul(cth,cth),cth),jcomp([1/(w*cu),-1.5_dp/(w*w*cu),3.75_dp/(w**3*cu),-13.125_dp/(w**4*cu)], &
            [w,-2*uj(0)*uj(1),-2*(uj(1)**2+uj(0)*uj(2)),-2*(3*uj(1)*uj(2)+uj(0)*uj(3))]))
         call f_derivs(xv,1,0,gx,2*u*cu,(cu-u)*(cu+u))  ! G and its x-derivatives
         gx=jcomp(gx,jcomp([xv,2/cu,2*u/cu**3,2*(1+2*u*u)/cu**5],uj))
         cj=[ch,sh,ch,sh]
         d=am*p+sy*(sd*yc*jmul(jmul(jmul(cj,cj),cj),gx)+g%a*cj)
         zx=cu/(2*s0*ch);th=sh/ch
      end subroutine

      ! f = 1/T' and f' in the v* coordinate (zeta, or x with the chain rule).
      subroutine tprime(z,f,f1)
         real(dp),intent(in)::z
         real(dp),intent(out)::f,f1
         real(dp)::d(0:3),zx,th,d1,d2
         call tjet(z,kz,d,zx,th)
         d1=d(1);d2=d(2)
         if(kz==1)then
            d2=d(2)*zx*zx-d(1)*th*(.25_dp+zx*zx);d1=d(1)*zx
         endif
         f=1/d1;f1=-d2*f*f
      end subroutine

      ! f = 1/(T-t) and f' in the rectangle's coordinate.
      subroutine tval(z,f,f1)
         real(dp),intent(in)::z
         real(dp),intent(out)::f,f1
         real(dp)::d(0:3),zx,th
         call tjet(z,kz,d,zx,th)
         if(kz>0)d(1)=d(1)*zx*(3-2*kz)
         f=1/(d(0)-t);f1=-d(1)*f*f
      end subroutine

   end subroutine mrev_edge

   ! Derivatives (orders 0-3) of a product, and of h(u(z)) from the outer
   ! derivatives h at u(z) (Faa di Bruno).
   pure function jmul(a,b) result(c)
      real(dp),intent(in)::a(0:3),b(0:3)
      real(dp)::c(0:3)
      c=[a(0)*b(0),a(1)*b(0)+a(0)*b(1),a(2)*b(0)+2*a(1)*b(1)+a(0)*b(2),a(3)*b(0)+3*(a(2)*b(1)+a(1)*b(2))+a(0)*b(3)]
   end function

   pure function jcomp(h,u) result(c)
      real(dp),intent(in)::h(0:3),u(0:3)
      real(dp)::c(0:3)
      c=[h(0),h(1)*u(1),h(2)*u(1)**2+h(1)*u(2),h(3)*u(1)**3+3*h(2)*u(1)*u(2)+h(1)*u(3)]
   end function

   ! Body limit Q(x) = (1+mu*(x+sin x))/cos(x/2)**3 = 1+e, x = v (short way,
   ! mu = 1/M) or pi-v (long way, mu = -1/M, side x > x_b of the minimum), from
   ! the bt tables (gen_tables.py): x/xq or (x-x_b)/(xq-x_b) below pi/2,
   ! cos(x/2)/c0, c0 = Q**(-1/3), above; within 3e-4 of min(x-x_b, pi-x).
   pure real(dp) function body_x(e,mu) result(x)
      real(dp),intent(in)::e,mu
      real(dp),parameter::mu_s=1/pi,mu_l=1/(2*pi)
      real(dp)::qh,a,db,em,q,qn,r,wt(4),ws(4),u
      integer::i,j,h,k
      qh=2*sqrt2*(1+mu*(pi/2+1))
      if(1+e>qh)then
         q=1/cbrt_c(1+e);a=(mu+mu_l)*((bt_nm-1)/(mu_s+mu_l));u=q*cbrt_c(qh)*(bt_na-1)
         call stencil(a,u,i,j,ws,wt)
         r=0
         do k=1,4
            r=r+ws(k)*(wt(1)*bt_b(i+1,j+k)+wt(2)*bt_b(i+2,j+k)+wt(3)*bt_b(i+3,j+k)+wt(4)*bt_b(i+4,j+k))
         enddo
         x=pi-2*asin(min(q*r,1._dp));return
      endif
      if(mu>=0)then
         q=2*e/(2*mu+sqrt(4*mu*mu+1.5_dp*e));qn=2*(qh-1)/(2*mu+sqrt(4*mu*mu+1.5_dp*(qh-1)))
         a=mu/mu_s;h=1;db=0
      else
         a=-mu;call body_min(a,db,em);a=a/mu_l;h=2
         q=sqrt(max(e-em,0._dp)*(8._dp/3));qn=sqrt((qh-1-em)*(8._dp/3))
      endif
      call stencil(a*(bt_nm-1),q/qn*(bt_na-1),i,j,ws,wt)
      r=0
      do k=1,4
         r=r+ws(k)*(wt(1)*bt_a(i+1,j+k,h)+wt(2)*bt_a(i+2,j+k,h)+wt(3)*bt_a(i+3,j+k,h)+wt(4)*bt_a(i+4,j+k,h))
      enddo
      x=db+q*r
   contains
      pure subroutine stencil(a,u,i,j,ws,wt)
         real(dp),intent(in)::a,u
         integer,intent(out)::i,j
         real(dp),intent(out)::ws(4),wt(4)
         j=min(max(int(a)-1,0),bt_nm-4);call lagrange4(a-j,ws)
         i=min(max(int(u)-1,0),bt_na-4);call lagrange4(u-i,wt)
      end subroutine
   end function

   ! Long-way body minimum for a = 1/M: x_b = (8a/3)(1 + a**2 poly(a**2)) (2e-8)
   ! and em = Q(x_b) - 1 < 0, cancellation-free (stationary, so exact to 1e-16).
   pure subroutine body_min(a,db,em)
      real(dp),intent(in)::a
      real(dp),intent(out)::db,em
      real(dp)::c,a2
      a2=a*a
      db=8*a/3*(1+a2*((((bt_db(1)*a2+bt_db(2))*a2+bt_db(3))*a2+bt_db(4))*a2+bt_db(5)))
      c=cos(db/2);em=(2*sin(db/4)**2*(1+c+c*c)-a*(db+sin(db)))/c**3
   end subroutine

   ! Roots of the fold model g(x) = (x-1)**2*(a2+2c/x+c/x**2) = dh (c = mc/3,
   ! a2 = 1-mc): in w = 1/x the quartic c*w**4+(a2-3c-dh)*w**2+2(c-a2)*w+a2 = 0
   ! is depressed, and Ferrari's resolvent puts both roots on w**2-s*w+k = 0:
   ! left root 1/w+, right root w+/k.
   pure real(dp) function g_root(mc,dh,right) result(x)
      real(dp),intent(in)::mc,dh
      logical,intent(in)::right
      real(dp)::c,a2,p,q,r,bb,y,s,k,w
      c=mc/3;a2=1-mc;p=(a2-3*c-dh)/c;q=2*(c-a2)/c;r=a2/c
      bb=(p*p-4*r)/4
      y=depressed_max(bb-p*p/3,2*p**3/27-p*bb/3-q*q/8)-p/3
      s=sqrt(2*y);k=p/2+y+q/(2*s)
      w=s/2+sqrt(max(s*s/4-k,0._dp))
      x=merge(w/k,1/w,right)
   end function

   ! Positive root of d**3 + e*d**2 = f (f > 0), cancellation-free: with
   ! d = f**(1/3)*X, X**3 + chi*X**2 = 1; for chi > 0 X = 1/W, W the largest
   ! root of W**3 - chi*W - 1 = 0, else X - chi/3 the largest root of the
   ! depressed cubic.
   pure real(dp) function cubic_pos(e,f) result(d)
      real(dp),intent(in)::e,f
      real(dp)::s,chi
      s=cbrt_c(f);chi=e/s
      if(chi>0)then
         d=s/depressed_max(-chi,-1._dp)
      elseif(chi<0)then
         d=s*(depressed_max(-chi*chi/3,2*chi**3/27-1)-chi/3)
      else
         d=s
      endif
   end function


   ! Batch of normalized solves, same arguments per case as excali_one, in
   ! chunks of nb cases through solve_chunk (the same stages as excali_one;
   ! results agree to rounding and statuses are identical). A batch of one
   ! case runs excali_one itself.
   subroutine excali_batch(n,r1,r2,tof,direction,v1,v2,status,half,nrev,branch)
   integer,intent(in)   :: n
   real(dp),intent(in)  :: r1(3,n),r2(3,n),tof(n)
   integer,intent(in)   :: direction(n),nrev(n),branch(n)
   real(dp),intent(out) :: v1(3,n),v2(3,n)
   integer,intent(out)  :: status(n),half(n)
   integer :: i0
   ! Chunks are independent; under OpenMP they run in parallel (dynamic
   ! schedule). Small batches stay serial outside any OpenMP region: entering
   ! one costs about 200-300 ns per call even when an if clause keeps it on
   ! one thread.
   if(n==1)then
      ! One case: the straight-line path (no block arrays), 7% faster.
      call excali_one(r1(:,1),r2(:,1),tof(1),direction(1),v1(:,1),v2(:,1),status(1),half(1),nrev(1),branch(1))
   elseif(n<=8*nb)then
      do i0=1,n,nb
         call one_chunk(i0)
      enddo
   else
      !$omp parallel do schedule(dynamic)
      do i0=1,n,nb
         call one_chunk(i0)
      enddo
      !$omp end parallel do
   endif

   contains
   subroutine one_chunk(i0)
      integer,intent(in)::i0
      integer::nc
      nc=min(nb,n-i0+1)
      call solve_chunk(nc,r1(:,i0:i0+nc-1),r2(:,i0:i0+nc-1),tof(i0:i0+nc-1),direction(i0:i0+nc-1), &
         v1(:,i0:i0+nc-1),v2(:,i0:i0+nc-1),status(i0:i0+nc-1),half(i0:i0+nc-1),nrev(i0:i0+nc-1),branch(i0:i0+nc-1))
   end subroutine
   end subroutine excali_batch

   ! nc <= nb normalized solves, stage by stage so independent cases overlap:
   ! geometry; table or closed-form starters; one 8-node rectangle at a fixed
   ! fraction of a lower bound on the distance to every singularity;
   ! mrev_edge for multiple revolutions off the tables; velocities. Every
   ! case takes exactly one path.
   subroutine solve_chunk(nc,r1c,r2c,tofc,dirc,v1c,v2c,statc,halfc,nrev,branch,pnormal)
   integer,intent(in)   :: nc
   real(dp),intent(in)  :: r1c(3,nc),r2c(3,nc),tofc(nc)
   integer,intent(in)   :: dirc(nc),nrev(nc),branch(nc)
   real(dp),intent(out) :: v1c(3,nc),v2c(3,nc)
   integer,intent(out)  :: statc(nc),halfc(nc)
   real(dp),optional,intent(in) :: pnormal(3)
   type(frame) :: gc(nb)
   real(dp) :: xc(nb),cc(nb),hsc(nb),hcc(nb),wbc(nb),lbc(nb),ubc(nb),ddc(nb),dec(nb),rbc(nb),elc,vsc,chord,semi
   integer :: kc(nb),kd(nb),krc(nb),j
   logical :: fast(nb),nearc(nb),finc(nb)
   ! Stage 1: geometry; the multi-revolution tables.
   do j=1,nc
      fast(j)=.false.;nearc(j)=.false.;halfc(j)=0;statc(j)=-1;rbc(j)=0;krc(j)=0
      if(abs(dirc(j))/=1.or.nrev(j)<0.or.abs(branch(j))/=1.or..not.(tofc(j)>0) &
         .or..not.ieee_is_finite(tofc(j)))cycle
      call geometry(r1c(:,j),r2c(:,j),tofc(j),dirc(j),gc(j),statc(j),halfc(j),pnormal)
      if(statc(j)/=-4)cycle
      lbc(j)=0;ubc(j)=pi;finc(j)=.false.
      if(nrev(j)>0)then
         ! Every ellipse through both points has a>=s/2, so T > m*2*pi*(s/2)**1.5
         ! (within 31% of the minimum time) rejects most infeasible times.
         chord=sqrt((gc(j)%rn1-gc(j)%rn2)**2+4*gc(j)%rn1*gc(j)%rn2*gc(j)%ss**2);semi=(gc(j)%rn1+gc(j)%rn2+chord)/2
         statc(j)=5
         if(gc(j)%t>=nrev(j)*pi*semi*sqrt(semi)/sqrt2*(1-1e-12_dp)) &
            call mrev_table(gc(j)%a/(sqrt2*semi)/(1+sqrt(chord/semi)),semi,gc(j)%t,nrev(j),branch(j),xc(j),kc(j),vsc,statc(j))
         if(statc(j)==2)then
            call mrev_fold(gc(j),nrev(j),branch(j),vsc,xc(j),statc(j));kc(j)=1
         endif
         fast(j)=statc(j)==1
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
         if(kc(j)==1.and.xc(j)>.75_dp*pi)then
            xc(j)=pi-xc(j);kc(j)=2
         endif
         if(kc(j)==2)then                           ! distance from the upper pole
            elc=lbc(j);lbc(j)=pi-ubc(j);ubc(j)=pi-elc
         endif
         if(statc(j)==1)call sincos_half(xc(j),hsc(j),hcc(j))
      endif
   enddo
   ! Stage 2: the zero-revolution region selector (its own loop: inside the
   ! geometry loop it costs 5-10 ns per solve in code generation), then the
   ! rectangle's coordinate, centre and half-width.
   do j=1,nc
      if(statc(j)/=-4.or.nrev(j)>0)cycle
      call zero_rev_start(gc(j),kc(j),xc(j),wbc(j),finc(j),nearc(j))
      fast(j)=xc(j)>0.or.kc(j)==0
   enddo
   do j=1,nc
      if(.not.fast(j).or.kc(j)==0)cycle
      call fast_coordinate(kc(j),xc(j),nrev(j),gc(j)%b,wbc(j),lbc(j),ubc(j),kd(j),cc(j),elc,rbc(j),hsc(j),hcc(j))
      fast(j)=elc>0;krc(j)=kexp(elc)
   enddo
   ! Stage 3: one rectangle each, then its acceptance.
   do j=1,nc
      ddc(j)=0;dec(j)=0
      if(fast(j).and.kc(j)/=0.and..not.finc(j))call rect_ratio(kd(j),cc(j),hsc(j),hcc(j),krc(j),gc(j)%y0,merge(wbc(j),gc(j)%yend,kd(j)==-2),gc(j)%a,gc(j)%b, &
         gc(j)%t,nrev(j),ddc(j),dec(j))
   enddo
   do j=1,nc
      if(fast(j).and.kc(j)/=0)call fast_accept(kd(j),cc(j),ddc(j),hsc(j),hcc(j),kc(j),xc(j),fast(j))
   enddo
   ! Stage 4: multiple revolutions off the tables (mrev_edge), then velocities.
   do j=1,nc
      if(fast(j))then
         ! Status 2: solved, low confidence. The error estimate de (relative to
         ! the half-width L) exceeds 1e-8 of the singularity bound R: a starter
         ! or width regression. Never observed; the tests assert it.
         statc(j)=merge(2,0,dec(j)>1e-8_dp*rbc(j)*pow2m(krc(j)))
         if(kc(j)==0)then
            kc(j)=1;xc(j)=0;hsc(j)=0;hcc(j)=1              ! parabola
         endif
      elseif(statc(j)==-4.and.nrev(j)>0)then
         call mrev_edge(gc(j),nrev(j),branch(j),kc(j),xc(j),hsc(j),hcc(j),statc(j))
      endif
      if(statc(j)>=0.and.statc(j)<=2)then
         call velocities(gc(j),r1c(:,j),r2c(:,j),dirc(j),kc(j),xc(j),hsc(j),hcc(j),nearc(j),v1c(:,j),v2c(:,j))
         if(.not.all(ieee_is_finite(v1c(:,j))).or..not.all(ieee_is_finite(v2c(:,j))))statc(j)=-4
      endif
      if(statc(j)<0.or.statc(j)>2)then
         v1c(:,j)=ieee_value(0.0_dp,ieee_quiet_nan);v2c(:,j)=v1c(:,j)
      endif
   enddo

   ! The stages' steps (internal, so the compiler inlines them).
   contains
   include 'excali_stages.inc'
   end subroutine solve_chunk

   ! Zero-height thin rectangle on [c-L,c+L], L=2**(-kr): Gauss-Legendre
   ! 8 nodes plus the two endpoints. dd=Q1/Q0 locates the
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

   pure subroutine rect_ratio(kind1,c,hs,hc,kr,y0,yend,a,b,t,m,dd,de)
      integer,intent(in)::kind1,kr,m
      real(dp),intent(in)::c,hs,hc,y0,yend,a,b,t
      real(dp),intent(out)::dd,de
      real(dp),dimension(n1+2)::uu,wq,ct,st,f,f1
      real(dp)::el,q0,q1,q2,r0,r1,r2,fa,fb,x,sh,ch,sn,cn,rs,yv,y1,ff,ff1,q,g,g1,ser,ser1,pk,y00,ab,mt,sb,cb
      integer::np,j,n
      logical::lw,us
      el=pow2m(kr);np=n1+2;uu=u8;wq=w8
      if(kr>kmax)then
         ! Below 2**-62 (xi ~ t**2 at the shortest times): cos=cosh=1 and
         ! sin=sinh=L*u/2 to rounding.
         ct=1;st=el*uu/2
      elseif(kind1<0)then
         ct=ch8f(:,kr);st=sh8f(:,kr)
      elseif(kind1>0)then
         ct=c8f(:,kr);st=s8f(:,kr)
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
      ! Q[g'] = g(c+L)-g(c-L) - L*sum w g'(x) for g = f, (x-c)f, (x-c)**2 f;
      ! dd=Q1/Q0 is the root's offset. de=|Q2/Q0-dd**2| vanishes for a single
      ! pole: it monitors the result (status 2 when large), it never re-solves.
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
      q0=fb-fa-el*(q0+r0)
      q1=el*(fb+fa)-el*(q1+r1)
      q2=el*el*(fb-fa)-el*(q2+r2)
      dd=q1/q0
      de=abs(q2/q0-dd*dd)
      if(.not.(dd==dd))dd=huge(dd)
      if(.not.(de==de))de=huge(de)
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

   ! Zero-revolution elliptic root v from tabe at tau in [-0.9,1-1e-6] (grid
   ! uniform in sign(tau)*(1-sqrt(1-|tau|)), clamped beyond 1-1e-6) and s in
   ! [0,1], by local 4x4 cubic Lagrange interpolation: side 1 (t>=T(pi/2))
   ! v=pi-s*sc*A, side 2 v=sqrt(B*sc) (sc from zero_rev_start). 0 outside.
   pure real(dp) function table_ve(tau,s,sc,side) result(v)
      real(dp),intent(in)::tau,s,sc
      integer,intent(in)::side
      real(dp)::u,wt(4),ws(4),r(4)
      integer::i,j
      v=0
      u=min(sign(1-sqrt(1-abs(tau)),tau),tabe_z1)
      if(.not.(u>=tabe_z0.and.s>=0.and.s<=1))return
      u=(u-tabe_z0)*((tabe_nt-1)/(tabe_z1-tabe_z0))
      i=min(max(int(u)-1,0),tabe_nt-4);call lagrange4(u-i,wt)
      u=s*(tabe_ns-1)
      j=min(max(int(u)-1,0),tabe_ns-4);call lagrange4(u-j,ws)
      ! Four independent column sums, then one combination (short chains).
      r=wt(1)*tabe(i+1,j+1:j+4,side)+wt(2)*tabe(i+2,j+1:j+4,side)+wt(3)*tabe(i+3,j+1:j+4,side)+wt(4)*tabe(i+4,j+1:j+4,side)
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
   ! linearized layer) when that gives d<0.1; else the body limit, a
   ! polynomial fit on each half of u (bb, within 2.4e-5).
   pure real(dp) function long_d(th,eta) result(d)
      real(dp),intent(in)::th,eta
      real(dp)::dl,x,u,h,z
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
      d=cubic_pos(2*dl,.75_dp*pi*eta*eta)
      if(d<.1_dp)return
      u=min(max((th-1._dp/3)/(pi/4-1._dp/3),0._dp),1._dp)
      if(u<=.5_dp)then
         z=4*u-1;h=bb(bb_n,1)
         do j=bb_n-1,1,-1
            h=h*z+bb(j,1)
         enddo
         d=pi-sqrt(max(h*u,0._dp))
      else
         z=4*u-3;h=bb(bb_n,2)
         do j=bb_n-1,1,-1
            h=h*z+bb(j,2)
         enddo
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
   pure subroutine mrev_table(tau,sl,t,m,br,v,kx,vs,stat)
      real(dp),intent(in)::tau,sl,t
      integer,intent(in)::m,br
      real(dp),intent(out)::v,vs
      integer,intent(out)::stat,kx
      real(dp),parameter::margin=1e-6_dp
      real(dp)::u,wt(4),wc(4),ws(4),wm(4),mu,cm,r,sig,dpl
      integer::i,ic,j,k,km(4),nk,ia,jb,kk,h
      stat=0;v=0;vs=0;kx=1
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
      u=(asin(tau*mr_tk)/mr_tc+1)*((mr_nt-1)/2._dp)   ! grid clustered toward |tau| = 0.9
      i=min(max(int(u)-1,0),mr_nt-4);call lagrange4(u-i,wt)
      u=sqrt(sig)*(mr_ns-1)                            ! grid uniform in sqrt|sigma|
      j=min(max(int(u)-1,0),mr_ns-4);call lagrange4(u-j,ws)
      do k=1,nk
         do jb=1,4
            v=v+wm(k)*ws(jb)*(wt(1)*mr_q(i+1,j+jb,h,km(k))+wt(2)*mr_q(i+2,j+jb,h,km(k)) &
               +wt(3)*mr_q(i+3,j+jb,h,km(k))+wt(4)*mr_q(i+4,j+jb,h,km(k)))
         enddo
      enddo
      ! The table holds g=(q-1)/(1-sig); the pole distance
      ! (1-sig)*(1-sig*g)*|P-v*|, 1-sig = r/(1+sig), keeps full precision as
      ! t -> infinity. The root comes back in v (kx=1), or in pi-v near pi (kx=2).
      dpl=r/(1+sig)*(1-sig*v)*merge(pi-vs,vs,br>0)
      kx=1
      if(br<0)then
         v=dpl
      elseif(dpl<.25_dp*pi)then
         v=dpl;kx=2
      else
         v=pi-dpl
      endif
      if(v>0.and.dpl>0.and.dpl<merge(pi-vs,vs,br>0).and.vs>0.and.vs<pi)stat=1
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

   ! Lagrange f,g velocities from the exact inputs: v1 = (r2 - f*r1)/g and
   ! v2 = (gd*r2 - r1)/g with f = 1-y/|r1|, gd = 1-y/|r2|, g = A*sqrt(y),
   ! A**2 = |r1||r2| + r1.r2, written as (r2-r1 + (y/|r1|)r1)/g and
   ! (r2-r1 - (y/|r2|)r2)/g. In plain double these are within 5 ulp of the
   ! double-double evaluation (20 at radius ratios near 150), whose extra
   ! precision is below the error y itself carries: errors against reference
   ! velocities are the same at p99, p99.9 and max, at 30 ns less per case.
   pure subroutine fg_velocities(r1,r2,yn,scale,d,v1,v2)
      real(dp),intent(in)::r1(3),r2(3),yn,scale
      integer,intent(in)::d
      real(dp),intent(out)::v1(3),v2(3)
      real(dp)::y,n1,n2,ig,f1,f2,s
      integer::k
      y=yn*scale
      n1=sqrt(dot_product(r1,r1));n2=sqrt(dot_product(r2,r2))
      ig=1/(d*sqrt(n1*n2+dot_product(r1,r2))*sqrt(y))
      f1=y/n1;f2=y/n2
      do k=1,3
         s=r2(k)-r1(k)
         v1(k)=(s+f1*r1(k))*ig;v2(k)=(s-f2*r2(k))*ig
      enddo
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
