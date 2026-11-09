   ! The solver's stages, included in the contains part of both solve_chunk
   ! (batches) and excali_one (single cases), so each host has its own
   ! internal copy that the compiler can inline; the formulas exist once.

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
         pt=dd(0)/g%t                               ! one quantity for both tests
         stat=5
         if(pt>1+16*epsilon(pt))return
         stat=3
         if(pt>=1-16*epsilon(pt))return
         ut=cbrt_c(pt);sg=br*sqrt((1-pt)/(1+ut+ut*ut))
         k=sqrt(6*dd(0)/dd(2))
         x=vs+k*sg*(1-dd(3)*k*sg/(6*dd(2)))
         stat=1
      end subroutine

      ! Zero-revolution region selector: the coordinate k and the starter x in
      ! it. wb is the short-way branch point; near marks the near-parabolic
      ! window (velocities use Lagrange f,g); fin marks a series that is
      ! already the root. Every starter is within 6.4% of the singularity
      ! bound R, so every rectangle is the same (fast_coordinate).
      !
      !   region                         k      starter
      !   parabola                       0      x = 0
      !   |t-tp| <= 3.1% tp (near)       +-1    series reversion in X
      !     2nd-order term < 1e-8 of X          the series itself (fin)
      !   ellipse, tau >= -0.9           1, 2   table_ve (tabe)
      !   ellipse, tau < -0.9            1, 2   long_d (closed form)
      !   hyperbola, t/tp >= 0.1         -1     table_w (tabh)
      !   hyperbola, t/tp < 0.1          -1,-2  large-w closed form
      !
      ! k=2 is pi-v (used when v > 3pi/4, or d < pi/4 from long_d); k=-2 is
      ! xi=wb-w.
      pure subroutine zero_rev_start(g,k,x,wb,fin,near)
         type(frame),intent(in)::g
         integer,intent(out)::k
         real(dp),intent(out)::x,wb
         logical,intent(out)::fin,near
         real(dp)::tp,u,tau,dc,sl,sw,r0,r1,dt
         tp=sqrt(g%y0)*(g%y0*sqrt2/3+g%a)
         k=0;x=0;wb=0;fin=.false.;near=.false.
         if(abs(g%t-tp)<=4*epsilon(tp)*tp)return    ! parabola
         near=abs(g%t-tp)<=.031_dp*tp      ! the series window (tabh starts at u=0.99)
         k=1
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
            ! to 0: the series is the root.
            call near_parabola_x(g%y0,g%a,g%b,g%t,tp,u,dt)
            if(u*k>0)x=sqrt(abs(u))
            fin=abs(dt)<1e-8_dp
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
            x=log(dc)
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
         if(tau<-.9_dp)then
            ! Long way toward tau -> -1: th=t/|b|**1.5, eta**2=2*yend/|b|; d=pi-v.
            x=long_d(g%t/((-g%b)*sqrt(-g%b)),sqrt(2*g%yend/(-g%b)))
            if(x>.25_dp*pi)then
               x=pi-x
            else
               k=2
            endif
            return
         endif
         ! Elliptic table tabe, split at Th=T(pi/2) (y(pi/2)=y0+b on both
         ! ways): s1=(Th/t)**(1/3) with scale (K/Th)**(1/3), K=pi*y(pi)**1.5/sqrt2
         ! (side 1), or s2=(t-tp)/(Th-tp) with scale s2*(s2+2*tp/(Th-tp)) (side 2).
         sw=g%y0+g%b
         sw=sqrt(sw)*(sw*(pi/(2*sqrt2))+g%a)             ! Th
         if(g%t>=sw)then
            u=rcbrt(sw/g%t);x=table_ve(tau,u,rcbrt(pi*g%ypi*sqrt(g%ypi)/(sqrt2*sw)),1)
         else
            u=(g%t-tp)/(sw-tp);x=table_ve(tau,u,u*(u+2*tp/(sw-tp)),2)
         endif
         if(x>.75_dp*pi)then
            x=pi-x;k=2
         endif
      end subroutine

      ! Coordinate, half-width el and centre half-angle functions of the
      ! rectangle. Every singularity of 1/(T-t) lies at least R from the root:
      ! min(v,pi-v) for ellipses (sqrt(y) branch points, mirror root, pole),
      ! min(w,pi) or min(2w,wb-w,pi) for hyperbolas in w, min(xi,2w,pi) in
      ! xi=wb-w, |X| near the parabola, the bracket [lo,hi] for m>0. Every
      ! rectangle has 8 nodes at L=R/4 (R/8 for m>0, where the other branch's
      ! root can be as close as 2R and t_min's rounding noise is largest),
      ! rounded down to a power of two by kexp (el<=1 keeps the node tables in
      ! range): error about (L/2R)**16 <= 8**-16.
      pure subroutine fast_coordinate(k,x,m,b,wb,lo,hi,kind1,c,el,rb,hs,hc)
         integer,intent(in)::k,m
         real(dp),intent(in)::x,b,wb,lo,hi
         integer,intent(out)::kind1
         real(dp),intent(out)::c,el,rb,hs,hc
         hs=0;hc=1;kind1=k;c=x
         if(m>0)then
            el=min(c-lo,hi-c);rb=el;el=el/8
            call sincos_half(c,hs,hc)
            return
         elseif(abs(k)==1.and.x<.3_dp)then
            kind1=0;c=k*x*x
            el=min(abs(c),.2_dp-abs(c))
            if(k>0)el=min(2*abs(c),.2_dp-abs(c))
            if(k<0.and.b>0)el=min(el,wb*wb-abs(c))
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
         rb=el;el=el/4
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
