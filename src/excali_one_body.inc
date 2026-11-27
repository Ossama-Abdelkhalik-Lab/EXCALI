   ! The body of excali_one after geometry: one (m, branch) solve in the frame g.
   ! Included by excali_one and, per revolution count, by excali_solve_all.
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
      call fast_coordinate(k,x,m,g%y0,g%b,wb,lb,ub,kd,c,el,rb,hs,hc)
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

