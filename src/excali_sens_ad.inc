! AD template for excali_sens.F90, included twice: ORD2 = 0 gives a
! gradient-only type (Jacobian), ORD2 = 1 adds the packed symmetric Hessian
! over the seven variables (r1, r2, alpha). Not compiled on its own.
   integer,parameter :: nv=7
#if ORD2
   integer,parameter :: nh=28
   integer,parameter :: hi(nh)=[1,1,2,1,2,3,1,2,3,4,1,2,3,4,5,1,2,3,4,5,6,1,2,3,4,5,6,7]
   integer,parameter :: hj(nh)=[1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7]
#endif
   type :: ad
      real(dp) :: v=0, g(nv)=0
#if ORD2
      real(dp) :: h(nh)=0
#endif
   end type
   interface operator(+)
      module procedure add_aa, add_ar, add_ra, add_ia
   end interface
   interface operator(-)
      module procedure sub_aa, sub_ar, sub_ra, neg_a, sub_ia
   end interface
   interface operator(*)
      module procedure mul_aa, mul_ar, mul_ra, mul_ia
   end interface
   interface operator(/)
      module procedure div_aa, div_ar, div_ra, div_ia
   end interface
contains
   elemental function var(x,i) result(r)
      real(dp),intent(in)::x
      integer,intent(in)::i
      type(ad)::r
      r%v=x;r%g=0;r%g(i)=1
#if ORD2
      r%h=0
#endif
   end function
   elemental function cst(x) result(r)
      real(dp),intent(in)::x
      type(ad)::r
      r%v=x;r%g=0
#if ORD2
      r%h=0
#endif
   end function
   elemental function add_aa(a,b) result(r)
      type(ad),intent(in)::a,b
      type(ad)::r
      r%v=a%v+b%v;r%g=a%g+b%g
#if ORD2
      r%h=a%h+b%h
#endif
   end function
   elemental function add_ar(a,b) result(r)
      type(ad),intent(in)::a
      real(dp),intent(in)::b
      type(ad)::r
      r=a;r%v=a%v+b
   end function
   elemental function add_ra(a,b) result(r)
      real(dp),intent(in)::a
      type(ad),intent(in)::b
      type(ad)::r
      r=b;r%v=a+b%v
   end function
   elemental function add_ia(a,b) result(r)
      integer,intent(in)::a
      type(ad),intent(in)::b
      type(ad)::r
      r=b;r%v=a+b%v
   end function
   elemental function sub_aa(a,b) result(r)
      type(ad),intent(in)::a,b
      type(ad)::r
      r%v=a%v-b%v;r%g=a%g-b%g
#if ORD2
      r%h=a%h-b%h
#endif
   end function
   elemental function sub_ar(a,b) result(r)
      type(ad),intent(in)::a
      real(dp),intent(in)::b
      type(ad)::r
      r=a;r%v=a%v-b
   end function
   elemental function neg_a(a) result(r)
      type(ad),intent(in)::a
      type(ad)::r
      r%v=-a%v;r%g=-a%g
#if ORD2
      r%h=-a%h
#endif
   end function
   elemental function sub_ra(a,b) result(r)
      real(dp),intent(in)::a
      type(ad),intent(in)::b
      type(ad)::r
      r=neg_a(b);r%v=a-b%v
   end function
   elemental function sub_ia(a,b) result(r)
      integer,intent(in)::a
      type(ad),intent(in)::b
      type(ad)::r
      r=neg_a(b);r%v=a-b%v
   end function
   elemental function mul_aa(a,b) result(r)
      type(ad),intent(in)::a,b
      type(ad)::r
      r%v=a%v*b%v;r%g=a%g*b%v+a%v*b%g
#if ORD2
      r%h=a%h*b%v+a%v*b%h+a%g(hi)*b%g(hj)+b%g(hi)*a%g(hj)
#endif
   end function
   elemental function mul_ar(a,b) result(r)
      type(ad),intent(in)::a
      real(dp),intent(in)::b
      type(ad)::r
      r%v=a%v*b;r%g=a%g*b
#if ORD2
      r%h=a%h*b
#endif
   end function
   elemental function mul_ra(a,b) result(r)
      real(dp),intent(in)::a
      type(ad),intent(in)::b
      type(ad)::r
      r=mul_ar(b,a)
   end function
   elemental function mul_ia(a,b) result(r)
      integer,intent(in)::a
      type(ad),intent(in)::b
      type(ad)::r
      r=mul_ar(b,real(a,dp))
   end function
   ! A function of one AD argument from its value and first two derivatives.
   elemental function lift(x,f0,f1,f2) result(r)
      type(ad),intent(in)::x
      real(dp),intent(in)::f0,f1,f2
      type(ad)::r
      r%v=f0;r%g=f1*x%g
#if ORD2
      r%h=f1*x%h+f2*x%g(hi)*x%g(hj)
#endif
   end function
   elemental function inv(a) result(r)
      type(ad),intent(in)::a
      type(ad)::r
      real(dp)::q
      q=1/a%v
      r=lift(a,q,-q*q,2*q*q*q)
   end function
   elemental function div_aa(a,b) result(r)
      type(ad),intent(in)::a,b
      type(ad)::r
      r=mul_aa(a,inv(b))
   end function
   elemental function div_ar(a,b) result(r)
      type(ad),intent(in)::a
      real(dp),intent(in)::b
      type(ad)::r
      r=mul_ar(a,1/b)
   end function
   elemental function div_ra(a,b) result(r)
      real(dp),intent(in)::a
      type(ad),intent(in)::b
      type(ad)::r
      r=mul_ar(inv(b),a)
   end function
   elemental function div_ia(a,b) result(r)
      integer,intent(in)::a
      type(ad),intent(in)::b
      type(ad)::r
      r=mul_ar(inv(b),real(a,dp))
   end function
   elemental function sqrta(a) result(r)
      type(ad),intent(in)::a
      type(ad)::r
      real(dp)::s
      s=sqrt(a%v)
      r=lift(a,s,.5_dp/s,-.25_dp/(s*a%v))
   end function
   pure function dot3(a,b) result(r)
      type(ad),intent(in)::a(3),b(3)
      type(ad)::r
      r=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
   end function
   pure function cross3(a,b) result(r)
      type(ad),intent(in)::a(3),b(3)
      type(ad)::r(3)
      r(1)=a(2)*b(3)-a(3)*b(2);r(2)=a(3)*b(1)-a(1)*b(3);r(3)=a(1)*b(2)-a(2)*b(1)
   end function
#if ORD2
   ! Packed Hessian over (r1, r2, alpha) -> d2/dz2 over (r1, r2, tof) (tof
   ! row/column zero), d2/dz dalpha, d2/dalpha2.
   pure subroutine unpack(h,zz,za,aa)
      real(dp),intent(in)::h(nh)
      real(dp),intent(out)::zz(7,7),za(7),aa
      integer::k
      zz=0;za=0
      do k=1,nh
         if(hj(k)<7)then
            zz(hi(k),hj(k))=h(k);zz(hj(k),hi(k))=h(k)
         elseif(hi(k)<7)then
            za(hi(k))=h(k)
         else
            aa=h(k)
         endif
      enddo
   end subroutine
#endif

   ! Derivatives at the solver's root: x in coordinate kind (0 X = xx,
   ! 1 v, 2 pi-v, -1 w, -2 xi), yv its y, fg its velocity form, h fixed to
   ! hfix when fixedh.
#if ORD2
   subroutine sens_eval(r1,r2,tof,smu,direction,nrev,kind,x,xx,yv,fg,fixedh,hfix,jac,hess)
      real(dp),intent(out)::hess(6,7,7)
#else
   subroutine sens_eval(r1,r2,tof,smu,direction,nrev,kind,x,xx,yv,fg,fixedh,hfix,jac)
#endif
      real(dp),intent(in)::r1(3),r2(3),tof,smu,x,xx,yv,hfix(3)
      integer,intent(in)::direction,nrev,kind
      logical,intent(in)::fg,fixedh
      real(dp),intent(out)::jac(6,7)
      type(ad)::rr1(3),rr2(3),e1(3),e2(3),nn(3),hh(3),t1(3),t2(3),pv(6), &
         n1,n2,al,sth,cth,cs,ss,s1,s2,pp,qq,a,b,y0,yend,y,c,s,f,gg,xc,sm,u,ig,kk,dc,vr1,vr2,vt1,vt2,wb,ww
      real(dp)::sq,fj(0:2),cj(0:2),gz(7),gza(7),gzz(7,7),az(7),azz(7,7),pz(7),pza(7),pzz(7,7),pa,paa,ga,gaa
      integer::i,j
      do i=1,3
         rr1(i)=var(r1(i),i);rr2(i)=var(r2(i),3+i)
      enddo
      select case(kind)
      case(0); al=var(xx,7)
      case default; al=var(x,7)
      end select
      n1=sqrta(dot3(rr1,rr1));n2=sqrta(dot3(rr2,rr2))
      do i=1,3
         e1(i)=rr1(i)/n1;e2(i)=rr2(i)/n2
      enddo
      nn=cross3(e1,e2)
      if(fixedh)then
         hh=cst(hfix)
         sth=dot3(nn,hh)
      else
         sth=sqrta(dot3(nn,nn))
         do i=1,3
            hh(i)=nn(i)*(direction/sth)
         enddo
         sth=sth*real(direction,dp)
      endif
      cth=dot3(e1,e2)
      if(cth%v>=0)then
         cs=sqrta((1+cth)*.5_dp);if(sth%v<0)cs=-cs
         ss=sth/(2*cs)
      else
         ss=sqrta((1-cth)*.5_dp)
         cs=sth/(2*ss)
      endif
      s1=sqrta(n1);s2=sqrta(n2);pp=s1*s2
      qq=(dot3(rr1-rr2,rr1+rr2)/(n1+n2))/(s1+s2)      ! s1 - s2
      a=sqrt2*pp*cs;b=sqrt2*a
      if(cs%v>0)then
         y0=qq*qq+2*pp*(ss*ss/(1+cs))
      else
         y0=qq*qq+2*pp*(1-cs)
      endif
      if(cs%v<0)then
         yend=qq*qq+2*pp*(ss*ss/(1-cs))
      else
         yend=qq*qq+2*pp*(1+cs)
      endif
      ! y, c = cos(alpha-angle), s = sin**2 (or -sinh**2), F in the coordinate.
      select case(kind)
      case(0)
         call fser(xx,fj,cj)
         xc=lift(al,xx*cj(0),cj(0)+xx*cj(1),2*cj(1)+xx*cj(2))   ! X*C(X) = 1 - cos
         y=y0+b*xc
         c=1-xc;s=xc*(2-xc)
         f=lift(al,fj(0),fj(1),fj(2))
      case(1)
         sq=sin(x)
         if(b%v>=0)then
            y=y0+2*b*lift(al,sin(x/2)**2,.5_dp*sq,.5_dp*cos(x))
         else
            y=yend-2*b*lift(al,cos(x/2)**2,-.5_dp*sq,-.5_dp*cos(x))
         endif
         c=lift(al,cos(x),-sq,-cos(x));s=lift(al,sq*sq,sin(2*x),2*cos(2*x))
         call fell(x,nrev,sq,cos(x),fj)
         f=lift(al,fj(0),fj(1),fj(2))
      case(2)
         sq=sin(x)                                      ! x = d = pi - v
         if(b%v>=0)then
            y=y0+2*b*lift(al,cos(x/2)**2,-.5_dp*sq,-.5_dp*cos(x))
         else
            y=yend-2*b*lift(al,sin(x/2)**2,.5_dp*sq,.5_dp*cos(x))
         endif
         c=lift(al,-cos(x),sq,cos(x));s=lift(al,sq*sq,sin(2*x),2*cos(2*x))
         call fell(pi-x,nrev,sq,-cos(x),fj)
         f=lift(al,fj(0),-fj(1),fj(2))
      case(-2)
         ! xi = wb - w, wb the branch point (y(wb) = 0): y = 2b sinh(wb-xi/2) sinh(xi/2)
         ! keeps full relative precision, and so do its derivatives.
         u=sqrta(y0/(2*b))
         wb=2*lift(u,asinh(u%v),1/sqrt(1+u%v**2),-u%v/(1+u%v**2)**1.5_dp)
         ww=wb-al
         sq=wb%v-.5_dp*x
         y=2*b*lift(wb-.5_dp*al,sinh(sq),cosh(sq),sinh(sq))*lift(al,sinh(.5_dp*x),.5_dp*cosh(.5_dp*x),.25_dp*sinh(.5_dp*x))
         c=lift(ww,cosh(ww%v),sinh(ww%v),cosh(ww%v));s=lift(ww,-sinh(ww%v)**2,-sinh(2*ww%v),-2*cosh(2*ww%v))
         if(ww%v>=.75_dp)then
            call fhyp(ww%v,fj);f=lift(ww,fj(0),fj(1),fj(2))
         else
            xc=-(ww*ww);call fser(xc%v,fj,cj);f=lift(xc,fj(0),fj(1),fj(2))
         endif
      case default
         y=y0-2*b*lift(al,sinh(x/2)**2,.5_dp*sinh(x),.5_dp*cosh(x))
         c=lift(al,cosh(x),sinh(x),cosh(x));s=lift(al,-sinh(x)**2,-sinh(2*x),-2*cosh(2*x))
         call fhyp(x,fj)
         f=lift(al,fj(0),fj(1),fj(2))
      end select
      y%v=yv                                         ! the solver's cancellation-free y
      u=sqrta(y)
      if(kind==-1.and.b%v<0.and.abs(yv*f%v)>4*abs(yv*f%v+a%v))then
         ! Long-way hyperbola where y*F + a cancels (y*F -> -a as w grows):
         ! the solver's form ((y0+b)cosh w - b)/sinh(w)**2 - y w/sinh(w)**3,
         ! all over sqrt2, with y0 + b = q**2 + 2p exactly, does not. (Where
         ! it cancels less than 4x the plain form is the more accurate.)
         sq=sinh(x);c=lift(al,cosh(x),sq,cosh(x))
         gg=u*(((qq*qq+2*pp)*c-b)*lift(al,1/sq**2,-2*cosh(x)/sq**3,-2/sq**2+6*cosh(x)**2/sq**4) &
            -y*al*lift(al,1/sq**3,-3*cosh(x)/sq**4,-3/sq**3+12*cosh(x)**2/sq**5))/sqrt2
         c=lift(al,cosh(x),sinh(x),cosh(x))
      else
         gg=u*(y*f+a)                                   ! G + sqrt(mu)*tof
      endif
      if(fg)then
         ig=1/(a*u)
         sm=y/n1
         do i=1,3
            pv(i)=(rr2(i)-rr1(i)+sm*rr1(i))*ig
         enddo
         sm=y/n2
         do i=1,3
            pv(3+i)=(rr2(i)-rr1(i)-sm*rr2(i))*ig
         enddo
      else
         kk=sqrta(2/y)
         dc=cs-c
         if(cs%v*c%v>0.and.abs(cs%v+c%v)>.5_dp)dc=(s-ss*ss)/(cs+c)
         vr1=kk*(dc-cs*qq/s1);vt1=kk*ss*s2/s1
         vr2=-kk*(dc+cs*qq/s2);vt2=kk*ss*s1/s2
         t1=cross3(hh,e1);t2=cross3(hh,e2)
         do i=1,3
            pv(i)=vr1*e1(i)+vt1*t1(i);pv(3+i)=vr2*e2(i)+vt2*t2(i)
         enddo
      endif
      ! Implicit derivatives of alpha w.r.t. z = [r1, r2, tof].
      gz(1:6)=gg%g(1:6);gz(7)=-1;ga=gg%g(7)              ! (in mu = 1 time units)
      az=-gz/ga;az(7)=az(7)*smu
#if ORD2
      call unpack(gg%h,gzz,gza,gaa)
      do j=1,7
         azz(:,j)=-(gzz(:,j)+gza*az(j)+az*gza(j)+gaa*az*az(j))/ga
      enddo
#endif
      do i=1,6
         pz=0;pz(1:6)=pv(i)%g(1:6);pa=pv(i)%g(7)
         jac(i,:)=(pz+pa*az)*smu
#if ORD2
         call unpack(pv(i)%h,pzz,pza,paa)
         do j=1,7
            hess(i,:,j)=(pzz(:,j)+pza*az(j)+az*pza(j)+paa*az*az(j)+pa*azz(:,j))*smu
         enddo
#endif
      enddo
   end subroutine
