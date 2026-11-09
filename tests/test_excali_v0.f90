program test_excali_v0
   use iso_fortran_env, only: real64
   use ieee_arithmetic
   use excali_v0_adapter
   implicit none
   integer,parameter::dp=real64,n=517
   real(dp)::r1(3,n),r2(3,n),t(n),v1(3,n),v2(3,n),w1(3),w2(3),a,th,inc,u(3),v(3),err
   integer::d(n),m(n),b(n),s(n),h(n),i,k,st,hf
   call excali_v0_initialize()
   u=[0.6_dp,0.8_dp,0._dp];v=[-0.48_dp,0.36_dp,0.8_dp]
   do i=1,n
      a=0.7_dp+3._dp*i/n;th=0.1_dp+6.0_dp*i/n
      d(i)=merge(1,-1,th<acos(-1._dp));m(i)=0;b(i)=-1
      r1(:,i)=a*u;r2(:,i)=a*(cos(th)*u+sin(th)*v);t(i)=th*a**1.5_dp
   enddo
   do k=0,0
      call excali_v0_batch(n,r1,r2,t,d,m,b,v1,v2,s,h)
      do i=1,n
         a=norm2(r1(:,i));th=t(i)/a**1.5_dp
         w1=v/sqrt(a);w2=(-sin(th)*u+cos(th)*v)/sqrt(a)
         err=max(norm2(v1(:,i)-w1),norm2(v2(:,i)-w2))
         ! v0's accuracy: 2.8e-6 worst on these orbits (its README: 1% of cases above 1.4e-5).
         if(s(i)/=0.or..not.ieee_is_finite(err).or.err>1e-5_dp)error stop 'v0 circular truth'
         call excali_v0_one(r1(:,i),r2(:,i),t(i),d(i),0,-1,w1,w2,st,hf)
         err=max(norm2(v1(:,i)-w1),norm2(v2(:,i)-w2))
         if(st/=0.or..not.ieee_is_finite(err).or.err>1e-10_dp)error stop 'scalar/batch agreement'
      enddo
      call excali_v0_one(r1(:,1),r2(:,1),t(1),d(1),1,-1,w1,w2,st,hf)
      if(st/=-6.or.any(ieee_is_finite(w1)))error stop 'unsupported multirev'
      call excali_v0_one(r1(:,1),r2(:,1),t(1)*0.001_dp,d(1),0,-1,w1,w2,st,hf)
      if(st>=0)error stop 'unsupported hyperbola'
   enddo
   print *, 'PASS: EXCALI v0, 517 rotated circular truths, both directions, scalar/batch, rejected regimes'
end program
