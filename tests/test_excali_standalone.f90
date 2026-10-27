! Analytic circular, parabolic and hyperbolic truth, extreme unit scales,
! both directions, multirevolution and exact antipodal transfer-plane API.
program test_excali_standalone
   use iso_fortran_env, only: real64
   use ieee_arithmetic
   use excali_solver
   implicit none
   integer,parameter::dp=real64
   real(dp),parameter::pi=acos(-1._dp)
   real(dp)::r1(3),r2(3),w1(3),w2(3),v1(3),v2(3),th,t,scale,mu,tu,vu,e,a,hh,nu,rr,mm
   integer::i,m,d,b,s,h,j,found
   do m=0,10
      do i=1,50
         th=.03_dp+(2*pi-.06_dp)*i/51
         d=merge(1,-1,th<pi)
         r1=[1._dp,0._dp,0._dp];r2=[cos(th),sin(th),0._dp];t=2*pi*m+th
         w1=[0._dp,1._dp,0._dp];w2=[-sin(th),cos(th),0._dp]
         found=0
         do b=-1,1,2
            if(m==0.and.b==1)cycle
            call excali_one(r1,r2,t,d,v1,v2,s,h,m,b)
            if(s/=0.and.s/=1)error stop 'circular status'
            if(max(norm2(v1-w1),norm2(v2-w2))<1e-10_dp)found=found+1
         enddo
         if(found<1)error stop 'circular truth branch missing'
      enddo
   enddo
   do j=1,3
      scale=1;mu=1
      if(j==2)then
         scale=1e150_dp;mu=1e300_dp
      elseif(j==3)then
         scale=1e-150_dp;mu=1e-300_dp
      endif
      tu=exp(1.5_dp*log(scale)-.5_dp*log(mu));vu=exp(.5_dp*(log(mu)-log(scale)))
      call excali_solve([scale,0._dp,0._dp],[-scale,0._dp,0._dp],pi*tu,mu,1,0,-1,v1,v2,s,[0._dp,0._dp,1._dp])
      if(s/=0)error stop 'antipodal status'
      if(norm2(v1/vu-[0._dp,1._dp,0._dp])>1e-10_dp.or.norm2(v2/vu-[0._dp,-1._dp,0._dp])>1e-10_dp) &
         error stop 'antipodal scaled truth'
   enddo
   r1=[1._dp,0._dp,0._dp];r2=[0._dp,2._dp,0._dp]
   call excali_solve(r1,r2,4*sqrt(2._dp)/3,1._dp,1,0,-1,v1,v2,s)
   if(s/=0.or.norm2(v1-[0._dp,sqrt(2._dp),0._dp])>1e-10_dp.or. &
      norm2(v2-[-1/sqrt(2._dp),1/sqrt(2._dp),0._dp])>1e-10_dp)error stop 'parabolic truth'
   do i=1,100
      e=1.01_dp+3._dp*i/100;a=1;hh=.001_dp+1.5_dp*i/100
      rr=a*(e*cosh(hh)-1);nu=2*atan(sqrt((e+1)/(e-1))*tanh(hh/2))
      r1=[a*(e-1),0._dp,0._dp];r2=rr*[cos(nu),sin(nu),0._dp]
      w1=[0._dp,sqrt((e+1)/(a*(e-1))),0._dp]
      w2=[-sin(nu),e+cos(nu),0._dp]/sqrt(a*(e*e-1))
      t=a**1.5_dp*(e*sinh(hh)-hh)
      call excali_one(r1,r2,t,1,v1,v2,s,h,0,-1)
      if(s/=0.or.norm2(v1-w1)>1e-10_dp.or.norm2(v2-w2)>1e-10_dp)error stop 'hyperbolic truth'
   enddo
   call excali_one(r1,r2,-1._dp,1,v1,v2,s,h)
   if(s/=-1)error stop 'invalid time'
   call excali_solve([1._dp,0._dp,0._dp],[-1._dp,0._dp,0._dp],pi,1._dp,1,0,-1,v1,v2,s)
   if(s/=-2)error stop 'plane required'
   print *, 'PASS: EXCALI independent conic truth, branches, antipodal geometry, units and input statuses'
end program
