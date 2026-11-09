! Adapters only: the imported EXCALI v0 kernel (original/excali_v0.F90, built
! with -DEXCALI_V0_STANDALONE) is unchanged apart from its names. Same input
! checks and statuses as excali_v1_adapter. v0 solves zero-revolution
! elliptic transfers only and returns NaN where it cannot solve (hyperbolic
! times, degenerate plane), reported here as status -4.
module excali_v0_adapter
   use iso_fortran_env, only: real64
   use ieee_arithmetic
   use excali_v0_mod, only: v0_solve=>excali_v0_solution, v0_init=>init_cache
   implicit none
   private
   integer,parameter::dp=real64
   public::excali_v0_initialize,excali_v0_one,excali_v0_batch
contains
   subroutine excali_v0_initialize()
      ! Call once, serially: v0 otherwise fills its cache lazily on the first
      ! solve, which races when cases are solved in parallel.
      call v0_init()
   end subroutine

   ! Shared checks; on success h is the unit normal times direction.
   subroutine check(r1,r2,tof,d,m,b,h,status,half,ok)
      real(dp),intent(in)::r1(3),r2(3),tof
      integer,intent(in)::d,m,b
      real(dp),intent(out)::h(3)
      integer,intent(out)::status,half
      logical,intent(out)::ok
      real(dp)::a,c
      ok=.false.;status=-1;half=0;h=0
      if(.not.all(ieee_is_finite(r1)).or..not.all(ieee_is_finite(r2)))return
      if(.not.ieee_is_finite(tof).or.tof<=0.or.abs(d)/=1.or.abs(b)/=1.or.m<0)return
      a=norm2(r1);c=norm2(r2)
      if(a<=0.or.c<=0)return
      h=cross(r1/a,r2/c);a=norm2(h)
      if(a<32*epsilon(a))then
         status=-2;half=-1;return
      endif
      if(a<1e-4_dp)half=1
      status=-6 ! zero-revolution ellipses only
      if(m/=0)return
      h=h*(d/a);ok=.true.
   end subroutine

   subroutine excali_v0_one(r1,r2,tof,d,m,b,v1,v2,status,half)
      real(dp),intent(in)::r1(3),r2(3),tof
      integer,intent(in)::d,m,b
      real(dp),intent(out)::v1(3),v2(3)
      integer,intent(out)::status,half
      real(dp)::h(3),x(3,1),y(3,1)
      logical::ok
      v1=ieee_value(0._dp,ieee_quiet_nan);v2=v1
      call check(r1,r2,tof,d,m,b,h,status,half,ok)
      if(.not.ok)return
      call v0_solve(reshape(r1,[3,1]),reshape(r2,[3,1]),[tof],1._dp,h,x,y,1)
      v1=x(:,1);v2=y(:,1);status=-4
      if(all(ieee_is_finite(v1)).and.all(ieee_is_finite(v2)))status=0
   end subroutine

   ! All valid cases rotated into their own xy plane (a shared +z normal then
   ! carries each case's direction) and solved in one v0 call, as v0 batches.
   subroutine excali_v0_batch(n,r1,r2,tof,d,m,b,v1,v2,status,half)
      integer,intent(in)::n,d(n),m(n),b(n)
      real(dp),intent(in)::r1(3,n),r2(3,n),tof(n)
      real(dp),intent(out)::v1(3,n),v2(3,n)
      integer,intent(out)::status(n),half(n)
      real(dp),allocatable::p(:,:),q(:,:),t(:),x(:,:),y(:,:),e(:,:),f(:,:)
      real(dp)::h(3),a
      logical::ok(n)
      integer::k
      allocate(p(3,n),q(3,n),t(n),x(3,n),y(3,n),e(3,n),f(3,n))
      do k=1,n
         p(:,k)=[1._dp,0._dp,0._dp];q(:,k)=[0._dp,1._dp,0._dp];t(k)=2._dp
         call check(r1(:,k),r2(:,k),tof(k),d(k),m(k),b(k),h,status(k),half(k),ok(k))
         if(.not.ok(k))cycle
         a=norm2(r1(:,k));e(:,k)=r1(:,k)/a;f(:,k)=cross(h,e(:,k))
         p(:,k)=[a,0._dp,0._dp]
         q(:,k)=[dot_product(r2(:,k),e(:,k)),dot_product(r2(:,k),f(:,k)),0._dp]
         t(k)=tof(k)
      enddo
      call v0_solve(p,q,t,1._dp,[0._dp,0._dp,1._dp],x,y,n)
      do k=1,n
         v1(:,k)=ieee_value(0._dp,ieee_quiet_nan);v2(:,k)=v1(:,k)
         if(.not.ok(k))cycle
         v1(:,k)=x(1,k)*e(:,k)+x(2,k)*f(:,k);v2(:,k)=y(1,k)*e(:,k)+y(2,k)*f(:,k)
         status(k)=-4
         if(all(ieee_is_finite(v1(:,k))).and.all(ieee_is_finite(v2(:,k))))status(k)=0
      enddo
   end subroutine

   pure function cross(a,b) result(c)
      real(dp),intent(in)::a(3),b(3)
      real(dp)::c(3)
      c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
   end function
end module
