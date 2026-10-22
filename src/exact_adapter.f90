! Adapters only: imported Exact v2.2 kernels remain byte-for-byte unchanged.
module exact_adapter
   use iso_fortran_env, only: real64
   use ieee_arithmetic
   use exact_original_serial, only: serial_solve=>rectsolution, serial_init=>init_cache
   use exact_original_vector, only: vector_solve=>rectsolution, vector_init=>init_cache
   implicit none
   private
   integer,parameter::dp=real64
   public::exact_initialize,exact_one,exact_batch
contains
   subroutine exact_initialize()
      ! Call once, serially, before entering any parallel solve region.
      call serial_init()
      call vector_init()
   end subroutine

   subroutine exact_one(vector,r1,r2,tof,d,m,b,v1,v2,status,half)
      logical,intent(in)::vector
      real(dp),intent(in)::r1(3),r2(3),tof
      integer,intent(in)::d,m,b
      real(dp),intent(out)::v1(3),v2(3)
      integer,intent(out)::status,half
      real(dp)::h(3),a,c,x(3,1),y(3,1)
      v1=ieee_value(0._dp,ieee_quiet_nan);v2=v1;status=-1;half=0
      if(.not.all(ieee_is_finite(r1)).or..not.all(ieee_is_finite(r2)))return
      if(.not.ieee_is_finite(tof).or.tof<=0.or.abs(d)/=1.or.abs(b)/=1.or.m<0)return
      a=norm2(r1);c=norm2(r2)
      if(a<=0.or.c<=0)return
      h=cross(r1/a,r2/c);a=norm2(h)
      if(a<32*epsilon(a))then
         status=-2;half=-1;return
      endif
      if(a<1e-4_dp)half=1
      status=-6 ! unsupported regime: never silently solve M=0 for M>0
      if(m/=0)return
      h=h*(d/a)
      if(vector)then
         call vector_solve(reshape(r1,[3,1]),reshape(r2,[3,1]),[tof],1._dp,h,x,y,1)
      else
         call serial_solve(reshape(r1,[3,1]),reshape(r2,[3,1]),[tof],1._dp,h,x,y,1)
      endif
      v1=x(:,1);v2=y(:,1);status=-4
      if(all(ieee_is_finite(v1)).and.all(ieee_is_finite(v2)))status=0
   end subroutine

   subroutine exact_batch(vector,n,r1,r2,tof,d,m,b,v1,v2,status,half)
      logical,intent(in)::vector
      integer,intent(in)::n,d(n),m(n),b(n)
      real(dp),intent(in)::r1(3,n),r2(3,n),tof(n)
      real(dp),intent(out)::v1(3,n),v2(3,n)
      integer,intent(out)::status(n),half(n)
      integer::i,j,k,last,nb
      integer,parameter::chunk=256
      ! Rotate each valid case into its own xy plane. A shared +z H then
      ! represents the requested orientation for every original batch row.
      real(dp)::p(3,chunk),q(3,chunk),t(chunk),x(3,chunk),y(3,chunk)
      real(dp)::e(3,chunk),f(3,chunk),h(3),a,c,hn
      logical::valid(chunk)
      !$omp parallel do private(i,j,k,last,nb,p,q,t,x,y,e,f,h,a,c,hn,valid) schedule(static)
      do i=1,n,chunk
         last=min(n,i+chunk-1);nb=last-i+1
         do j=1,nb
            k=i+j-1;status(k)=-1;half(k)=0;valid(j)=.false.
            p(:,j)=[1._dp,0._dp,0._dp];q(:,j)=[0._dp,1._dp,0._dp];t(j)=2._dp
            if(.not.all(ieee_is_finite(r1(:,k))).or..not.all(ieee_is_finite(r2(:,k))))cycle
            if(.not.ieee_is_finite(tof(k)).or.tof(k)<=0.or.abs(d(k))/=1.or.abs(b(k))/=1.or.m(k)<0)cycle
            a=norm2(r1(:,k));c=norm2(r2(:,k));if(a<=0.or.c<=0)cycle
            e(:,j)=r1(:,k)/a;h=cross(e(:,j),r2(:,k)/c);hn=norm2(h)
            if(hn<32*epsilon(a))then
               status(k)=-2;half(k)=-1;cycle
            endif
            if(hn<1e-4_dp)half(k)=1
            status(k)=-6;if(m(k)/=0)cycle
            h=h*(d(k)/hn);f(:,j)=cross(h,e(:,j))
            p(:,j)=[a,0._dp,0._dp]
            q(:,j)=[dot_product(r2(:,k),e(:,j)),dot_product(r2(:,k),f(:,j)),0._dp]
            t(j)=tof(k);valid(j)=.true.
         enddo
         if(vector)then
            call vector_solve(p(:,:nb),q(:,:nb),t(:nb),1._dp,[0._dp,0._dp,1._dp],x(:,:nb),y(:,:nb),nb)
         else
            call serial_solve(p(:,:nb),q(:,:nb),t(:nb),1._dp,[0._dp,0._dp,1._dp],x(:,:nb),y(:,:nb),nb)
         endif
         do j=1,nb
            k=i+j-1;v1(:,k)=ieee_value(0._dp,ieee_quiet_nan);v2(:,k)=v1(:,k)
            if(.not.valid(j))cycle
            v1(:,k)=x(1,j)*e(:,j)+x(2,j)*f(:,j);v2(:,k)=y(1,j)*e(:,j)+y(2,j)*f(:,j)
            status(k)=-4
            if(all(ieee_is_finite(v1(:,k))).and.all(ieee_is_finite(v2(:,k))))status(k)=0
         enddo
      enddo
      !$omp end parallel do
   end subroutine
   pure function cross(a,b) result(c)
      real(dp),intent(in)::a(3),b(3)
      real(dp)::c(3)
      c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
   end function
end module
