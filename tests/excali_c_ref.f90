! Reference for tests/test_excali_c.c: the same requests through direct
! Fortran calls of the excali_solver module (default kinds, Fortran arrays).
module excali_c_ref
   use iso_c_binding, only: c_int, c_double
   use iso_fortran_env, only: real64
   use excali_solver, only: excali_solve, excali_batch
   implicit none
contains
   subroutine ref_solve(r1,r2,tof,mu,direction,nrev,branch,v1,v2,status) bind(C,name='ref_solve')
      real(c_double),intent(in)::r1(3),r2(3)
      real(c_double),value::tof,mu
      integer(c_int),value::direction,nrev,branch
      real(c_double),intent(out)::v1(3),v2(3)
      integer(c_int),intent(out)::status
      real(real64)::a(3),b(3),u(3),w(3)
      integer::s
      a=r1;b=r2
      call excali_solve(a,b,real(tof,real64),real(mu,real64),int(direction),int(nrev),int(branch),u,w,s)
      v1=u;v2=w;status=s
   end subroutine
   subroutine ref_solve_plane(r1,r2,tof,mu,direction,nrev,branch,pn,v1,v2,status) bind(C,name='ref_solve_plane')
      real(c_double),intent(in)::r1(3),r2(3),pn(3)
      real(c_double),value::tof,mu
      integer(c_int),value::direction,nrev,branch
      real(c_double),intent(out)::v1(3),v2(3)
      integer(c_int),intent(out)::status
      real(real64)::a(3),b(3),p(3),u(3),w(3)
      integer::s
      a=r1;b=r2;p=pn
      call excali_solve(a,b,real(tof,real64),real(mu,real64),int(direction),int(nrev),int(branch),u,w,s,p)
      v1=u;v2=w;status=s
   end subroutine
   subroutine ref_batch(n,r1,r2,tof,direction,nrev,branch,v1,v2,status,half) bind(C,name='ref_batch')
      integer(c_int),value::n
      real(c_double),intent(in)::r1(3,n),r2(3,n),tof(n)
      integer(c_int),intent(in)::direction(n),nrev(n),branch(n)
      real(c_double),intent(out)::v1(3,n),v2(3,n)
      integer(c_int),intent(out)::status(n),half(n)
      real(real64),allocatable::a(:,:),b(:,:),t(:),u(:,:),w(:,:)
      integer,allocatable::d(:),m(:),br(:),s(:),h(:)
      a=r1;b=r2;t=tof;d=direction;m=nrev;br=branch
      allocate(u(3,n),w(3,n),s(n),h(n))
      call excali_batch(int(n),a,b,t,d,u,w,s,h,m,br)
      v1=u;v2=w;status=s;half=h
   end subroutine
end module
