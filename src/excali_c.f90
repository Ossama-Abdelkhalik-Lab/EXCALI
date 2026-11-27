! C interface to EXCALI (bind(C)). Arrays of 3-vectors are column-major
! (3,n): the same bytes as a C array double[n][3] or a C-contiguous NumPy
! array of shape (n,3), so no copy is made in either direction. Integer
! arguments are C int. Every routine writes all outputs; check status first
! (see README for the codes). include/excali.h declares these routines.
module excali_c
   use iso_c_binding, only: c_int, c_double, c_intptr_t, c_loc
   use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use excali_solver, only: excali_solve, excali_batch
   implicit none
   private
   ! Integer arrays are passed straight to excali_batch, so C int must be the
   ! default integer kind (true for gfortran); this fails to compile otherwise.
   integer,parameter::int_kind_check(merge(1,-1,c_int==kind(0)))=0
   public :: excali_solve_c, excali_solve_plane_c, excali_batch_c, excali_batch_mu_c, excali_c_data_address

contains

   ! One dimensional case: positions in any distance unit, tof in the time
   ! unit and mu in distance**3/time**2 (the Fortran excali_solve).
   subroutine excali_solve_c(r1,r2,tof,mu,direction,nrev,branch,v1,v2,status) bind(C,name='excali_solve')
      real(c_double),intent(in)::r1(3),r2(3)
      real(c_double),value,intent(in)::tof,mu
      integer(c_int),value,intent(in)::direction,nrev,branch
      real(c_double),intent(out)::v1(3),v2(3)
      integer(c_int),intent(out)::status
      integer::s
      call excali_solve(r1,r2,tof,mu,int(direction),int(nrev),int(branch),v1,v2,s)
      status=int(s,c_int)
   end subroutine

   ! As excali_solve, with the transfer-plane normal for exactly antipodal
   ! endpoints (perpendicular to both position vectors).
   subroutine excali_solve_plane_c(r1,r2,tof,mu,direction,nrev,branch,plane_normal,v1,v2,status) &
         bind(C,name='excali_solve_plane')
      real(c_double),intent(in)::r1(3),r2(3),plane_normal(3)
      real(c_double),value,intent(in)::tof,mu
      integer(c_int),value,intent(in)::direction,nrev,branch
      real(c_double),intent(out)::v1(3),v2(3)
      integer(c_int),intent(out)::status
      integer::s
      call excali_solve(r1,r2,tof,mu,int(direction),int(nrev),int(branch),v1,v2,s,plane_normal)
      status=int(s,c_int)
   end subroutine

   ! n normalized cases (mu = 1, any distance unit consistent with tof):
   ! the Fortran excali_batch. half(i) is the geometry flag (-1/0/1).
   subroutine excali_batch_c(n,r1,r2,tof,direction,nrev,branch,v1,v2,status,half) bind(C,name='excali_batch')
      integer(c_int),value,intent(in)::n
      real(c_double),intent(in)::r1(3,*),r2(3,*),tof(*)
      integer(c_int),intent(in)::direction(*),nrev(*),branch(*)
      real(c_double),intent(out)::v1(3,*),v2(3,*)
      integer(c_int),intent(out)::status(*),half(*)
      if(n<=0)return
      ! Assumed-size actuals to explicit-shape dummies: sequence association,
      ! never a copy.
      call excali_batch(int(n),r1,r2,tof,direction,v1,v2,status,half,nrev,branch)
   end subroutine

   ! n cases with one gravitational parameter mu: time is rescaled by
   ! sqrt(mu) (a temporary copy of tof) and velocities by sqrt(mu) in place;
   ! mu = 1 calls excali_batch directly. Results agree with per-case
   ! excali_solve to rounding (statuses identical).
   subroutine excali_batch_mu_c(n,r1,r2,tof,mu,direction,nrev,branch,v1,v2,status,half) bind(C,name='excali_batch_mu')
      integer(c_int),value,intent(in)::n
      real(c_double),value,intent(in)::mu
      real(c_double),intent(in)::r1(3,*),r2(3,*),tof(*)
      integer(c_int),intent(in)::direction(*),nrev(*),branch(*)
      real(c_double),intent(out)::v1(3,*),v2(3,*)
      integer(c_int),intent(out)::status(*),half(*)
      real(c_double),allocatable::tt(:)
      real(c_double)::sm
      integer::i
      if(n<=0)return
      if(.not.(mu>0.and.mu<huge(mu)))then
         ! Invalid mu: every case is invalid input, outputs NaN.
         v1(:,1:n)=ieee_value(1._c_double,ieee_quiet_nan);v2(:,1:n)=v1(:,1:n)
         status(1:n)=-1;half(1:n)=0
         return
      endif
      if(mu==1)then
         call excali_batch_c(n,r1,r2,tof,direction,nrev,branch,v1,v2,status,half)
         return
      endif
      sm=sqrt(mu)
      allocate(tt(n))
      tt=tof(1:n)*sm
      call excali_batch(int(n),r1,r2,tt,direction,v1,v2,status,half,nrev,branch)
      do i=1,n
         if(status(i)>=0.and.status(i)<=2)then
            v1(:,i)=v1(:,i)*sm;v2(:,i)=v2(:,i)*sm
         endif
      enddo
   end subroutine

   ! Diagnostic: the address Fortran sees for x (used by the tests to check
   ! that arrays are passed without a copy).
   integer(c_intptr_t) function excali_c_data_address(x) bind(C,name='excali_c_data_address')
      real(c_double),intent(in),target::x(*)
      excali_c_data_address=transfer(c_loc(x(1)),excali_c_data_address)
   end function

end module excali_c
