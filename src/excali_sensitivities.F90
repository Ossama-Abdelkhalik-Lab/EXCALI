! First- and second-order sensitivities of EXCALI's solution,
!   p = [v1, v2](z),  z = [r1, r2, tof]  (mu fixed),
! as the 6x7 Jacobian dp/dz and the 6x7x7 Hessian d2p/dz2.
!
! Method: implicit differentiation at the converged root. EXCALI solves the
! time equation G(alpha, r1, r2) = sqrt(y)*(y*F(alpha) + a) - sqrt(mu)*tof = 0
! for alpha in one of its coordinates; y depends on the positions only through
! the two invariants r1+r2 (in y0, yend) and a = sqrt(2)*sqrt(r1*r2)*cs_s.
! G and the velocity formulas are evaluated once at the root in second-order
! forward-mode automatic differentiation over the seven variables
! (r1, r2, alpha), using the same cancellation-free forms as the solver; then
!   alpha_z  = -G_z/G_alpha,
!   alpha_zz = -(G_zz + G_za alpha_z' + alpha_z G_za' + G_aa alpha_z alpha_z')/G_alpha,
!   J = p_z + p_a alpha_z,
!   H = p_zz + p_za alpha_z' + alpha_z p_za' + p_aa alpha_z alpha_z' + p_a alpha_zz.
!
! Coordinates: X = +-x**2 for v or w below 1.6 (m = 0), where the time
! equation is regular through X = 0 and the series for F keeps its
! derivatives exact (the closed form's F', F'' cancel like 1/x**2 there); v; pi - v near the pole; w for
! hyperbolas, and the solver's xi = wb - w next to the short-way branch point
! wb, where y = 2b sinh(wb - xi/2) sinh(xi/2) and its derivatives keep full
! precision. y's value is always the solver's (cancellation-free). Half-angle functions of the
! transfer angle are signed about the transfer-plane normal h: cs_s = cos(theta/2)
! (negative on the long way), ss = sin(theta/2) >= 0, from 1+-cos(theta) on the
! side where that does not cancel and sin(theta)/(2*the other) otherwise, so
! they are smooth through 180 degrees when h is fixed (plane_normal). With
! plane_normal at exactly collinear endpoints h is held fixed; out-of-plane
! derivatives then refer to that fixed plane.
module excali_sens_util
   use iso_fortran_env, only: real64
   implicit none
   private
   public :: dp, pi, sqrt2, fell, fhyp, fser
   integer,parameter :: dp=real64
   real(dp),parameter :: pi=acos(-1._dp), sqrt2=sqrt(2._dp)
   ! Parabolic series of (x - sin x cos x)/sin(x)**3 in X = x**2, 40 terms
   ! (exact rationals rounded; highest first), accurate to rounding for |X| < 2.56.
   real(dp),parameter :: pall(40)=[ &
      1.07005216020305206e-36_dp, 1.00371021500316345e-35_dp, 9.40233813648184982e-35_dp, 8.79542327801632922e-34_dp, &
      8.21556968602699418e-33_dp, 7.66200749482881528e-32_dp, 7.13398449677639024e-31_dp, 6.63076578420779214e-30_dp, &
      6.15163339809479518e-29_dp, 5.69588598328062163e-28_dp, 5.26283845018836362e-27_dp, 4.85182164288525860e-26_dp, &
      4.46218201338992527e-25_dp, 4.09328130211158755e-24_dp, 3.74449622431223251e-23_dp, 3.41521816248448676e-22_dp, &
      3.10485286453982803e-21_dp, 2.81282014770355301e-20_dp, 2.53855360801467948e-19_dp, 2.28150033533065582e-18_dp, &
      2.04112063373833767e-17_dp, 1.81688774727389515e-16_dp, 1.60828759085431170e-15_dp, 1.41481848631905828e-14_dp, &
      1.23599090346256394e-13_dp, 1.07132720587347656e-12_dp, 9.20361401176818234e-12_dp, 7.82638894546463940e-11_dp, &
      6.57716241989020179e-10_dp, 5.45160892374639085e-09_dp, 4.44550883834929939e-08_dp, 3.55474389799148566e-07_dp, &
      2.77528805919328775e-06_dp, 2.10318506449458828e-05_dp, 1.53448954044192130e-04_dp, 1.06511544011544013e-03_dp, &
      6.90476190476190486e-03_dp, 4.04761904761904781e-02_dp, 2.00000000000000011e-01_dp, 6.66666666666666630e-01_dp]

contains
   ! ------------------------------------------------- scalar coordinate jets
   ! F(x) = (x + m*pi - sin x cos x)/(sqrt2 sin(x)**3) and two derivatives
   ! (s'F + 3cF = sqrt2), at x with s = sin x, c = cos x.
   pure subroutine fell(x,m,s,c,f)
      real(dp),intent(in)::x,s,c
      integer,intent(in)::m
      real(dp),intent(out)::f(0:2)
      real(dp)::r
      r=1/s
      f(0)=(x+m*pi-s*c)/sqrt2*r**3
      f(1)=(sqrt2-3*c*f(0))*r
      f(2)=(3*s*f(0)-4*c*f(1))*r
   end subroutine
   ! Hyperbolic F(w) = (sinh w cosh w - w)/(sqrt2 sinh(w)**3), w >= 0.75.
   pure subroutine fhyp(w,f)
      real(dp),intent(in)::w
      real(dp),intent(out)::f(0:2)
      real(dp)::s,c
      s=sinh(w);c=cosh(w)
      f(0)=(s*c-w)/(sqrt2*s**3)
      f(1)=(sqrt2-3*c*f(0))/s
      f(2)=-(4*c*f(1)+3*s*f(0))/s
   end subroutine
   ! F(X) from the parabolic series and C(X) = (1-cos sqrt X)/X, |X| < 2.56.
   pure subroutine fser(xx,f,cc)
      real(dp),intent(in)::xx
      real(dp),intent(out)::f(0:2),cc(0:2)
      real(dp)::t
      integer::j
      f=0
      do j=1,40
         f(2)=f(2)*xx+2*f(1);f(1)=f(1)*xx+f(0);f(0)=f(0)*xx+pall(j)
      enddo
      f=f/sqrt2
      ! C(X) = sum (-1)**j X**j/(2j+2)!, highest first
      cc=0
      do j=24,0,-1
         t=(-1)**j/gamma(2._dp*j+3)
         cc(2)=cc(2)*xx+2*cc(1);cc(1)=cc(1)*xx+cc(0);cc(0)=cc(0)*xx+t
      enddo
   end subroutine

end module excali_sens_util

module excali_sens_d1
   use excali_sens_util
   implicit none
   private
   public :: sens_eval
#define ORD2 0
#include "excali_sens_ad.inc"
#undef ORD2
end module excali_sens_d1

module excali_sens_d2
   use excali_sens_util
   implicit none
   private
   public :: sens_eval
#define ORD2 1
#include "excali_sens_ad.inc"
#undef ORD2
end module excali_sens_d2

module excali_sens
   use excali_sens_util, only: dp
   use excali_sens_d1, only: eval1 => sens_eval
   use excali_sens_d2, only: eval2 => sens_eval
   use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use excali_solver, only: excali_one
   implicit none
   private
   public :: excali_sens_solve
   real(dp),parameter :: xlim=1.6_dp   ! zero-rev roots below this v or w use X = +-x**2
contains
   !> Solves like excali_solve and returns the velocities with their 6x7
   !! Jacobian w.r.t. [r1, r2, tof] and, if hess is present, the 6x7x7
   !! Hessian (hess(i,j,k) = d2 p_i/dz_j dz_k). status as excali_solve;
   !! jac/hess are NaN unless 0 <= status <= 2.
   subroutine excali_sens_solve(r1,r2,tof,mu,direction,nrev,branch,v1,v2,jac,status,hess,plane_normal)
      real(dp),intent(in)::r1(3),r2(3),tof,mu
      integer,intent(in)::direction,nrev,branch
      real(dp),intent(out)::v1(3),v2(3),jac(6,7)
      integer,intent(out)::status
      real(dp),optional,intent(out)::hess(6,7,7)
      real(dp),optional,intent(in)::plane_normal(3)
      real(dp)::rs(4),x,xx,smu,e1(3),e2(3),hn,hfix(3)
      integer::half,kind
      logical::fixedh
      jac=ieee_value(0._dp,ieee_quiet_nan)
      if(present(hess))hess=jac(1,1)
      smu=sqrt(mu)
      call excali_one(r1,r2,tof*smu,direction,v1,v2,status,half,nrev,branch,plane_normal,rs)
      if(status<0.or.status>2)return
      v1=v1*smu;v2=v2*smu
      x=rs(2);xx=0
      kind=nint(rs(1))                     ! 0 X, 1 v, 2 pi-v, -1 w, -2 xi = wb - w
      if(nrev==0.and.abs(kind)==1.and.x<xlim)then
         xx=merge(x*x,-x*x,kind>0);kind=0
      endif
      ! Collinear endpoints with plane_normal: the solver's fixed plane.
      e1=r1/sqrt(dot_product(r1,r1));e2=r2/sqrt(dot_product(r2,r2))
      hn=norm2([e1(2)*e2(3)-e1(3)*e2(2),e1(3)*e2(1)-e1(1)*e2(3),e1(1)*e2(2)-e1(2)*e2(1)])
      fixedh=hn<32*epsilon(hn).and.present(plane_normal)
      hfix=0
      if(fixedh)hfix=plane_normal*(direction/norm2(plane_normal))
      if(present(hess))then
         call eval2(r1,r2,tof,smu,direction,nrev,kind,x,xx,rs(3),rs(4)>.5_dp,fixedh,hfix,jac,hess)
      else
         call eval1(r1,r2,tof,smu,direction,nrev,kind,x,xx,rs(3),rs(4)>.5_dp,fixedh,hfix,jac)
      endif
   end subroutine
end module excali_sens
