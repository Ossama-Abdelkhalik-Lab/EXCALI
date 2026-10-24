! Training workload for profile-guided builds (ELAM_PGO=GENERATE): random
! transfers across the solver's regimes. Times are drawn relative to the
! parabolic time, so short/long-way ellipses, hyperbolas, near-parabolic
! arcs and multi-revolution branches (feasible or not) all occur.
program pgo_train
   use iso_fortran_env, only: real64
   use elam_solver, only: elam_solve
   implicit none
   integer,parameter::dp=real64
   real(dp),parameter::pi=acos(-1._dp)
   real(dp)::r1(3),r2(3),v1(3),v2(3),u(6),th,rho,tp,tof,s,c
   integer::i,status,d,m,br,solved
   call random_seed(put=[(12345+i,i=1,64)])
   solved=0
   do i=1,400000
      call random_number(u)
      th=0.02_dp+(2*pi-0.04_dp)*u(1)
      rho=0.05_dp+0.95_dp*u(2)
      r1=[1._dp,0._dp,0._dp];r2=rho*[cos(th),sin(th),0._dp]
      d=1;if(th>pi)d=-1
      ! Parabolic time scale from Euler/Lambert: s and c of the triangle.
      c=norm2(r2-r1);s=(1+rho+c)/2
      tp=sqrt(2._dp)/3*(s**1.5_dp-sign(1._dp,pi-th)*(s-c)**1.5_dp)
      m=0;br=-1
      if(u(3)<0.3_dp)then
         m=1+int(4*u(4));if(u(5)<0.5_dp)br=1
         tof=(m*2*pi+th)*(0.5_dp+u(6))*s**1.5_dp
      else
         tof=tp*exp(4*u(4)-1.5_dp)
      endif
      call elam_solve(r1,r2,tof,1._dp,d,m,br,v1,v2,status)
      if(status==0.or.status==1)solved=solved+1
   enddo
   print '(a,i0,a)','pgo_train: ',solved,' of 400000 transfers solved'
end program
