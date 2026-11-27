! excali_solve_all against excali_one: every (m, branch) entry identical,
! N_max + 1 infeasible on both branches, entries in the documented order.
program test_excali_solve_all
   use iso_fortran_env, only: real64
   use excali_solver
   implicit none
   integer,parameter::dp=real64,cap=300
   real(dp)::r1(3),r2(3),tof,u(6),th,rr,v1(3,2*cap+1),v2(3,2*cap+1),w1(3),w2(3)
   integer::nrev(2*cap+1),br(2*cap+1),st(2*cap+1),ns,nmax,h,i,j,dir,sw,bad,tot
   call random_seed(put=[(17,i=1,64)])
   bad=0;tot=0
   do i=1,3000
      call random_number(u)
      th=0.02_dp+6.24_dp*u(1);rr=0.2_dp+0.8_dp*u(2);dir=merge(1,-1,u(3)<.5_dp)
      r1=[1._dp,0._dp,0._dp];r2=rr*[cos(th),sin(th)*cos(u(6)),sin(th)*sin(u(6))]
      tof=10**(-0.5_dp+3*u(4))
      call excali_solve_all(r1,r2,tof,dir,cap,ns,nmax,nrev,br,v1,v2,st,h)
      if(nmax<0.or.ns/=2*nmax+1)bad=bad+1
      do j=1,ns
         tot=tot+1
         if(j>1.and.(nrev(j)/=(j)/2.or.br(j)/=merge(-1,1,mod(j,2)==0)))bad=bad+1
         call excali_one(r1,r2,tof,dir,w1,w2,sw,h,nrev(j),br(j))
         if(sw/=st(j))bad=bad+1
         if(sw>=0.and.sw<=2)then
            if(any(w1/=v1(:,j)).or.any(w2/=v2(:,j)))bad=bad+1
         endif
         if(j>1.and..not.(sw>=0.and.sw<=2))bad=bad+1
      enddo
      do j=-1,1,2
         call excali_one(r1,r2,tof,dir,w1,w2,sw,h,nmax+1,j)
         if(sw/=5)bad=bad+1
      enddo
   enddo
   print '(a,i0,a,i0)',' solutions checked: ',tot,'  failures: ',bad
   if(bad>0)error stop 'excali_solve_all disagrees with excali_one'
end program
