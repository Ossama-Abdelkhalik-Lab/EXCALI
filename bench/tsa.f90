! Time excali_solve_all against a plain excali_one loop over the same (m, branch) list.
program tsa
  use excali_solver
  use iso_fortran_env, only: real64,int64
  implicit none
  integer,parameter::dp=real64,nc=2000
  real(dp)::r1(3,nc),r2(3,nc),tof(nc),u(4),th,rr,vv1(3,401),vv2(3,401),w1(3),w2(3),pme,c,s,acc,ta,tn
  integer::nrev(401),br(401),st(401),ns,nmax,h,i,j,k,dir(nc),sw,nlist(nc),ntarget(5)=[1,2,5,10,50],it,rep,ntot
  integer(int64)::c0,c1,cr
  call random_seed(put=[(4,i=1,64)])
  do it=1,5
    do i=1,nc
      call random_number(u)
      th=0.02_dp+6.24_dp*u(1);rr=0.2_dp+0.8_dp*u(2);dir(i)=merge(1,-1,u(3)<.5_dp)
      r1(:,i)=[1._dp,0._dp,0._dp];r2(:,i)=rr*[cos(th),sin(th),0._dp]
      c=sqrt(1+rr*rr-2*rr*cos(th));s=(1+rr+c)/2;pme=acos(-1._dp)*s**1.5_dp/sqrt(2._dp)
      tof(i)=(ntarget(it)+0.3_dp+0.6_dp*u(4))*pme
    enddo
    ta=huge(1._dp);tn=huge(1._dp);acc=0
    do rep=1,7
      call system_clock(c0,cr);ntot=0
      do i=1,nc
        call excali_solve_all(r1(:,i),r2(:,i),tof(i),dir(i),200,ns,nmax,nrev,br,vv1,vv2,st,h)
        nlist(i)=ns;ntot=ntot+ns;acc=acc+vv1(1,ns)
      enddo
      call system_clock(c1);ta=min(ta,real(c1-c0,dp)/cr/ntot*1e9_dp)
      call system_clock(c0,cr)
      do i=1,nc
        call excali_one(r1(:,i),r2(:,i),tof(i),dir(i),w1,w2,sw,h,0,-1);acc=acc+w1(1)
        do k=1,(nlist(i)-1)/2+1          ! the same list plus the n_up check excali_solve_all makes
          do j=-1,1,2
            call excali_one(r1(:,i),r2(:,i),tof(i),dir(i),w1,w2,sw,h,k,j);acc=acc+w1(1)
          enddo
        enddo
      enddo
      call system_clock(c1);tn=min(tn,real(c1-c0,dp)/cr/ntot*1e9_dp)
    enddo
    print '(a,i3,a,f6.1,a,f8.1,a,f8.1,a,f6.3,a,es9.2)','N~',ntarget(it),'  sols/case',real(ntot,dp)/nc,'  solve_all ns/sol',ta,'  naive ns/sol',tn,'  speedup',tn/ta,'  chk',acc
  enddo
end program
