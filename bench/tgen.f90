! Timing on generated cases (no data files): zero-rev mix and multi-rev, single
! (excali_one) and batch (excali_batch), best of 9, ns per case.
program tgen
  use excali_solver
  use iso_fortran_env, only: real64,int64
  implicit none
  integer,parameter::dp=real64,nz=50000,nm=10000
  real(dp),allocatable::r1(:,:),r2(:,:),tof(:),v1(:,:),v2(:,:)
  integer,allocatable::dir(:),nr(:),br(:),st(:),hf(:)
  real(dp)::u(7),th,rr,c,s,cs,p,a,q,y0,tp,best1,bestb,chk,rot(3,3),ph
  integer::i,set,n,rep
  integer(int64)::c0,c1,cr
  character(8)::lab
  call random_seed(put=[(21,i=1,64)])
  do set=1,2
    n=merge(nz,nm,set==1); lab=merge('zero_rev','multirev',set==1)
    if(allocated(r1))deallocate(r1,r2,tof,v1,v2,dir,nr,br,st,hf)
    allocate(r1(3,n),r2(3,n),tof(n),v1(3,n),v2(3,n),dir(n),nr(n),br(n),st(n),hf(n))
    do i=1,n
      call random_number(u)
      th=0.02_dp+6.24_dp*u(1);rr=0.2_dp+0.8_dp*u(2);dir(i)=merge(1,-1,u(3)<.5_dp)
      ph=6.283_dp*u(6)
      r1(:,i)=[cos(ph),sin(ph),0._dp]
      r2(:,i)=rr*[cos(ph+th),sin(ph+th)*cos(u(7)),sin(ph+th)*sin(u(7))]
      c=sqrt(1+rr*rr-2*rr*cos(th));s=(1+rr+c)/2
      if(set==1)then
        ! parabolic time from the transfer frame, then t/tp in 10**[-1.5,1.5]
        cs=sqrt((1+rr+2*rr*cos(th))/(4*rr))*0+abs(cos(th/2));p=sqrt(rr)
        a=dir(i)*sqrt(2._dp)*p*cs;q=(1-rr)/(1+sqrt(rr));y0=q*q+2*p*(1-cs*cs)/(1+cs);if(dir(i)<0)y0=q*q+2*p*(1+cs)
        tp=sqrt(y0)*(y0*sqrt(2._dp)/3+a)
        tof(i)=tp*10**(-1.5_dp+3*u(4));nr(i)=0;br(i)=-1
      else
        nr(i)=1+int(10*u(4)**2);br(i)=merge(1,-1,u(5)<.5_dp)
        tof(i)=(nr(i)+0.3_dp+2.7_dp*u(5))*acos(-1._dp)*s**1.5_dp/sqrt(2._dp)
      endif
    enddo
    best1=huge(1._dp);bestb=huge(1._dp);chk=0
    do rep=1,9
      call system_clock(c0,cr)
      do i=1,n
        call excali_one(r1(:,i),r2(:,i),tof(i),dir(i),v1(:,i),v2(:,i),st(i),hf(i),nr(i),br(i))
      enddo
      call system_clock(c1);best1=min(best1,real(c1-c0,dp)/cr/n*1e9_dp);chk=chk+sum(v1(1,:),mask=st==0)
      call system_clock(c0,cr)
      call excali_batch(n,r1,r2,tof,dir,v1,v2,st,hf,nr,br)
      call system_clock(c1);bestb=min(bestb,real(c1-c0,dp)/cr/n*1e9_dp);chk=chk+sum(v1(1,:),mask=st==0)
    enddo
    print '(a8,a,f8.2,a,f8.2,a,i6,a,es10.2)',lab,' single',best1,' batch',bestb,' solved',count(st>=0.and.st<=2),' chk',chk
  enddo
end program
