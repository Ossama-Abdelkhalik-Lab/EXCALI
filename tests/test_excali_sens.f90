! excali_sens_solve against 100-digit references (tests/data/sens_reference.txt,
! from mpmath central differences of an independent Lambert solver): the
! 6x7 Jacobian and 6x7x7 Hessian on stress cases, relative Frobenius error per
! block within max(1e-14, 5 x the physical floor), the floor being the change
! of the exact J (H) when the double inputs move by one ulp. Exactly collinear
! endpoints with a plane normal: J only (out-of-plane second derivatives
! depend on convention there).
program test_excali_sens
   use iso_fortran_env, only: real64
   use excali_sens
   implicit none
   integer,parameter::dp=real64
   character(1024)::fn
   character(64)::name
   real(dp)::r1(3),r2(3),tof,mu,pn(3),v1(3),v2(3),jac(6,7),hes(6,7,7),jr(6,7),hr(6,7,7),gj,gh,ej,eh
   integer::u,n,i,d,m,br,pf,st,bad
   call get_command_argument(1,fn)
   open(newunit=u,file=fn,status='old')
   read(u,*)n
   bad=0
   do i=1,n
      read(u,'(a)')name
      read(u,*)r1,r2,tof,mu,d,m,br,pf,pn
      read(u,*)gj,gh
      read(u,*)jr
      read(u,*)hr
      if(pf==1)then
         call excali_sens_solve(r1,r2,tof,mu,d,m,br,v1,v2,jac,st,hes,pn)
      else
         call excali_sens_solve(r1,r2,tof,mu,d,m,br,v1,v2,jac,st,hes)
      endif
      ej=norm2(jac-jr)/norm2(jr);eh=0
      if(gh>0)eh=norm2(hes-hr)/norm2(hr)
      print '(a28,a,i2,a,es9.2,a,es9.2,a,es9.2,a,es9.2,a)',name,' status',st,'  J',ej,' (gate',gj,')  H',eh,' (gate',gh,')'
      if(st<0.or.st>2.or..not.(ej<=gj).or..not.(eh<=gh))bad=bad+1
   enddo
   if(bad>0)error stop 'sensitivities outside their gates'
end program
