program transfer
   use iso_fortran_env, only: real64
   use excali_solver, only: excali_solve
   implicit none
   real(real64)::r1(3),r2(3),v1(3),v2(3),mu,tof
   integer::status
   mu=398600.4418_real64 ! km^3/s^2
   r1=[7000._real64,0._real64,0._real64]
   r2=[0._real64,7000._real64,0._real64]
   tof=acos(-1._real64)/2*sqrt(7000._real64**3/mu)
   call excali_solve(r1,r2,tof,mu,1,0,-1,v1,v2,status)
   print *, 'status:',status
   print *, 'departure velocity (km/s):',v1
   print *, 'arrival velocity (km/s):',v2
   if(status/=0)error stop 'transfer failed'
end program
