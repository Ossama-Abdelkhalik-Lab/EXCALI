! Fortran reference for the CUDA port's host test: same rows, same output format.
program fortran_driver
   use excali_solver, only: excali_one
   use iso_fortran_env, only: real64
   implicit none
   integer,parameter::dp=real64
   character(4096)::line
   character(512)::fi,fo
   real(dp)::v(13),r1(3),r2(3),v1(3),v2(3)
   integer::u,uo,ios,st,h
   call get_command_argument(1,fi); call get_command_argument(2,fo)
   open(newunit=u,file=fi,status='old'); open(newunit=uo,file=fo,status='replace')
   read(u,'(a)')line
   do
      read(u,'(a)',iostat=ios)line; if(ios/=0)exit
      v=0; read(line,*,iostat=ios)v
      if(nint(v(11))/=0)cycle
      r1=v(3:5);r2=v(6:8)
      call excali_one(r1,r2,v(9),nint(v(10)),v1,v2,st,h,0,-1)
      write(uo,'(i0,1x,i0,6(1x,z16.16))')nint(v(1)),st,v1,v2
   enddo
end program
