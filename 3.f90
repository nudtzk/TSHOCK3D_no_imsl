program main
use vars
use subs
implicit none
integer i,j,k
real::time_begin,time_end
integer::number1,number2,number3



l=0.0
h=0.0
t=0.0
st=0
sn=0

call input_and_init()
call generate_symm_boundary()
!call generate_symm_boundary_cylinder()

do j=1,nj
	do k=1,3
		a(k,j)=(f_ext(k,j)-f_int(k,j))/mass(j)
    end do
end do
open(888,file="test_dp_fun.txt")
write(888,"(3x,'t',3x,'p',3x,'e',3x,'q',3x,'b1_or_x1=',3x,'b2_or_x2',3x,'mod',3x,'condition',3x,'strain_vol',3x,'x1',3x,'x1_1',3x,'p_dst',3x,'p_dsp')")

open(999,file="test_strain.txt")
write(999,"(3x,'t',3x,'p',3x,'plastic_index',3x,'strain_p_eff',3x,'strain_vol')")
call cpu_time(time_begin)
do while (t<t_total)




	 	   
   t=t+dt
   st=st+1
   sn=sn+1
   if(problem_mode==1)then
   call compute_energy_xray()
   end if


   call element_compute()
   if(st==nscr) then
   call out_put()
   call out_put_element()
     st=0
!-----------------------------------------------------
	 print*,t,'us'
	 print*,int(t/t_total*100),'% has been caculated'

number1=0
number2=0
number3=0
do i=1,ne

if(fail(i)==1) then
number1=number1+1
end if

if(gas(i)==1) then
number2=number2+1
end if

if(plastic_index(i)==1) then
number3=number3+1
end if

end do
call cpu_time(time_end)
print*,'Time of operation is',time_end-time_begin
print*,'fail_number=',number1
print*,'gas_number=',number2
print*,'plastic_index=',number3
!-----------------------------------------------------
   end if
   call out_put_energy()
!   call out_put_history()
  if(problem_mode==1) call compute_blowoff()
   end do





end                        
