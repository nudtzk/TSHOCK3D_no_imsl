module subs
 !  use imsl
   use vars
   use out_of_imsl
contains

subroutine input_and_init()
   implicit none
	
   integer i,j,k,n1,n2,n3,n4,n5,n6,n7,n8
   real(8) lamda
   namelist/compute_control/&
           problem,problem_mode,shape,eos_mode,flyer_length,lmax,lmax_iso,lmax_aniso,length_1,length_2,length_3,dl_1,dl_2,dl_3,    &
		   q1,q2,qhg,dt,t_total,nscr,ne,nj,vf,     &
		   num_elefix_x,num_elefix_y,num_elefix_z,ne_his_x,ne_his_y,ne_his_z,num_time_record,time_record
   
   namelist/layer_setting/mat_model,yield_mode,length_layer

   namelist/model_0/layer_no_iso,rho0_iso,c0_iso,s_iso,gama_iso,      &
                    h0_iso,es_iso,ex_iso,ey_iso,ez_iso,pxy_iso,pxz_iso,pyz_iso,gxy_iso,gxz_iso,gyz_iso,strs_yld_0,strs_frct_0

   namelist/model_0_yield_1/ax_iso,nx_iso,beta0_iso,strain_rate0_iso


   namelist/model_1/layer_no_aniso,rho0_aniso,c0_aniso,s_aniso,gama_aniso,h0_aniso,es_aniso,     &
                    ex_aniso,ey_aniso,ez_aniso,pxy_aniso,pxz_aniso,pyz_aniso,gxy_aniso,gxz_aniso,gyz_aniso,     &
                    strs_yld_110,strs_yld_220,strs_yld_330,strs_yld_120,strs_yld_130,strs_yld_230,       &
					strs_frct_110,strs_frct_220,strs_frct_330,strs_frct_120,strs_frct_130,strs_frct_230
					        

					    
   namelist/model_1_yield_1/ax_aniso,nx_aniso,ay_aniso,ny_aniso,az_aniso,nz_aniso,       &
                            axy_aniso,nxy_aniso,axz_aniso,nxz_aniso,ayz_aniso,nyz_aniso,beta0_aniso,strain_rate0_aniso

   

   namelist/xray_output/cm_r1,cm_r2,cm_theta_0,cm_theta_10,cm_theta_20,cm_theta_30,cm_theta_40,cm_theta_45,    &
                        cm_theta_50,cm_theta_60,cm_theta_70,cm_theta_80,cm_theta_90

   namelist/model_1_paper_1/parameter_A1,parameter_A2,parameter_A3,parameter_B0,parameter_T1,                        &
                             a_11,a_22,a_33,a_12,a_13,a_23,a_44,a_55,a_66,sigma_well_1,sigma_well_10,mean_plastic    
						 
   namelist/model_1_paper_2/tensil_11,tensil_22,tensil_33,tensil_12,tensil_23,tensil_31,                             &
						     fracture_e11,fracture_e22,fracture_e33,fracture_e12,fracture_e23,fracture_e31,           &
						     damage_coupling

   namelist/xray/n_blackbody,t_kev,b_energy,wl_short,wl_long,rad_time_mode,rad_time_width,rad_energy_total,xlayer
            

 



      
! initial globle logical variables fo material model
   mat0=.false.
   mat1=.false.  
   open(100,file="in1.txt")
   open(101,file="ele.dat")
   open(102,file="node.dat")
 

   read(100,nml=compute_control)
   call allocate_mode()
   



   h_dt=0.5*dt
    !在compute_volume中决定顺序
   read(100,nml=layer_setting)
 
   if(lmax>max_layer)   stop "Some input parameters not right, program terminated."

   if(shape==0) then
       length_x=length_1
	   length_y=length_2
	   length_z=length_3
	   dx=dl_1
	   dy=dl_2
	   dz=dl_3
	   lx=length_layer
   else if(shape==1) then
       r_in=length_1
	   r_out=length_2
	   r_height=length_3
	   dr=dl_1
	   d_theta=dl_2
	   d_height=dl_3
	   lr=length_layer
	   
   end if
 
   do i=1,lmax
 
       select case(mat_model(i))
	   case(0) 
	       mat0=.true.
	   case(1)
	       mat1=.true.
		  
	   case default
	       write(*,"('there is not such material model No.',    &
		             i2,'in the layer No.',i2,'.')")    &
					 mat_model(i),i
		   stop "program terminated."
	   end select
   end do
   
   if(mat0) then

       read(100,nml=model_0)

	   do i=1,lmax_iso
	       rho0(layer_no_iso(i))=rho0_iso(i)
		   c0(layer_no_iso(i))=c0_iso(i)
           s(layer_no_iso(i))=s_iso(i)
		   gama(layer_no_iso(i))=gama_iso(i)
           
		   h0(layer_no_iso(i))=h0_iso(i)
		   es(layer_no_iso(i))=es_iso(i)
		   cn(layer_no_iso(i))=c0_iso(i)**2/(gama_iso(i)*es_iso(i))
           
		   ex_0(layer_no_iso(i))=ex_iso(i)
		   ey_0(layer_no_iso(i))=ex_iso(i)
		   ez_0(layer_no_iso(i))=ex_iso(i)
		   pxy(layer_no_iso(i))=pxy_iso(i)
		   pxz(layer_no_iso(i))=pxz_iso(i)
		   pyz(layer_no_iso(i))=pyz_iso(i)

		   gxy_0(layer_no_iso(i))=gxy_iso(i)
		   gxz_0(layer_no_iso(i))=gxz_iso(i)
		   gyz_0(layer_no_iso(i))=gyz_iso(i)

		   strs_yld_0(layer_no_iso(i))=strs_yld_0(i)
           strs_frct_0(layer_no_iso(i))=strs_frct_0(i)

		   write(2,*)layer_no_iso(i),rho0(layer_no_iso(i)), ex_0(layer_no_iso(i)),ey_0(layer_no_iso(i)),ez_0(layer_no_iso(i)),pxy(layer_no_iso(i)),pxz(layer_no_iso(i)),pyz(layer_no_iso(i))
	   end do
       
       if(yield_mode(layer_no_iso(1))==1) then   !相同类型的材料选择相同的屈服模式，因此由各类材料中的第一层即可确定屈服模式
	 
		    read(100,nml=model_0_yield_1)
		    	   
			do i=1,lmax_iso
			    ax(layer_no_iso(i))=ax_iso(i)
		        nx(layer_no_iso(i))=nx_iso(i)
				beta0(layer_no_iso(i))=beta0_iso(i)
		        strain_rate0(layer_no_iso(i))=strain_rate0_iso(i)
			end do
            write(2,"('material layer',i2,'is isotropic and rate-related yield mode is used')")layer_no_iso(1)
	   end if

   end if

   if(mat1) then

       read(100,nml=model_1)
   
	!   read(100,nml=model_1_paper_1)
	!   read(100,nml=model_1_paper_2)


       do i=1,lmax_aniso
	       rho0(layer_no_aniso(i))=rho0_aniso(i)
		   c0(layer_no_aniso(i))=c0_aniso(i)
           s(layer_no_aniso(i))=s_aniso(i)
		   gama(layer_no_aniso(i))=gama_aniso(i)
         
		   h0(layer_no_aniso(i))=h0_aniso(i)
		   es(layer_no_aniso(i))=es_aniso(i)
		   cn(layer_no_aniso(i))=c0_aniso(i)**2/(gama_aniso(i)*es_aniso(i))
           
           ex_0(layer_no_aniso(i))=ex_aniso(i)
           ey_0(layer_no_aniso(i))=ey_aniso(i)
		   ez_0(layer_no_aniso(i))=ez_aniso(i)

		   pxy(layer_no_aniso(i))=pxy_aniso(i)
		   pxz(layer_no_aniso(i))=pxz_aniso(i)
		   pyz(layer_no_aniso(i))=pyz_aniso(i)

		   gxy_0(layer_no_aniso(i))=gxy_aniso(i)
		   gxz_0(layer_no_aniso(i))=gxz_aniso(i)
           gyz_0(layer_no_aniso(i))=gyz_aniso(i)


           strs_yld_110(layer_no_aniso(i))=strs_yld_110(i)
		   strs_yld_220(layer_no_aniso(i))=strs_yld_220(i)
		   strs_yld_330(layer_no_aniso(i))=strs_yld_330(i)
           strs_yld_120(layer_no_aniso(i))=strs_yld_120(i)
		   strs_yld_130(layer_no_aniso(i))=strs_yld_130(i)
		   strs_yld_230(layer_no_aniso(i))=strs_yld_230(i)
 
           strs_frct_110(layer_no_aniso(i))=strs_frct_110(i)
		   strs_frct_220(layer_no_aniso(i))=strs_frct_220(i)
		   strs_frct_330(layer_no_aniso(i))=strs_frct_330(i)
		   strs_frct_120(layer_no_aniso(i))=strs_frct_120(i)
 		   strs_frct_130(layer_no_aniso(i))=strs_frct_130(i)
		   strs_frct_230(layer_no_aniso(i))=strs_frct_230(i)
		          
		write(2,*)layer_no_aniso(i),rho0(layer_no_aniso(i)),ex_0(layer_no_aniso(i)),ey_0(layer_no_aniso(i)),ez_0(layer_no_aniso(i)),pxy(layer_no_aniso(i)),pxz(layer_no_aniso(i)),pyz(layer_no_aniso(i))
		       
        end do
 
		if(yield_mode(layer_no_aniso(1))==1) then

       read(100,nml=model_1_yield_1)
        


		    do i=1,lmax_aniso
			    ax(layer_no_aniso(i))=ax_aniso(i)
		        nx(layer_no_aniso(i))=nx_aniso(i)
		        ay(layer_no_aniso(i))=ay_aniso(i)
		        ny(layer_no_aniso(i))=ny_aniso(i)
		        az(layer_no_aniso(i))=az_aniso(i)	   
		        nz(layer_no_aniso(i))=nz_aniso(i)
                axy(layer_no_aniso(i))=axy_aniso(i)
			    nxy(layer_no_aniso(i))=nxy_aniso(i)
                axz(layer_no_aniso(i))=axz_aniso(i)
			    nxz(layer_no_aniso(i))=nxz_aniso(i)            
                ayz(layer_no_aniso(i))=ayz_aniso(i)
			    nyz(layer_no_aniso(i))=nyz_aniso(i)						
				beta0(layer_no_aniso(i))=beta0_aniso(i)
                strain_rate0(layer_no_aniso(i))=strain_rate0_aniso(i)
				write(2,"('material layer',i2,1x,'is anisotropic and rate-related yield mode is used')")layer_no_aniso(i)
		    end do
			

		end if

 end if

   read(101,*)ia
   read(102,*)xyz0
   
   do i=1,ne
       ne_layer(i)=ia(9,i)   !fortran里读入数据，先在数组第一个位置。结合前处理文件，是按行扫（1-8，i）
   end do
 
   if(shape==0) then

   else if(shape==1) then

  !     read(100,nml=xray_output)
	   
   end if

   if (problem_mode==0) then  

   else if(problem_mode==1) then
  
       read(100,nml=xray)

   end if


   open(900,file="peak_axis_time.txt")
   write(900,"(4x,'t=',10x,'i=',8x,'peak_axis_x=',8x,'peak_axis_y=',7x,'peak_axis_z=',7x,'peak_axis_stress_xx_time=',7x)")
   open(950,file="peak_axis_line.txt")
   write(950,"(4x,'i=',5x,'peak_axis_x=',8x,'peak_axis_y=',8x,'peak_axis_stress_xx_line=',7x,'t=',10x)")
   open(1000,file="strs_peak.txt")
   write(1000,"(4x,'t=',10x,'i=',8x,'peak_x=',8x,'peak_y=',7x,'peak_z=',7x,'peak_stress_xx=',7x)")
   open(1020,file="blowoff.txt")
   write(1020,"('t=',10x,'blowoff=')") 

   open(901,file="failure_condtion.txt")
   write(901,"('t=',10x,'failure_number=',10x,'failure_rate')")


!--------------------------------------------------------------------------------------------------------------
!  初始化节点坐标、位移、速度、加速度、内力分布、外力分布，单元应力分布
   do i=1,ne
       fail(i)=0
   end do
      



   do j=1,nj
       xyz(1,j)=xyz0(1,j)
       xyz(2,j)=xyz0(2,j)
	   xyz(3,j)=xyz0(3,j)

   end do


   do i=1,ne
       n1=ia(1,i)
       n2=ia(2,i)
       n3=ia(3,i)
       n4=ia(4,i)
	   n5=ia(5,i)
       n6=ia(6,i)
       n7=ia(7,i)
       n8=ia(8,i)
       rho(i)=rho0(ne_layer(i))
       call volume_caculate(i)
       volume0(i)=volume(i)          ! 得到初始体积

   end do 


   do j=1,nj
       u(:,j)=0.0                    !节点位移
       a(:,j)=0.0                    !节点加速度
       f_ext(:,j)=0.0
       f_int(:,j)=0.0
       f_hg(:,j)=0.0
   end do 


   do i=1,ne
       gas(i)=0
       e(i)=0.0
       p(i)=0.0
       q(i)=0.0
       strain_p_rate_eff(i)=0.0
	   rate(i)=1.0
	   plastic_index(i)=0
   end do


   do i=1,ne


       strain_xyz(:,:,i)=0.0
       stress_xyz(:,:,i)=0.0
       stress_d_123(:,:,i)=0.0
	   stress_123(:,:,i)=0.0

	if(mat_model(ia(9,i))==1) then
	   strs_yld_123(1,1,i)=strs_yld_110(ne_layer(i))
       strs_yld_123(2,2,i)=strs_yld_220(ne_layer(i))
       strs_yld_123(3,3,i)=strs_yld_330(ne_layer(i))
       strs_yld_123(1,2,i)=strs_yld_120(ne_layer(i))
	   strs_yld_123(1,3,i)=strs_yld_130(ne_layer(i))
	   strs_yld_123(2,3,i)=strs_yld_230(ne_layer(i))

	   strs_frct_123(1,1,i)=strs_frct_110(ne_layer(i))
       strs_frct_123(2,2,i)=strs_frct_220(ne_layer(i))
	   strs_frct_123(3,3,i)=strs_frct_330(ne_layer(i))
	   strs_frct_123(1,2,i)=strs_frct_120(ne_layer(i))
       strs_frct_123(1,3,i)=strs_frct_130(ne_layer(i))
	   strs_frct_123(2,3,i)=strs_frct_230(ne_layer(i))
	else if(mat_model(ia(9,i))==0) then
	   strs_yld(i)=strs_yld_0(ne_layer(i))
	   strs_frct(i)=strs_frct_0(ne_layer(i))
	end if


      !是否需要矩阵化还要看看

   end do
   
   if(problem_mode==0) then
       do j=1,nj
	       v(:,j)=0.0
	   end do
       do j=1,nj
	       if(xyz0(1,j)<(flyer_length-dx/2))then
			v(1,j)=vf
			else if(xyz0(1,j)>(flyer_length-dx/2) .and. xyz0(1,j)<(flyer_length+dx/2))then
			v(1,j)=vf/2
			end if
	   end do

!完全是针对陈主席的k文件格式修改的速度初始方案

   else if(problem_mode==1) then

   end if



   shape_n(1)=0.125*(1-l)*(1+h)*(1+m)
   shape_n(2)=0.125*(1-l)*(1-h)*(1+m)
   shape_n(3)=0.125*(1-l)*(1-h)*(1-m)
   shape_n(4)=0.125*(1-l)*(1+h)*(1-m)
   shape_n(5)=0.125*(1+l)*(1+h)*(1+m)
   shape_n(6)=0.125*(1+l)*(1-h)*(1+m)
   shape_n(7)=0.125*(1+l)*(1-h)*(1-m)
   shape_n(8)=0.125*(1+l)*(1+h)*(1-m)

   shape_n0(1)=0.125
   shape_n0(2)=0.125
   shape_n0(3)=0.125
   shape_n0(4)=0.125
   shape_n0(5)=0.125
   shape_n0(6)=0.125
   shape_n0(7)=0.125
   shape_n0(8)=0.125

   dn0_l(1)=-0.125
   dn0_l(2)=-0.125
   dn0_l(3)=-0.125
   dn0_l(4)=-0.125
   dn0_l(5)=0.125
   dn0_l(6)=0.125
   dn0_l(7)=0.125
   dn0_l(8)=0.125

   dn0_h(1)=0.125
   dn0_h(2)=-0.125
   dn0_h(3)=-0.125
   dn0_h(4)=0.125
   dn0_h(5)=0.125
   dn0_h(6)=-0.125
   dn0_h(7)=-0.125
   dn0_h(8)=0.125

   dn0_m(1)=0.125
   dn0_m(2)=0.125
   dn0_m(3)=-0.125
   dn0_m(4)=-0.125
   dn0_m(5)=0.125
   dn0_m(6)=0.125
   dn0_m(7)=-0.125
   dn0_m(8)=-0.125

   do i=1,ne
       x0(i)=0.0
       y0(i)=0.0
	   z0(i)=0.0
       do j=1,8
           k=ia(j,i)                          !第i个单元第j个节点对应的整体节点号
           x0(i)=x0(i)+shape_n0(j)*xyz0(1,k)   !第i个单元任一点的x坐标,当l,h,m为(0,0,0)时,表示每单元中心的初始整体坐标
           y0(i)=y0(i)+shape_n0(j)*xyz0(2,k)
		   z0(i)=z0(i)+shape_n0(j)*xyz0(3,k)		  
       end do
!环形问题旋转计算初始角度，这个留在以后解决
	  
   end do
!回去查阅下输出格式问题
 
 
  n_cmup=0
  if(shape==0)then
	 do i=1,ne
		if(y0(i)<0.9*dy .and. z0(i)<0.9*dz) then
		  n_cmup=n_cmup+1
		 cmup(n_cmup)=i
		 end if
	 end do
	 print*,'axis_number is',n_cmup
ne_his_x(1)=cmup(int(n_cmup/8*1)) 
ne_his_x(2)=cmup(int(n_cmup/8*2)) 
ne_his_x(3)=cmup(int(n_cmup/8*3)) 
ne_his_x(4)=cmup(int(n_cmup/8*4)) 
ne_his_x(5)=cmup(int(n_cmup/8*5)) 
ne_his_x(6)=cmup(int(n_cmup/8*6)) 
ne_his_x(7)=cmup(int(n_cmup/8*7)) 
ne_his_x(8)=cmup(int(n_cmup/8*8)) 
ne_his_x(9)=cmup(n_cmup)
  else
    n_cmup=0
  end if

call output_filename_definition()


!计算对角化质量矩阵,采用高斯单点积分法 在(l,h,m)=(0,0,0)处积分,每个节点分担1/8质量
   do k=1,ne
      do j=1,8
          i=ia(j,k)
		  mass(i)=mass(i)+rho0(ne_layer(k))*0.125*volume0(k)	!节点质量 
      end do
   end do !实际上分了三个方向的m一维向量有点多余的，集中质量在节点上。

   open(190,file="mass.txt")
   write(190,*) "m(i)="
   write(190,10) (mass(i),i=1,nj)
   10 format (6(f12.8),3x)

   do i=1,lmax
       if(mat_model(i)==0)  then                                   !完全可以统一，各向同性是各向异性特殊情况
	       lamda=1-pxy(i)-2*pxy(i)*pxy(i)                          !除了多了三维剪切这一部分是一样的
           c123_0(1,1,i)=ex_0(i)*(1-pxy(i))/lamda                  !层属性
		   c123_0(2,2,i)=c11_0(i)
		   c123_0(3,3,i)=c11_0(i)
		   c123_0(1,2,i)=ex_0(i)*pxy(i)/lamda
		   c123_0(1,3,i)=c12_0(i)
		   c123_0(2,3,i)=c12_0(i)
		   c123_0(4,4,i)=gxy_0(i)
           c123_0(5,5,i)=gxz_0(i)
		   c123_0(6,6,i)=gyz_0(i)

	   else if(mat_model(i)==1) then

	       lamda=1-pxy(i)*pxy(i)*ey_0(i)/ex_0(i)-pxz(i)*pxz(i)*ez_0(i)/ex_0(i)        &
	             -pyz(i)*pyz(i)*ez_0(i)/ey_0(i)-2*pxy(i)*pyz(i)*pxz(i)*ez_0(i)/ex_0(i)


		   c123_0(1,1,i)=ex_0(i)*(1-pyz(i)*pyz(i)*ez_0(i)/ey_0(i))/lamda    !层属性
           c123_0(2,2,i)=ey_0(i)*(1-pxz(i)*pxz(i)*ez_0(i)/ex_0(i))/lamda
           c123_0(3,3,i)=ez_0(i)*(1-pxy(i)*pxy(i)*ey_0(i)/ex_0(i))/lamda
           c123_0(1,2,i)=ey_0(i)*(pxy(i)+pyz(i)*pxz(i)*ez_0(i)/ey_0(i))/lamda
           c123_0(1,3,i)=ez_0(i)*(pxz(i)+pxy(i)*pyz(i))/lamda
           c123_0(2,3,i)=ez_0(i)*(pyz(i)+pxy(i)*pxz(i)*ey_0(i)/ez_0(i))/lamda
           c123_0(4,4,i)=gxy_0(i)
		   c123_0(5,5,i)=gxz_0(i)
		   c123_0(6,6,i)=gyz_0(i)
           
	   end if

	   open(120,file="c_compute.txt")
           write(120,"('c11_0=',f10.6,'c22_0=',f10.6,'c33_0=',f10.6,'c12_0=',f10.6,'c13_0=',f10.6,'c23_0=',f10.6,'c44_0=',f10.6)")    &
		              c11_0(i),c22_0(i),c33_0(i),c12_0(i),c13_0(i),c23_0(i),c44_0(i)


       
   end do


   if(problem_mode==1)  then
   
   call xray_init()
   end if

         
   do i=1,num_elefix_x
       open(600+i,file=fne_his_x(i))
       write(600+i,"('t=',10x,'x=',10x,'y=',10x,'z=',10x,'v_cen_x=',4x,'v_cen_y=',4x,'v_cen_z=',4x,'strain_xx=',2x,'strain_yy=',2x,'strain_zz=',2x,'strain_xy=',2x,'strain_xz=',2x,'strain_yz=',2x,        &
	                 'stress_xx=',2x,'stress_yy=',2x,'stress_zz=',2x,'stress_xy=',2x,'stress_xz=',2x,'stress_yz=',2x,'p=',10x,'e=')")  
					   
   end do

   do i=1,num_elefix_y
       j=ne_his_y(i)
	   open(700+i,file=fne_his_y(i))
       write(700+i,"('t=',10x,'x=',10x,'y=',10x,'z=',10x,'v_cen_x=',4x,'v_cen_y=',4x,'v_cen_z=',4x,'strain_xx=',2x,'strain_yy=',2x,'strain_zz=',2x,'strain_xy=',2x,'strain_xz=',2x,'strain_yz=',2x,        &
	                 'stress_xx=',2x,'stress_yy=',2x,'stress_zz=',2x,'stress_xy=',2x,'stress_xz=',2x,'stress_yz=',2x,'p=',10x,'e=')")  
   end do

   do i=1,num_elefix_z
       j=ne_his_z(i)
		open(800+i,file=fne_his_z(i))
        write(800+i,"('t=',10x,'x=',10x,'y=',10x,'z=',10x,'v_cen_x=',4x,'v_cen_y=',4x,'v_cen_z=',4x,'strain_xx=',2x,'strain_yy=',2x,'strain_zz=',2x,'strain_xy=',2x,'strain_xz=',2x,'strain_yz=',2x,        &
	                 'stress_xx=',2x,'stress_yy=',2x,'stress_zz=',2x,'stress_xy=',2x,'stress_xz=',2x,'stress_yz=',2x,'p=',10x,'e=')")  
  end do




  
   return 
end subroutine input_and_init

subroutine generate_symm_boundary()
implicit none
integer::i,j

integer::nj_xz
nj_xz=0
do j=1,nj
if(xyz0(2,j)<0.14/2)then
 nj_xz=nj_xz+1
 symm_xz(nj_xz)=j
end if

end do
end subroutine generate_symm_boundary

!--------------------------------------------------------------------------------------------------------------------
subroutine compute_c(i)
implicit none
integer i
real(8)::D(6,6)
real(8)::D_inverse(6,6)
real(8) lamda   

   !刚度系数  
   if(mat_model(ne_layer(i))==0)  then                          !实际上可以统一
	   if(yield_mode(ne_layer(i))==0) then
		   ex(i)=ex_0(ne_layer(i))
           ey(i)=ex(i)
		   ez(i)=ex(i)
		   gxy(i)=gxy_0(ne_layer(i))
		   gxz(i)=gxy(i)
           gyz(i)=gxy(i)
		   pxy(i)=pxy(ne_layer(i))
		   pxz(i)=pxy(i)
           pyz(i)=pxy(i)
		   
            else if(yield_mode(ne_layer(i))==1) then
		   ex(i)=ex_0(ne_layer(i))
		   ey(i)=ex(i)
		   ez(i)=ex(i)
		   gxy(i)=gxy_0(ne_layer(i))
		   gxz(i)=gxy(i)
           gyz(i)=gxy(i)
		   pxy(i)=pxy(ne_layer(i))
		   pxz(i)=pxy(i)
           pyz(i)=pxy(i)		  
             end if      
   	   
   else if(mat_model(ne_layer(i))==1) then
	   if(yield_mode(ne_layer(i))==0) then
		   ex(i)=ex_0(ne_layer(i))
	       ey(i)=ey_0(ne_layer(i))
		   ez(i)=ez_0(ne_layer(i))
		   gxy(i)=gxy_0(ne_layer(i))
		   gxz(i)=gxz_0(ne_layer(i))
		   gyz(i)=gyz_0(ne_layer(i))
           pxy(i)=pxy(ne_layer(i))
		   pxz(i)=pxz(ne_layer(i))
		   pyz(i)=pyz(ne_layer(i))
		   
	   else if(yield_mode(ne_layer(i))==1) then
		   ex(i)=ex_0(ne_layer(i))
		   ey(i)=ey_0(ne_layer(i))
		   ez(i)=ez_0(ne_layer(i))
		   gxy(i)=gxy_0(ne_layer(i))
		   gxz(i)=gxz_0(ne_layer(i))
		   gyz(i)=gyz_0(ne_layer(i))
		   pxy(i)=pxy(ne_layer(i))
		   pxz(i)=pxz(ne_layer(i))
		   pyz(i)=pyz(ne_layer(i))

           end if  
    end if
	D=0
    D(1,1)=1/ex(i)
	D(1,2)=-pxy(i)/ex(i)
	D(1,3)=-pxz(i)/ex(i)
	D(2,1)=D(1,2)
	D(2,2)=1/ey(i)
	D(2,3)=-pyz(i)/ey(i)
	D(3,1)=D(1,3)
	D(3,2)=D(2,3)
	D(3,3)=1/ez(i)
	D(4,4)=1/gxy(i)
	D(5,5)=1/gxz(i)
	D(6,6)=1/gyz(i)
	
	call solve_inverse(D,D_inverse,6)
	c123(:,:,i)=D_inverse

    c11(i)=c123(1,1,i)
	c12(i)=c123(1,2,i)
	c13(i)=c123(1,3,i)
	c21(i)=c123(2,1,i)
	c22(i)=c123(2,2,i)
	c23(i)=c123(2,3,i)
	c31(i)=c123(3,1,i)
	c32(i)=c123(3,2,i)
	c33(i)=c123(3,3,i)
	c44(i)=c123(4,4,i)
	c55(i)=c123(5,5,i)
	c66(i)=c123(6,6,i)




end subroutine compute_c    !c矩阵6*6
!--------------------------------------------------------------------------------------------------------------

subroutine velocity_compute()
   implicit none 
   integer(4) i,j,k,n1,n2,n3,n4,n5,n6,n7,n8

   do j=1,nj

       v(1,j)=v(1,j)+h_dt*a(1,j)  !n+1/2时刻的速度
       v(2,j)=v(2,j)+h_dt*a(2,j)
       v(3,j)=v(3,j)+h_dt*a(3,j)

       v_mid(1,j)=v(1,j)        !每个节点中间时刻的速度（即第一次部分更新节点速度）
       v_mid(2,j)=v(2,j) 
       v_mid(3,j)=v(3,j) 
   
       du(1,j)=dt*v_mid(1,j)        !位移增量
       du(2,j)=dt*v_mid(2,j)
       du(3,j)=dt*v_mid(3,j)
       
       u(1,j)=u(1,j)+dt*v_mid(1,j)    !用中间时刻的速度*dt计算位移,匀加速运动的中点速度为全程平均速度
       u(2,j)=u(2,j)+dt*v_mid(2,j) 
       u(3,j)=u(3,j)+dt*v_mid(3,j)
   end do  

   do j=1,nj                !求各节点更新后的坐标
       xyz(1,j)=xyz(1,j)+du(1,j)
       xyz(2,j)=xyz(2,j)+du(2,j) 
	   xyz(3,j)=xyz(3,j)+du(3,j)    
   end do 

   do i=1,ne
       if (fail(i)==0) then
           x(i)=0.0
           y(i)=0.0
		   z(i)=0.0
		  




           do j=1,8
               k=ia(j,i)                      !第i个单元第j个节点对应的整体节点号
               x(i)=x(i)+shape_n(j)*xyz(1,k)   !第i个单元任一点的x坐标,当l,h,m为(0,0,0)时,表示每单元中心的现时整体坐标,问题是l h m什么时候出现了？？？ans,系统默认为0
               y(i)=y(i)+shape_n(j)*xyz(2,k)
			   z(i)=z(i)+shape_n(j)*xyz(3,k) 
           end do

       end if
	
	     	       
   end do


!计算单元体积和密度及比体积

   do i=1,ne
       if (fail(i)==0) then
           n1=ia(1,i)
           n2=ia(2,i)
           n3=ia(3,i)
           n4=ia(4,i)
           n5=ia(5,i)
           n6=ia(6,i)
           n7=ia(7,i)
           n8=ia(8,i)


           volume_old(i)=volume(i)    !中间过渡量

           call volume_caculate(i)
		 
           dvolume(i)=volume(i)-volume_old(i)                             		                                                                 
           volume_average(i)=0.5*(volume(i)+volume_old(i))
           rho_old(i)=rho(i)
           rho(i)=rho0(ne_layer(i))*volume0(i)/volume(i)
		   vol0(i)=1.0/rho0(ne_layer(i))
           vol_old(i)=1.0/rho_old(i)     !计算比体积
           vol(i)=1.0/rho(i)
           dvol(i)=vol(i)-vol_old(i)
           vol_average(i)=(vol(i)+vol_old(i))/2.0
       end if
   end do


!计算人工粘性项,改为三维是否存在问题？？？？？注意人工体积粘性只针对压缩状态
   do i=1,ne
       if (fail(i)==0) then
           q_old(i)=q(i) 
           if (dvol(i)>=0) then
               q(i)=0.0
           else
               q(i)=q1*q1*volume(i)*rho(i)*(dvol(i)/(vol(i)*dt))**2       &
                  +q2*sqrt(volume(i))*c0(ne_layer(i))*rho(i)*abs(dvol(i)/(vol(i)*dt))  !黄霞版本粘性力Von-Neumann
       !       q(i)=q1*rho(i)*volume(i)**(2/3)*(dvol(i)/(vol(i)*dt))**2+q2*rho(i)*volume(i)**(1/3)*c0(ne_layer(i))*abs(dvol(i)/(vol(i)*dt))

           end if		
           q_av(i)=0.5*(q_old(i)+q(i))
	   end if
   end do

end subroutine velocity_compute


!计算更新后的单元的应变、应力
subroutine  element_compute()  
                                  !local属性为element_compute带出来的
   implicit none
   integer i,j,k
   integer::num_bound_xy,num_bound_xz
   integer::n1,n2,n3,n4,n5,n6,n7,n8
   integer::np

   real(8) velocity_xyz_local(3,8)
   real(8) dstrain_xyz_local(3,3)
   real(8) strain_xyz_local(3,3)
   real(8) strain_xyz_rate_local(3,3)
   real(8) dstrain_123_local(3,3)
   real(8) strain_123_local(3,3)
   real(8) dstrain_d_123_local(3,3)
   real(8) stress_123_old_local(3,3)
   real(8) stress_d_123_old(3,3)
   real(8) dstress_d_123(3,3)
   real(8) stress_d_123_old_local(3,3)
   real(8) strain_123_rate_local(3,3)

   real(8) try_stress_123_local(3,3)
   real(8) try_p
   real(8) try_sd_123(3,3)
   real(8) e_old,p_old
  ! real(8) strain_123_old_local(3,3)
   real(8) f1_mises,f_mises
   real(8) xx
   real(8) f_hill,f1_hill

   real::temp
   integer::temp_index
     print*,1,t
   call velocity_compute()  
   do i=1,ne
	

   if(i==ne_his_x(6))then
   
  
   temp=0.0  
   do j=1,3
	do k=1,3
		temp=temp+strain_p_123(j,k,i)*strain_p_123(j,k,i)
	end do
   end do    
	temp=sqrt(temp)
  
   open(999,file="test_strain.txt")
   write(999,99)t,p(i),plastic_index(i),temp,(strain_123(1,1,i)+strain_123(2,2,i)+strain_123(3,3,i))
    
  
   99 format (2(f8.4,3x),i2,3x,7(f8.6,3x))





   	
   end if
	   if (fail(i)==0) then

           call jacb(i)
 
!           strain_xx(i)=(dn_x(1)*u(2*n1-1)+dn_x(2)*u(2*n2-1)+dn_x(3)*u(2*n3-1)+dn_x(4)*u(2*n4-1))   !以拉为正
!           strain_yy(i)=(dn_y(1)*u(2*n1)+dn_y(2)*u(2*n2)+dn_y(3)*u(2*n3)+dn_y(4)*u(2*n4))
!           strain_xy(i)=(0.5*(dn_x(1)*u(2*n1)+dn_x(2)*u(2*n2)+dn_x(3)*u(2*n3)+dn_x(4)*u(2*n4)    &
!                             +dn_y(1)*u(2*n1-1)+dn_y(2)*u(2*n2-1)+dn_y(3)*u(2*n3-1)+dn_y(4)*u(2*n4-1)))
		do j=1,8      !初始化velocity_xyz(3,8)matrix
          velocity_xyz_local(1,j)=v(1,ia(j,i))
		  velocity_xyz_local(2,j)=v(2,ia(j,i))
		  velocity_xyz_local(3,j)=v(3,ia(j,i))
       end do
	                                                                    ! velocity_xyz(:,:,i)=velocity_xyz_local(:,:)
       
           dstrain_xyz_local=matmul(velocity_xyz_local,dn_xyz)*dt
		   dstrain_xyz_local=(dstrain_xyz_local+transpose(dstrain_xyz_local))/2

	                    
		   strain_xyz_rate_local=dstrain_xyz_local/dt
           strain_xyz_rate(:,:,i)=strain_xyz_rate_local(:,:)

           strain_xyz(:,:,i)=strain_xyz(:,:,i)+dstrain_xyz_local(:,:)
		   strain_xyz_local(:,:)=strain_xyz(:,:,i)                                     !更新i号单元在整个体系中的应变，meanwhile，将其带入计算器中


           call rotate_matrix(i)			!需要坐标，在velocity_compute中更新了
		
           if(shape==0) then
				rotate_R=rotate_R

	       else if(shape==1) then  !需要不同的初始角度
	          
			   call original_rotation(i)			
			   rotate_R=matmul(rotate_R,rotate_original) !非正六面体需要这个，文献里有如何rotate的，三个角度罢了,留以后来解决。
	
		   end if

		 !  strain_123_old_local(:,:)=strain_123(:,:,i)    !计算开始前保存下上一步计算得到的123应变
	
		   !12=Rt*xy*R

		   strain_123_rate_local=matmul(matmul(transpose(rotate_R),strain_xyz_rate_local),rotate_R)


           ! !其实说很多余很多余先隐藏strain_123_rate(:,:,i)=strain_123_rate_local(:,:)   
           !由xx旋转到主轴系11，用矩阵乘法，RT*strain_rate*R即可，不用表示为如此复杂的形式，而且黄的做法体现出了很强的理论基础，但真的略显啰嗦
        
           
		   dstrain_123_local=strain_123_rate_local*dt

           strain_123(:,:,i)=strain_123(:,:,i)+dstrain_123_local(:,:)
           strain_123_local=strain_123(:,:,i)
		 
           !这个可以作为专门的来求解球偏分解用的
		   diag_vector=0
		   diag_sum=0
		   diag_matrix=0

		   call diag_extract(dstrain_123_local,diag_vector,3)
		   do j=1,3
				diag_sum=diag_sum+diag_vector(j)
		   end do
				diag_sum=diag_sum/3
				diag_vector=diag_sum
		   call	diag_enchase(diag_matrix,diag_vector,3)
           dstrain_d_123_local=dstrain_123_local-diag_matrix


!           dstrain_d_11=2.0*(vol(i)-vol_old(i))/vol(i)/3.0
!		    dstrain_d_22=-(vol(i)-vol_old(i))/vol(i)/3.0
!           dstrain_d_33=-(vol(i)-vol_old(i))/vol(i)/3.0
!           dstrain_d_12=dstrain_12


!          
            strain_rate_eff(i)=0
         do j=1,3
			do k=1,3
				strain_rate_eff(i)=strain_rate_eff(i)+dstrain_d_123_local(j,k)**2
			end do
		  end do
             	strain_rate_eff(i)=sqrt(2.0*strain_rate_eff(i)/3.0)/dt !等效应变率，利用偏斜表示为根号下2/3偏斜矩数值的平方的遍历,dt除在最后了。
	   
		   e_old=e(i)
           p_old=p(i)                                                  !记录下上一步的压力和能量

		   	
		   stress_123_old_local=stress_123(:,:,i)

		   diag_vector=0
		   diag_sum=0
		   diag_matrix=0
           call diag_extract(stress_123(:,:,i),diag_vector,3)
		   do j=1,3
				diag_sum=diag_sum+diag_vector(j)
		   end do
				diag_sum=diag_sum/3
				diag_vector=diag_sum
		   call	diag_enchase(diag_matrix,diag_vector,3) 
          stress_d_123(:,:,i)=stress_123(:,:,i)-diag_matrix


		   stress_d_123_old_local=stress_d_123(:,:,i)             !记录下上一步的123下的应力和应力偏量张量在塑性计算中
         	
           if(gas(i)==1) then                                     !已经汽化怎么说

               if(rho(i)>=rho0(ne_layer(i))) then
	               call p_eos_gruneisen(i,e_old,p_old)
	           else
	               call p_eos_puff(i,e_old,p_old)
	           end if


			   stress_123(:,:,i)=0
			   stress_123(1,1,i)=-p(i)
			   stress_123(2,2,i)=-p(i)
			   stress_123(3,3,i)=-p(i)	

               stress_xyz(:,:,i)=stress_123(:,:,i)    !无轴旋转问题
		       energy_xray(i)=0.0 
            else
		       if(yield_mode(ne_layer(i))==1) then	!在做之前算好了等效应变率	       
           !弹性试算
		           if (strain_rate_eff(i)<=strain_rate0(ne_layer(i)))  then
                       rate(i)=1.0
                   else
                       rate(i)=1.0+beta0(ne_layer(i))*dlog(strain_rate_eff(i)/strain_rate0(ne_layer(i)))!这不是应变率，而是率相关里那个率因子，涉及到strain-rate撒。目前看rate是不需要修改的            

				   end if 
			   else 
			       rate(i)=1.0 !都是率形式，只不顾率无关率为1.0
			   end if

			   call compute_c(i) !这个没问题了

!第一步导入的应变完全没问题

       call elastic_try(i,dstrain_123_local,dstrain_d_123_local ,        &
                                try_p,try_stress_123_local,&
                                stress_d_123_old_local,  &
							    e_old,p_old)

   
               if(yield_mode(ne_layer(i))==0) then

		           if(mat_model(ne_layer(i))==0) then !各向同性的mises准则
                 
			          call f1_mises_yield_unrelated(i,try_stress_123_local,      &
			                                           f_mises,f1_mises)
				
				     
                       if(f1_mises>0.0) then
					       plastic_index(i)=1					   	
				           call plastic_compute_unrelated_mises(i,stress_123_old_local,      & 
				                                                  dstrain_123_local,dstrain_d_123_local  ,      &
									                              stress_d_123_old_local,    &
												                  e_old,p_old)            
				       else  
						   stress_123(:,:,i)=try_stress_123_local
                           p(i)=try_p

	                   end if
!-----------------------------------------------------------------------------------------------------------
			       else	if(mat_model(ne_layer(i))==1) then
			   
                       call f1_hill_yield_unrelated(i,xx,try_stress_123_local,      &
			                                        f_hill,f1_hill)
	                   if (f1_hill>0.0) then
					       plastic_index(i)=1

  




                           call plastic_compute_unrelated_hill(i,stress_123_old_local,      & !用事先存下的而不是try
				                                                 dstrain_123_local,dstrain_d_123_local  ,      &
									                             stress_d_123_old_local,    &
												                 e_old,p_old)
	                   else
                           
							stress_123(:,:,i)=try_stress_123_local
                            p(i)=try_p
	                   end if

		           end if
				   		   	     
               else if (yield_mode(ne_layer(i))==1) then

                   if(mat_model(ne_layer(i))==0) then

			           call f1_mises_yield_related(i,try_stress_123_local,      &
			                                       f_mises,f1_mises)
                       if(f1_mises>0.0) then
					       plastic_index(i)=1
				           call plastic_compute_related_mises(i,stress_123_old_local,      & !用事先存下的而不是try
				                                                 dstrain_123_local,dstrain_d_123_local  ,      &
									                             stress_d_123_old_local,    &
												                 e_old,p_old)

					       strs_yld(i)=strs_yld_0(ne_layer(i))*(1+ax(ne_layer(i))*strain_p_eff(i)**nx(ne_layer(i)))
				       else 
				          
							stress_123(:,:,i)=try_stress_123_local
                            p(i)=try_p
	                   end if

                   else	if(mat_model(ne_layer(i))==1) then
                           strs_yld_123(1,1,i)=strs_yld_110(ne_layer(i))*rate(i)*(1+ax(ne_layer(i))*strain_p_eff(i)**nx(ne_layer(i)))*rate(i)
                           strs_yld_123(2,2,i)=strs_yld_220(ne_layer(i))*rate(i)*(1+ay(ne_layer(i))*strain_p_eff(i)**ny(ne_layer(i)))*rate(i)
                           strs_yld_123(3,3,i)=strs_yld_330(ne_layer(i))*rate(i)*(1+az(ne_layer(i))*strain_p_eff(i)**nz(ne_layer(i)))*rate(i)
	                       strs_yld_123(1,2,i)=strs_yld_120(ne_layer(i))*rate(i)*(1+axy(ne_layer(i))*strain_p_eff(i)**nxy(ne_layer(i)))*rate(i)
						   strs_yld_123(1,3,i)=strs_yld_130(ne_layer(i))*rate(i)*(1+axz(ne_layer(i))*strain_p_eff(i)**nxz(ne_layer(i)))*rate(i)
					       strs_yld_123(2,3,i)=strs_yld_230(ne_layer(i))*rate(i)*(1+ayz(ne_layer(i))*strain_p_eff(i)**nyz(ne_layer(i)))*rate(i)
				 
                       call f1_hill_yield_related(i,xx,try_stress_123_local,      &
			                                      f_hill,f1_hill) 

                       !write(*,*) f1_hill
                       if (f1_hill>0.0) then
					       plastic_index(i)=1

if(i==ne_his_x(6))then
  temp_index=plastic_index(i)+temp_index
  if(temp_index==1 .and. plastic_index(i)==1)then
   open(998,file="limit.txt")
   write(998,*),(strain_123(1,1,i)+strain_123(2,2,i)+strain_123(3,3,i))
  end if  
 
   end if





                           call plastic_compute_related_hill(i,stress_123_old_local,      & !用事先存下的而不是try
				                                                 dstrain_123_local,dstrain_d_123_local  ,      &
									                             stress_d_123_old_local,    &
												                 e_old,p_old)

   

                       else
                    
							stress_123(:,:,i)=try_stress_123_local
                            p(i)=try_p
                       end if

		           end if

               end if

               !旋转回来，这里同样使用imsl比较方便

               stress_xyz(:,:,i)=matmul(matmul(rotate_R,stress_123(:,:,i)),transpose(rotate_R))
               p(i)=p(i)

               energy_xray(i)=0.0

		       if(e(i)>=es(ne_layer(i)))  gas(i)=1   !判定材料是否汽化
			   if(t>rad_time_width) then

                      call inner_element_frct(i)  !判定材料单元是否失效	
				   !	  call inner_element_frct_paper(i)		  
			   end if

           end if


       else    !失效的话

		   stress_xyz(:,:,i)=0
		   p(i)=0.0
		   e(i)=e(i)+energy_xray(i)
		   
		   energy_xray(i)=0.0

       end if

   end do   
                                                                                                                    

   do j=1,nj
	   do k=1,3	
       f_int_old(k,j)=f_int(k,j)
       f_ext_old(k,j)=f_ext(k,j)
	    end do
   end do
  
   call f_int_total() 
   call f_hourglass()                     
   call f_ext_total()    
   do j=1,nj
		do k=1,3
		if(problem_mode==0)then
	   a(k,j)=(f_ext(k,j)-f_int(k,j)-f_hg(k,j))/mass(j)
	   else if(problem_mode==1)then    
	   a(k,j)=(f_ext(k,j)-f_int(k,j)-f_hg(k,j))/mass(j)
	   end if	       
       v(k,j)=v_mid(k,j)+h_dt*a(k,j)         !再次更新节点速度,n+1时刻的速度
	   
!	   if(any(symm_xy==j))then
!	   v(3,j)=0.0
!	   end if
	   if(any(symm_xz==j))then
	   v(2,j)=0.0
	   end if
	   v(3,j)=0.0


	    end do	
   end do
!采用中心差分法变步长的中心差分法                          
     	       

   do i=1,ne 
       if (fail(i)==0) then      
           v_cen_x(i)=0.0
           v_cen_y(i)=0.0
           do j=1,8
               k=ia(j,i)     !第i个单元第j个节点对应的整体节点号
               v_cen_x(i)=v_cen_x(i)+shape_n(j)*v(1,k)   !第i个单元任一点的速度,当l,h,m为(0,0,0)时,表示每单元中心的现时速度
	           v_cen_y(i)=v_cen_y(i)+shape_n(j)*v(2,k)
               v_cen_z(i)=v_cen_z(i)+shape_n(j)*v(3,k)

           end do
	   end if


   end do

   

!计算能量
   call energy_compute()

   return

end subroutine element_compute



subroutine elastic_try(i,dstrain_123_local,dstrain_d_123_local,      & !增量形式
                       try_p,try_stress_123_local,   &
					   stress_d_123_old_local, &  
					   e_old,p_old)
 
   implicit none



   integer(4) i,j,k
   real(8) dstrain_123_local(3,3)!在增量形式eos中
   real(8) dstrain_d_123_local(3,3)!在增量形式eos中
   real(8) try_p
   real(8) try_stress_123_local(3,3)
   real(8) stress_d_123_old_local(3,3)!在增量形式eos中
   real(8) e_old,p_old


   real(8) dstress_123_local(3,3)    !仅限于本函数
   real(8) dstress_d_123_local(3,3)
   real(8) try_sd_123(3,3)

!计算应力和偏应力
  
  dstress_123_local=0           

  do j=1,3
	 do k=1,3
			dstress_123_local(j,j)=dstress_123_local(j,j)+c123(j,k,i)*dstrain_123_local(k,k)
	 end do
  end do
	        dstress_123_local(1,2)=c123(4,4,i)*dstrain_123_local(1,2)	
			dstress_123_local(1,3)=c123(5,5,i)*dstrain_123_local(1,3)
			dstress_123_local(2,3)=c123(6,6,i)*dstrain_123_local(2,3)

			dstress_123_local(2,1)=dstress_123_local(1,2)
			dstress_123_local(3,1)=dstress_123_local(1,3)
			dstress_123_local(3,2)=dstress_123_local(2,3)
  


!算出的应力增量也没问题
		   diag_vector=0
		   diag_sum=0
		   diag_matrix=0

		   call diag_extract(dstress_123_local,diag_vector,3)
		   do j=1,3
				diag_sum=diag_sum+diag_vector(j)
		   end do
				diag_sum=diag_sum/3
				diag_vector=diag_sum
		   call	diag_enchase(diag_matrix,diag_vector,3)

          dstress_d_123_local=dstress_123_local-diag_matrix



 
!求出了偏应力在弹性段的增量，利用本构关系

!	stress_d_123(:,:,i)=stress_d_123(:,:,i)+dstress_d_123_local                    !i下123系中应力偏量增量，还没有确定是否塑性，这一步是否恰当？？？ans:只是做了下桥梁
!	stress_123(:,:,i)=stress_123(:,:,i)+dstress_123_local    
!------------------------------------------------
  if(problem_mode==0) then !求压力

    call dp_gruneisen(i,e_old,p_old,dstrain_d_123_local,&
                        stress_d_123_old_local,&
						 0.0,0.0,0.0)!碰撞问题只有压缩

  
  else if(problem_mode==1) then
       if(rho(i)>=rho0(ne_layer(i)))  then            
           call dp_gruneisen(i,e_old,p_old,dstrain_d_123_local, &
                             stress_d_123_old_local,&
						     0.0,0.0,0.0)!碰撞问题只有压缩
       else
           call dp_puff_solid(i,e_old,p_old,dstrain_d_123_local, &
                             stress_d_123_old_local,&
						      0.0,0.0,0.0)
       end if
   end if
   
 
  ! try_sd_123=stress_d_123(:,:,i)!try_sd_123为偏应力的中间量
  try_p=p(i)
  diag_matrix=0
  diag_vector=p(i)+q(i)
  call	diag_enchase(diag_matrix,diag_vector,3)

!------------------------------------------------

   try_stress_123_local=stress_d_123_old_local+dstress_d_123_local-diag_matrix     !绕了一大圈，完成了球，偏分解

 !-------------------------------------------------
   return
end subroutine elastic_try  !目标就是计算try_stress_123_local 

subroutine f1_mises_yield_unrelated(i,strs_123 ,     &
			                        f_mises,f1_mises)!需要看看三维的mises屈服准则是如何滴。
   implicit none
   integer(4) i
   real(8) strs_123(3,3)
   real(8) f_mises,f1_mises

   f_mises=((strs_123(1,1)-strs_123(2,2))**2+(strs_123(2,2)-strs_123(3,3))**2+(strs_123(3,3)-strs_123(1,1))**2       &
             +6.0*strs_123(1,2)**2+6.0*strs_123(1,3)**2+6.0*strs_123(2,3)**2)/strs_yld(i)+0.0000001       !strs_yld(i)拿0.5e-2代替了下，改回来记得
                 
   if(f_mises<0.0) stop "f_mises < 0.0, sqrt(f) error."
   f1_mises=sqrt(f_mises)-1
   f1_yield(i)=f1_mises

!   write(115,*) f1_yield(i),i,t
   return 
   end subroutine f1_mises_yield_unrelated



subroutine f1_hill_yield_unrelated(i,xx,strs_123,f_hill,f1_hill)
   implicit none
   integer(4) i
   real(8) xx
   real(8) y11,y22,y33
   real(8) strs_123(3,3)
   real(8) f_hill,f1_hill

   y11=1.0/(strs_yld_123(1,1,i)**2)-1.0/(strs_yld_123(2,2,i)**2)-1.0/(strs_yld_123(3,3,i)**2)     
   y22=1.0/(strs_yld_123(2,2,i)**2)-1.0/(strs_yld_123(3,3,i)**2)-1.0/(strs_yld_123(1,1,i)**2)
   y33=1.0/(strs_yld_123(3,3,i)**2)-1.0/(strs_yld_123(1,1,i)**2)-1.0/(strs_yld_123(2,2,i)**2)

   xx=(strs_yld_123(1,1,i)+strs_yld_123(2,2,i)+strs_yld_123(3,3,i))/3.0

   f_hill=(strs_123(1,1)**2)/(strs_yld_123(1,1,i)**2)+(strs_123(2,2)**2)/(strs_yld_123(2,2,i)**2)+(strs_123(3,3)**2)/(strs_yld_123(3,3,i)**2)        &
            +y11*strs_123(2,2)*strs_123(3,3)+y22*strs_123(1,1)*strs_123(3,3)+y33*strs_123(1,1)*strs_123(2,2)                &
		    +(strs_123(1,2)**2)/(strs_yld_123(1,2,i)**2)+(strs_123(1,3)**2)/(strs_yld_123(1,3,i)**2)+(strs_123(2,3)**2)/(strs_yld_123(2,3,i)**2)+0.00001
 !  write(*,*)f_hill,i,xyz(1,i),t
   !边界处有问题-------------------------
   if(f_hill<0.0) then
    f_hill=0.0001
	end if
	!--------------------------
   !stop "f_hill < 0.0, sqrt(f) error."
   f1_hill=xx*sqrt(f_hill/6.0)-xx/(6**0.5)
   f1_yield(i)=f1_hill



   return
end subroutine f1_hill_yield_unrelated

subroutine f1_mises_yield_related(i,strs_123,f_mises,f1_mises)
implicit none
integer(4) i  
!real(8) rate
real(8) strs_123(3,3)
real(8) f_mises,f1_mises

!    if(strain_p_rate_eff(i)<=strain_rate0(ne_layer(i)))  then
!      rate(i)=1.0
!    else
!      rate(i)=1+beta(ne_layer(i))*dlog(strain_p_rate_eff(i)/strain_rate0(ne_layer(i)))
!    end if
	
!	f_mises=((strs_11**2)+(strs_22**2)+(strs_33**2)-strs_22*strs_33-strs_11*strs_33       &
!             -strs_11*strs_22+3*strs_12**2 )/(strs_yld(i)**2)-1.0
    f_mises=((strs_123(1,1)-strs_123(2,2))**2+(strs_123(2,2)-strs_123(3,3))**2+(strs_123(3,3)-strs_123(1,1))**2       &
             +6.0*strs_123(1,2)**2+6.0*strs_123(1,3)**2+6.0*strs_123(2,3)**2)/(2.0*strs_yld(i)**2)+0.0000001 
	 
    f1_mises=f_mises/(rate(i)**2)-1.0                           !是否和率无关一样需要sqrt，这个很关键的问题，或者说sqrt是否有必要，更进一步是否在plastic_mises中有2反制措施
	f1_yield(i)=f1_mises
    
    return
end subroutine f1_mises_yield_related


subroutine f1_hill_yield_related(i,xx,strs_123,f_hill,f1_hill)
   implicit none
   integer(4) i
   real(8) xx
   real(8) y11,y22,y33
   real(8) strs_123(3,3)
   real(8) f_hill,f1_hill

   y11=1.0/(strs_yld_123(1,1,i)**2)-1.0/(strs_yld_123(2,2,i)**2)-1.0/(strs_yld_123(3,3,i)**2)     
   y22=1.0/(strs_yld_123(2,2,i)**2)-1.0/(strs_yld_123(3,3,i)**2)-1.0/(strs_yld_123(1,1,i)**2)
   y33=1.0/(strs_yld_123(3,3,i)**2)-1.0/(strs_yld_123(1,1,i)**2)-1.0/(strs_yld_123(2,2,i)**2)

   xx=(strs_yld_123(1,1,i)+strs_yld_123(2,2,i)+strs_yld_123(3,3,i))/3.0

   f_hill=(strs_123(1,1)**2)/(strs_yld_123(1,1,i)**2)+(strs_123(2,2)**2)/(strs_yld_123(2,2,i)**2)+(strs_123(3,3)**2)/(strs_yld_123(3,3,i)**2)        &
            +y11*strs_123(2,2)*strs_123(3,3)+y22*strs_123(1,1)*strs_123(3,3)+y33*strs_123(1,1)*strs_123(2,2)                &
		    +(strs_123(1,2)**2)/(strs_yld_123(1,2,i)**2)+(strs_123(1,3)**2)/(strs_yld_123(1,3,i)**2)+(strs_123(2,3)**2)/(strs_yld_123(2,3,i)**2)+0.00001
 !  write(*,*)f_hill,i,xyz(1,i),t
   !边界处有问题-------------------------
   if(f_hill<0.0) then
    f_hill=0.0001
	end if
	!--------------------------
   !stop "f_hill < 0.0, sqrt(f) error."
   f1_hill=xx*sqrt(f_hill/6.0)-xx/(6**0.5)
   f1_yield(i)=f1_hill
   return
end subroutine f1_hill_yield_related  

 
subroutine plastic_compute_unrelated_mises(i,strs_123,      & !弹性试计算应力的带入就是为了在子程序中再次带入求解mises值
				                             dstrain_123_local,dstrain_d_123_local,        & !只用了应变增量和上一步计算的应力偏量
									         stress_d_123_old_local,    &
											 e_old,p_old)

implicit none
integer(4) i,j
real(8) dlamda(1,1)                 !   是个数值！
real(8) F_sigma(1,6)
real(8) d123(3,3)
real(8) f_mises,f1_mises


real(8) strs_123(3,3)                 !试计算应力
real(8) dstrain_123_local(3,3) !这里有问题
real(8) dstrain_d_123_local(3,3)
real(8) stress_d_123_old_local(3,3)
real(8) e_old,p_old

real(8) dstress_d_123(3,3)
real(8) dstress_d_column(6,1)
 

real dstrain_p_123(3,3)
real(8) dstrain_e_123(3,3)     

real(8) dstrain_123_column(6,1)
real(8) dstrain_e_column(6,1)



    call f1_mises_yield_unrelated(i,strs_123,f_mises,f1_mises)



d123(1,1)=(2*strs_123(1,1)-strs_123(2,2)-strs_123(3,3))/(strs_yld(i)**2)/(2*sqrt(f_mises))!计算导数
d123(2,2)=(2*strs_123(2,2)-strs_123(1,1)-strs_123(3,3))/(strs_yld(i)**2)/(2*sqrt(f_mises))
d123(3,3)=(2*strs_123(3,3)-strs_123(2,2)-strs_123(1,1))/(strs_yld(i)**2)/(2*sqrt(f_mises))
d123(1,2)=6*strs_123(1,2)/(strs_yld(i)**2)/(2*sqrt(f_mises))
d123(1,3)=6*strs_123(1,3)/(strs_yld(i)**2)/(2*sqrt(f_mises))
d123(2,3)=6*strs_123(2,3)/(strs_yld(i)**2)/(2*sqrt(f_mises))
d123(2,1)=d123(1,2)
d123(3,1)=d123(1,3)
d123(3,2)=d123(2,3)


F_sigma(1,1)=d123(1,1)
F_sigma(1,2)=d123(2,2)
F_sigma(1,3)=d123(3,3)
F_sigma(1,4)=d123(1,2)
F_sigma(1,5)=d123(1,3)
F_sigma(1,6)=d123(2,3)

dstrain_123_column(1,1)=dstrain_123_local(1,1)
dstrain_123_column(2,1)=dstrain_123_local(2,2)
dstrain_123_column(3,1)=dstrain_123_local(3,3)
dstrain_123_column(4,1)=dstrain_123_local(1,2)
dstrain_123_column(5,1)=dstrain_123_local(1,3)
dstrain_123_column(6,1)=dstrain_123_local(2,3)



 !f1_mises关于6个应力分量的导数F_sigma   

    dlamda=matmul(matmul(F_sigma,c123(:,:,i)),dstrain_123_column)/matmul(matmul(F_sigma,c123(:,:,i)),transpose(F_sigma))
   

    dstrain_p_123=dlamda(1,1)*d123                                         !应变增量的塑性部分利用增量理论

    strain_p_123(:,:,i)=strain_p_123(:,:,i)+dstrain_p_123
	dstrain_e_123=dstrain_123_local-dstrain_p_123                           !应变增量的弹性部分利用总的应变增量减去塑性部分                    
    
	dstrain_e_column(1,1)=dstrain_e_123(1,1)
	dstrain_e_column(2,1)=dstrain_e_123(2,2)
	dstrain_e_column(3,1)=dstrain_e_123(3,3)
	dstrain_e_column(4,1)=dstrain_e_123(1,2)
	dstrain_e_column(5,1)=dstrain_e_123(1,3)
	dstrain_e_column(6,1)=dstrain_e_123(2,3)


    dstress_d_column=matmul(c123(:,:,i),dstrain_e_column)
	
    dstress_d_123=0
	dstress_d_123(1,2)=dstress_d_column(4,1)
	dstress_d_123(1,3)=dstress_d_column(5,1)
	dstress_d_123(2,3)=dstress_d_column(6,1)
    dstress_d_123=dstress_d_123+transpose(dstress_d_123)
    dstress_d_123(1,1)=dstress_d_column(1,1)
	dstress_d_123(2,2)=dstress_d_column(2,1)
	dstress_d_123(3,3)=dstress_d_column(3,1)

		  diag_vector=0
		  diag_sum=0
		  diag_matrix=0
		  call diag_extract(dstress_d_123,diag_vector,3)
		   do j=1,3
				diag_sum=diag_sum+diag_vector(j)
		   end do
				diag_sum=diag_sum/3
				diag_vector=diag_sum
		   call	diag_enchase(diag_matrix,diag_vector,3)

    dstress_d_123=dstress_d_123-diag_matrix
	stress_d_123(:,:,i)=stress_d_123_old_local+dstress_d_123										!其实可以跳过stress_d_123_old_local这个超级麻烦的

    if(problem_mode==0) then
       call dp_gruneisen(i,e_old,p_old,dstrain_d_123_local,      &
                         stress_d_123_old_local,                 &
                         dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
	else if(problem_mode==1) then
	    if(rho(i)>=rho0(ne_layer(i))) then
            call dp_gruneisen(i,e_old,p_old,dstrain_d_123_local, &
                              stress_d_123_old_local,&
                              dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
        else
            call dp_puff_solid(i,e_old,p_old,dstrain_d_123_local, &
                              stress_d_123_old_local,&
                              dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
        end if
	end if

   
    diag_vector=0
	diag_vector=p(i)+q(i)
	call diag_enchase(diag_matrix,diag_vector,3)
	stress_123(:,:,i)=stress_d_123(:,:,i)-diag_matrix                !完成了应力的计算

end subroutine plastic_compute_unrelated_mises 
    
 


subroutine plastic_compute_unrelated_hill(i,strs_123,      & !弹性试计算应力的带入就是为了在子程序中再次带入求解mises值
				                            dstrain_123_local,dstrain_d_123_local,        & !只用了应变增量和上一步计算的应力偏量
									        stress_d_123_old_local,    &
											e_old,p_old)
implicit none 
integer(4) i,j
real  dlamda(1,1)
real(8) F_sigma(1,6)
real(8) xx                            !xx是个什么玩意啊！！！
real(8) y11,y22,y33
real(8) d123(3,3)
real(8) f_hill,f1_hill

real(8) strs_123(3,3)                 !试计算应力
real(8) dstrain_123_local(3,3)        
real(8) dstrain_d_123_local(3,3)
real(8) stress_d_123_old_local(3,3)
real(8) e_old,p_old

real(8) dstress_d_123(3,3)
real(8) dstress_d_column(6,1)
 

real dstrain_p_123(3,3)
real(8) dstrain_e_123(3,3)     

real(8) dstrain_123_column(6,1)
real(8) dstrain_e_column(6,1)


   y11=1.0/(strs_yld_123(1,1,i)**2)-1.0/(strs_yld_123(2,2,i)**2)-1.0/(strs_yld_123(3,3,i)**2)     
   y22=1.0/(strs_yld_123(2,2,i)**2)-1.0/(strs_yld_123(3,3,i)**2)-1.0/(strs_yld_123(1,1,i)**2)
   y33=1.0/(strs_yld_123(3,3,i)**2)-1.0/(strs_yld_123(1,1,i)**2)-1.0/(strs_yld_123(2,2,i)**2)
  
   call f1_hill_yield_unrelated(i,xx,strs_123,f_hill,f1_hill)

	d123(1,1)=(xx/(6**0.5))*(2*strs_123(1,1)/(strs_yld_123(1,1,i)**2)+y33*strs_123(2,2)+y22*strs_123(3,3))/(2*sqrt(f_hill))!计算导数
	d123(2,2)=(xx/(6**0.5))*(2*strs_123(2,2)/(strs_yld_123(2,2,i)**2)+y11*strs_123(3,3)+y33*strs_123(1,1))/(2*sqrt(f_hill))
	d123(3,3)=(xx/(6**0.5))*(2*strs_123(3,3)/(strs_yld_123(3,3,i)**2)+y11*strs_123(2,2)+y22*strs_123(1,1))/(2*sqrt(f_hill))
	d123(1,2)=(xx/(6**0.5))*(2*strs_123(1,2)/(strs_yld_123(1,2,i)**2))/(2*sqrt(f_hill))
	d123(1,3)=(xx/(6**0.5))*(2*strs_123(1,3)/(strs_yld_123(1,3,i)**2))/(2*sqrt(f_hill))
	d123(2,3)=(xx/(6**0.5))*(2*strs_123(2,3)/(strs_yld_123(2,3,i)**2))/(2*sqrt(f_hill))
	d123(2,1)=d123(1,2)
	d123(3,1)=d123(1,3)
	d123(3,2)=d123(2,3)

	F_sigma(1,1)=d123(1,1)
	F_sigma(1,2)=d123(2,2)
	F_sigma(1,3)=d123(3,3)
	F_sigma(1,4)=d123(1,2)
	F_sigma(1,5)=d123(1,3)
	F_sigma(1,6)=d123(2,3)

	dstrain_123_column(1,1)=dstrain_123_local(1,1)
	dstrain_123_column(2,1)=dstrain_123_local(2,2)
	dstrain_123_column(3,1)=dstrain_123_local(3,3)
	dstrain_123_column(4,1)=dstrain_123_local(1,2)
	dstrain_123_column(5,1)=dstrain_123_local(1,3)
	dstrain_123_column(6,1)=dstrain_123_local(2,3)

  !d11,d22,d33,d12为f1对应力的偏导数
    dlamda=matmul(matmul(F_sigma,c123(:,:,i)),dstrain_123_column)/matmul(matmul(F_sigma,c123(:,:,i)),transpose(F_sigma))

   dstrain_p_123=dlamda(1,1)*d123

     
  strain_p_123(:,:,i)=strain_p_123(:,:,i)+dstrain_p_123
	dstrain_e_123=dstrain_123_local-dstrain_p_123
 
                           !应变增量的弹性部分利用总的应变增量减去塑性部分                    
    
	dstrain_e_column(1,1)=dstrain_e_123(1,1)
	dstrain_e_column(2,1)=dstrain_e_123(2,2)
	dstrain_e_column(3,1)=dstrain_e_123(3,3)
	dstrain_e_column(4,1)=dstrain_e_123(1,2)
	dstrain_e_column(5,1)=dstrain_e_123(1,3)
	dstrain_e_column(6,1)=dstrain_e_123(2,3)
	
    dstress_d_column=matmul(c123(:,:,i),dstrain_e_column)
	
    dstress_d_123=0
	dstress_d_123(1,2)=dstress_d_column(4,1)
	dstress_d_123(1,3)=dstress_d_column(5,1)
	dstress_d_123(2,3)=dstress_d_column(6,1)
    dstress_d_123=dstress_d_123+transpose(dstress_d_123)
    dstress_d_123(1,1)=dstress_d_column(1,1)
	dstress_d_123(2,2)=dstress_d_column(2,1)
	dstress_d_123(3,3)=dstress_d_column(3,1)


		  diag_vector=0
		  diag_sum=0
		  diag_matrix=0
		   call diag_extract(dstress_d_123,diag_vector,3)
		   do j=1,3
				diag_sum=diag_sum+diag_vector(j)
		   end do
				diag_sum=diag_sum/3
				diag_vector=diag_sum
		   call	diag_enchase(diag_matrix,diag_vector,3)

    dstress_d_123=dstress_d_123-diag_matrix
	stress_d_123(:,:,i)=stress_d_123_old_local+dstress_d_123										!其实可以跳过stress_d_123_old_local这个超级麻烦的

    if(problem_mode==0) then
       call dp_gruneisen(i,e_old,p_old,dstrain_d_123_local,      &
                         stress_d_123_old_local,                 &
                         dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
	else if(problem_mode==1) then
	    if(rho(i)>=rho0(ne_layer(i))) then
            call dp_gruneisen(i,e_old,p_old,dstrain_d_123_local, &
                              stress_d_123_old_local,&
                              dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
        else
            call dp_puff_solid(i,e_old,p_old,dstrain_d_123_local, &
                              stress_d_123_old_local,&
                              dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
        end if
	end if

    diag_vector=0
	diag_vector=p(i)+q(i)
	call diag_enchase(diag_matrix,diag_vector,3)
	stress_123(:,:,i)=stress_d_123(:,:,i)-diag_matrix                !完成了应力的计算
end subroutine plastic_compute_unrelated_hill

subroutine plastic_compute_related_mises(i,strs_123,      & !弹性试计算应力的带入就是为了在子程序中再次带入求解mises值
				                           dstrain_123_local,dstrain_d_123_local,        & !只用了应变增量和上一步计算的应力偏量
									       stress_d_123_old_local,    &
										   e_old,p_old)

implicit none
integer i,j,k
real(8) dlamda(1,1),elta
real(8) F_sigma(1,6)
real(8) d123(3,3)
real(8) f_mises,f1_mises

real(8) strs_123(3,3)                 !试计算应力
real(8) dstrain_123_local(3,3) !这里有问题
real(8) dstrain_d_123_local(3,3)
real(8) stress_d_123_old_local(3,3)
real(8) e_old,p_old

real(8) dstress_d_123(3,3)
real(8) dstress_d_column(6,1)
 
real dstrain_p_123(3,3)
real(8) dstrain_e_123(3,3)     
real(8) dstrain_e_column(6,1)

real(8)::dstrain_p_eff

real(8) dstrain_123_column(6,1)





    call f1_mises_yield_related(i,strs_123,f_mises,f1_mises)   !再次计算f是为了代入上一时刻的应力值而不是试算应力值
    if(f1_mises< 0.0) stop "f< 0.0, error."

	d123(1,1)=(2*strs_123(1,1)-strs_123(2,2)-strs_123(3,3))/(strs_yld(i)**2)!计算导数
	d123(2,2)=(2*strs_123(2,2)-strs_123(1,1)-strs_123(3,3))/(strs_yld(i)**2)
	d123(3,3)=(2*strs_123(3,3)-strs_123(2,2)-strs_123(1,1))/(strs_yld(i)**2)
	d123(1,2)=6*strs_123(1,2)/(strs_yld(i)**2)
	d123(1,3)=6*strs_123(1,3)/(strs_yld(i)**2)
	d123(2,3)=6*strs_123(2,3)/(strs_yld(i)**2)
	d123(2,1)=d123(1,2)
	d123(3,1)=d123(1,3)
	d123(3,2)=d123(2,3)


    elta=1.0/(2*beta0(ne_layer(i))*sqrt(f_mises))
	dlamda=0   
	do k=1,3
		do j=1,3
		dlamda(1,1)=dlamda(1,1)+d123(i,j)*d123(i,j)
		end do
	end do                                          
    dlamda(1,1)=dt/(elta*sqrt(2.0/3.0*dlamda(1,1))) 

	F_sigma(1,1)=d123(1,1)
	F_sigma(1,2)=d123(2,2)
	F_sigma(1,3)=d123(3,3)
	F_sigma(1,4)=d123(1,2)
	F_sigma(1,5)=d123(1,3)
	F_sigma(1,6)=d123(2,3)
	               
	dstrain_p_123=dlamda(1,1)*strain_rate0(ne_layer(i))*exp((sqrt(f_mises)-1.0)/beta0(ne_layer(i)))*elta*d123

	strain_p_123(:,:,i)=strain_p_123(:,:,i)+dstrain_p_123
	dstrain_e_123=dstrain_123_local-dstrain_p_123   

	dstrain_e_column(1,1)=dstrain_e_123(1,1)
	dstrain_e_column(2,1)=dstrain_e_123(2,2)
	dstrain_e_column(3,1)=dstrain_e_123(3,3)
	dstrain_e_column(4,1)=dstrain_e_123(1,2)
	dstrain_e_column(5,1)=dstrain_e_123(1,3)
	dstrain_e_column(6,1)=dstrain_e_123(2,3)


	dstrain_p_eff=0
	do k=1,3
		do j=1,3
		dstrain_p_eff=dstrain_p_eff+dstrain_p_123(i,j)*dstrain_p_123(i,j)
		end do
	end do  
    dstrain_p_eff=sqrt(2.0/3.0*dstrain_p_eff)        !计算等效塑性应变增量
	                   
    strain_p_eff(i)=strain_p_eff(i)+dstrain_p_eff

    strain_p_rate_eff(i)=dstrain_p_eff/dt

 dstress_d_column=matmul(c123(:,:,i),dstrain_e_column)
	
    dstress_d_123=0
	dstress_d_123(1,2)=dstress_d_column(4,1)
	dstress_d_123(1,3)=dstress_d_column(5,1)
	dstress_d_123(2,3)=dstress_d_column(6,1)
    dstress_d_123=dstress_d_123+transpose(dstress_d_123)

    dstress_d_123(1,1)=dstress_d_column(1,1)
	dstress_d_123(2,2)=dstress_d_column(2,1)
	dstress_d_123(3,3)=dstress_d_column(3,1)

		  diag_vector=0
		  diag_sum=0
		  diag_matrix=0
		   call diag_extract(dstress_d_123,diag_vector,3)
		   do j=1,3
				diag_sum=diag_sum+diag_vector(j)
		   end do
				diag_sum=diag_sum/3
				diag_vector=diag_sum
			call diag_enchase(diag_matrix,diag_vector,3)	

    dstress_d_123=dstress_d_123-diag_matrix
	stress_d_123(:,:,i)=stress_d_123_old_local+dstress_d_123										!其实可以跳过stress_d_123_old_local这个超级麻烦的

    if(problem_mode==0) then
       call dp_gruneisen(i,e_old,p_old,dstrain_d_123_local,      &
                         stress_d_123_old_local,                 &
                         dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
	else if(problem_mode==1) then
	    if(rho(i)>=rho0(ne_layer(i))) then
            call dp_gruneisen(i,e_old,p_old,dstrain_d_123_local, &
                              stress_d_123_old_local,&
                              dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
        else
            call dp_puff_solid(i,e_old,p_old,dstrain_d_123_local, &
                              stress_d_123_old_local,&
                              dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
        end if
	end if

    diag_vector=0
	diag_vector=p(i)+q(i)
	call diag_enchase(diag_matrix,diag_vector,3)
	stress_123(:,:,i)=stress_d_123(:,:,i)-diag_matrix                !完成了应力的计算
end subroutine plastic_compute_related_mises   


!计算塑性应变增量,应力 
subroutine plastic_compute_related_hill(i,strs_123,      & !弹性试计算应力的带入就是为了在子程序中再次带入求解mises值
				                          dstrain_123_local,dstrain_d_123_local,        & !只用了应变增量和上一步计算的应力偏量
									      stress_d_123_old_local,    &
										  e_old,p_old)
implicit none 
integer(4) i,j
real  dlamda(1,1)
real(8) F_sigma(1,6)
real(8) xx                            !xx是个什么玩意啊！！！
real(8) y11,y22,y33
real(8) d123(3,3)
real(8) f_hill,f1_hill

real(8) strs_123(3,3)                 !试计算应力
real(8) dstrain_123_local(3,3)        
real(8) dstrain_d_123_local(3,3)
real(8) stress_d_123_old_local(3,3)
real(8) e_old,p_old

real(8) dstress_d_123(3,3)
real(8) dstress_d_column(6,1)
 

real dstrain_p_123(3,3)
real(8) dstrain_e_123(3,3)     

real(8) dstrain_123_column(6,1)
real(8) dstrain_e_column(6,1)


   y11=1.0/(strs_yld_123(1,1,i)**2)-1.0/(strs_yld_123(2,2,i)**2)-1.0/(strs_yld_123(3,3,i)**2)     
   y22=1.0/(strs_yld_123(2,2,i)**2)-1.0/(strs_yld_123(3,3,i)**2)-1.0/(strs_yld_123(1,1,i)**2)
   y33=1.0/(strs_yld_123(3,3,i)**2)-1.0/(strs_yld_123(1,1,i)**2)-1.0/(strs_yld_123(2,2,i)**2)
  
   call f1_hill_yield_unrelated(i,xx,strs_123,f_hill,f1_hill)

	d123(1,1)=(xx/(6**0.5))*(2*strs_123(1,1)/(strs_yld_123(1,1,i)**2)+y33*strs_123(2,2)+y22*strs_123(3,3))/(2*sqrt(f_hill))!计算导数
	d123(2,2)=(xx/(6**0.5))*(2*strs_123(2,2)/(strs_yld_123(2,2,i)**2)+y11*strs_123(3,3)+y33*strs_123(1,1))/(2*sqrt(f_hill))
	d123(3,3)=(xx/(6**0.5))*(2*strs_123(3,3)/(strs_yld_123(3,3,i)**2)+y11*strs_123(2,2)+y22*strs_123(1,1))/(2*sqrt(f_hill))
	d123(1,2)=(xx/(6**0.5))*(2*strs_123(1,2)/(strs_yld_123(1,2,i)**2))/(2*sqrt(f_hill))
	d123(1,3)=(xx/(6**0.5))*(2*strs_123(1,3)/(strs_yld_123(1,3,i)**2))/(2*sqrt(f_hill))
	d123(2,3)=(xx/(6**0.5))*(2*strs_123(2,3)/(strs_yld_123(2,3,i)**2))/(2*sqrt(f_hill))
	d123(2,1)=d123(1,2)
	d123(3,1)=d123(1,3)
	d123(3,2)=d123(2,3)

	F_sigma(1,1)=d123(1,1)
	F_sigma(1,2)=d123(2,2)
	F_sigma(1,3)=d123(3,3)
	F_sigma(1,4)=d123(1,2)
	F_sigma(1,5)=d123(1,3)
	F_sigma(1,6)=d123(2,3)

	dstrain_123_column(1,1)=dstrain_123_local(1,1)
	dstrain_123_column(2,1)=dstrain_123_local(2,2)
	dstrain_123_column(3,1)=dstrain_123_local(3,3)
	dstrain_123_column(4,1)=dstrain_123_local(1,2)
	dstrain_123_column(5,1)=dstrain_123_local(1,3)
	dstrain_123_column(6,1)=dstrain_123_local(2,3)

  !d11,d22,d33,d12为f1对应力的偏导数
    dlamda=matmul(matmul(F_sigma,c123(:,:,i)),dstrain_123_column)/matmul(matmul(F_sigma,c123(:,:,i)),transpose(F_sigma))

   dstrain_p_123=dlamda(1,1)*d123

  strain_p_123(:,:,i)=strain_p_123(:,:,i)+dstrain_p_123
	dstrain_e_123=dstrain_123_local-dstrain_p_123                           !应变增量的弹性部分利用总的应变增量减去塑性部分                    
    
	dstrain_e_column(1,1)=dstrain_e_123(1,1)
	dstrain_e_column(2,1)=dstrain_e_123(2,2)
	dstrain_e_column(3,1)=dstrain_e_123(3,3)
	dstrain_e_column(4,1)=dstrain_e_123(1,2)
	dstrain_e_column(5,1)=dstrain_e_123(1,3)
	dstrain_e_column(6,1)=dstrain_e_123(2,3)
	
    dstress_d_column=matmul(c123(:,:,i),dstrain_e_column)
	
    dstress_d_123=0
	dstress_d_123(1,2)=dstress_d_column(4,1)
	dstress_d_123(1,3)=dstress_d_column(5,1)
	dstress_d_123(2,3)=dstress_d_column(6,1)
    dstress_d_123=dstress_d_123+transpose(dstress_d_123)
    dstress_d_123(1,1)=dstress_d_column(1,1)
	dstress_d_123(2,2)=dstress_d_column(2,1)
	dstress_d_123(3,3)=dstress_d_column(3,1)


		  diag_vector=0
		  diag_sum=0
		  diag_matrix=0
		   call diag_extract(dstress_d_123,diag_vector,3)
		   do j=1,3
				diag_sum=diag_sum+diag_vector(j)
		   end do
				diag_sum=diag_sum/3
				diag_vector=diag_sum
		   call	diag_enchase(diag_matrix,diag_vector,3)

    dstress_d_123=dstress_d_123-diag_matrix
	stress_d_123(:,:,i)=stress_d_123_old_local+dstress_d_123										!其实可以跳过stress_d_123_old_local这个超级麻烦的

    if(problem_mode==0) then
       call dp_gruneisen(i,e_old,p_old,dstrain_d_123_local,      &
                         stress_d_123_old_local,                 &
                         dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
	else if(problem_mode==1) then
	    if(rho(i)>=rho0(ne_layer(i))) then
            call dp_gruneisen(i,e_old,p_old,dstrain_d_123_local, &
                              stress_d_123_old_local,&
                              dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
        else
            call dp_puff_solid(i,e_old,p_old,dstrain_d_123_local, &
                              stress_d_123_old_local,&
                              dstrain_p_123(1,1),dstrain_p_123(2,2),dstrain_p_123(3,3))
        end if
	end if

    diag_vector=0
	diag_vector=p(i)+q(i)
	call diag_enchase(diag_matrix,diag_vector,3)
	stress_123(:,:,i)=stress_d_123(:,:,i)-diag_matrix     

end subroutine plastic_compute_related_hill    

  
!弹塑性条件下压力与内能的求解   
subroutine dp_gruneisen(i,e_old,p_old,dstrain_d_123,    & !核心还是能量守恒,实际上都没用到dstrain_123_local啊
                         stress_d_123_old,                                      &
						 dstrain_p_11,dstrain_p_22,dstrain_p_33)
implicit none
integer i,j,k
real(8) a1,a2,a3,a1_1,ea
real(8) b1,b2

real(8) sd_av_123(3,3)


real(8) strain_vol,dstrain_vol
real(8) e_old,p_old
real(8) dstrain_d_123(3,3)
real(8) stress_d_123_old(3,3)
real dstrain_p_11,dstrain_p_22,dstrain_p_33
real(8) dst,dsp 


   strain_vol=(vol(i)-vol0(i))/vol(i)                    !体应变, 注意这里的theta同文章中的符号取反的
   dstrain_vol=dvol(i)/vol(i)                            !体应变的改变量，增量形式中需要
  

         
   a1=rho0(ne_layer(i))*c0(ne_layer(i))**2.0
   a2=a1*(2.0*s(ne_layer(i))-1.0)
   a3=a1*(3.0*s(ne_layer(i))**2.0-4.0*s(ne_layer(i))+1.0)
   
   if(mat_model(ne_layer(i))==0) then
       a1_1=a1
   else if (mat_model(ne_layer(i))==1) then
       a1_1=(c11(i)+c22(i)+c33(i)+2.0*c12(i)+2.0*c13(i)+2.0*c23(i))/9.0

   end if

   sd_av_123=0.5*(stress_d_123(:,:,i)+stress_d_123_old)
  
   !对于形式的p=b1+b2*e的eos方程，引入中间变量ea,可得到p=(b1+b2*ea)/(1+0.5*b2*dv)这个需要好好研究下
if(eos_mode==1)then
    dst=-(c11(i)+c12(i)+c13(i))*dstrain_d_123(1,1)/3.0-(c12(i)+c22(i)+c23(i))*dstrain_d_123(2,2)/3.0-(c13(i)+c23(i)+c33(i))*dstrain_d_123(3,3)/3.0
	dsp=+(c11(i)+c12(i)+c13(i))*dstrain_p_11/3.0+(c12(i)+c22(i)+c23(i))*dstrain_p_22/3.0+(c13(i)+c23(i)+c33(i))*dstrain_p_33/3.0
else if(eos_mode==0)then
	dst=0.0
	dsp=0.0
	a1_1=a1
end if




   b1=p_old+(-a1_1+2*(a2-gama(ne_layer(i))/2.0*a1_1)*strain_vol                   &
        -3*(a3-gama(ne_layer(i))/2.0*a2)*(strain_vol**2))*dstrain_vol                    & 
	    -(rho0(ne_layer(i))*gama(ne_layer(i))-rho0(ne_layer(i))*gama(ne_layer(i))*strain_vol)*e_old+dst+dsp
   b2=rho0(ne_layer(i))*gama(ne_layer(i))-rho0(ne_layer(i))*gama(ne_layer(i))*strain_vol-rho0(ne_layer(i))*gama(ne_layer(i))*dstrain_vol




        ea=0
		do j=1,3
			do k=1,3
				ea=ea+sd_av_123(j,k)*dstrain_d_123(j,k)!算形变能增量                             !
             end do
		end do
		ea=ea*vol_average(i)  
		ea=ea+e_old-(0.5*p_old+q_av(i))*dvol(i)        !请见系列张雄   没有任何问题
!    ea=e_old+(sd_av_11*dstrain_d_11                 &
!        +sd_av_22*dstrain_d_22                            &
!		+2*sd_av_12*dstrain_d_12)*vol_average(i)         &
!		-(0.5*p_old+q_av(i))*dvol(i)

!   dp(i)=(-a1_1+2*(a2-gama(ne_layer(i))/2.0*a1)*strain_vol                   &
!        -3*(a3-gama(ne_layer(i))/2.0*a2)*(strain_vol**2))*dstrain_vol                    &
!        -(c11(i)+c13(i)-c22(i)-c23(i))*dstrain_d_11/3.0                                  &
!	    +(rho0(ne_layer(i))*gama(ne_layer(i))-rho0(ne_layer(i))*gama(ne_layer(i))*strain_vol)*de(i)                           &
!	    -rho0(ne_layer(i))*gama(ne_layer(i))*e_old*dstrain_vol       &
!		+(c11(i)+c12(i)+c13(i))*dstrain_p_11/3.0            &   !塑性变形时平均压力的修正比弹性时多塑性应变相关部分
!        +(c12(i)+c22(i)+c23(i))*dstrain_p_22/3.0                  &
!        +(c13(i)+c23(i)+c33(i))*dstrain_p_33/3.0

   if(problem_mode==1) then
       ea=ea+energy_xray(i)     
   end if
	
             
   p(i)=(b1+b2*ea)/(1+0.5*b2*dvol(i))
   e(i)=(p(i)-b1)/b2       !单位质量上的内能
   dp(i)=p(i)-p_old
   de(i)=e(i)-e_old 	
!   p(i)=p_old+dp(i)
!   e(i)=e_old+de(i)

   if(i==ne_his_x(6))then          
    p_dst=p_dst+dst
	p_dsP=p_dsp+dsP

   open(888,file="test_dp_fun.txt")
   write(888,88)t,p(i),e(i),q_av(i),b1,b2,1,1,strain_vol,a1,a1_1,p_dst,p_dsp 
   88 format (6(f8.4,3x),i1,3x,i1,3x,5(f8.4,3x))
   
	end if


end subroutine dp_gruneisen 



!dp_puff_solid 为既考虑体积变化非线性，又考虑各向异性的修正PUFF物态方程
subroutine dp_puff_solid (i,e_old,p_old,dstrain_d_123,  & !核心还是能量守恒,实际上都没用到dstrain_123_local啊
                         stress_d_123_old,&
					     dstrain_p_11,dstrain_p_22,dstrain_p_33)
   implicit none
   integer(4) i,j,k

   real(8) sd_av_123(3,3)!求偏应力平均值

   
   real(8) x1,x2,x3
   real(8) x1_1
   real(8) c1,c2
   real(8) ea
   real(8) strain_vol,dstrain_vol
   real(8) e_old,p_old
   real(8) dstrain_d_123(3,3)
   real(8) stress_d_123_old(3,3)
   real dstrain_p_11,dstrain_p_22,dstrain_p_33
   real(8) dst,dsp 

!   strain_vol=strain_11(i)+strain_22(i)
!   dstrain_vol=dstrain_11+dstrain_22
   strain_vol=(vol(i)-vol0(i))/vol(i)
   dstrain_vol=dvol(i)/vol(i)

   x1=rho0(ne_layer(i))*c0(ne_layer(i))*c0(ne_layer(i))
   x2=-0.5*x1-0.5*h0(ne_layer(i))*x1/gama(ne_layer(i))+0.5*x1*cn(ne_layer(i))
   x3=5*x1/24.0+5.0*h0(ne_layer(i))*x1/(8.0*gama(ne_layer(i)))-1.25*x1*cn(ne_layer(i))          &
      -h0(ne_layer(i))*x1*cn(ne_layer(i))/(4.0*gama(ne_layer(i)))+x1*cn(ne_layer(i))*cn(ne_layer(i))/6.0
   
   x1_1=(c11(i)+c22(i)+c33(i)+2.0*c12(i)+2.0*c13(i)+2.0*c23(i))/9.0 
   
    sd_av_123=0.5*(stress_d_123(:,:,i)+stress_d_123_old)
 if(mat_model(ne_layer(i))==0) then
       x1_1=x1
   else if (mat_model(ne_layer(i))==1) then
       x1_1=x1_1
   end if          

if(eos_mode==1)then
    dst=-(c11(i)+c12(i)+c13(i))*dstrain_d_123(1,1)/3.0-(c12(i)+c22(i)+c23(i))*dstrain_d_123(2,2)/3.0-(c13(i)+c23(i)+c33(i))*dstrain_d_123(3,3)/3.0
	dsp=+(c11(i)+c12(i)+c13(i))*dstrain_p_11/3.0+(c12(i)+c22(i)+c23(i))*dstrain_p_22/3.0+(c13(i)+c23(i)+c33(i))*dstrain_p_33/3.0
else if(eos_mode==0)then
	dst=0.0
	dsp=0.0
	x1_1=x1
end if



   c1=p_old+(-x1_1+2.0*x2*strain_vol-3.0*x3*strain_vol*strain_vol)*dstrain_vol  &
     -(rho0(ne_layer(i))*gama(ne_layer(i))-1.5*rho0(ne_layer(i))*gama(ne_layer(i))*strain_vol    &
	 +0.5*rho0(ne_layer(i))*h0(ne_layer(i))*strain_vol)*e_old+dst+dsp        

  		 
   c2=rho0(ne_layer(i))*gama(ne_layer(i))-1.5*rho0(ne_layer(i))*gama(ne_layer(i))*strain_vol+0.5*rho0(ne_layer(i))*h0(ne_layer(i))*strain_vol      &
     -1.5*rho0(ne_layer(i))*gama(ne_layer(i))*dstrain_vol+0.5*rho0(ne_layer(i))*h0(ne_layer(i))*dstrain_vol

    
		
		ea=0
		do j=1,3
			do k=1,3
				ea=ea+sd_av_123(j,k)*dstrain_d_123(j,k)                             !
             end do
		end do   
		ea=0.5*ea+e_old-(0.5*p_old+q_av(i))*dvol(i)
			
!    ea=e_old+(sd_av_11*dstrain_d_11                 &
!        +sd_av_22*dstrain_d_22                            &
!		+2*sd_av_12*dstrain_d_12)*vol_average(i)         &
!		-(0.5*p_old+q_av(i))*dvol(i)


      
   if(problem_mode==1) then
       ea=ea+energy_xray(i)    
   end if


   p(i)=(c1+c2*ea)/(1+0.5*dvol(i)*c2)
   e(i)=(p(i)-c1)/c2
   
   dp(i)=p(i)-p_old
   de(i)=e(i)-e_old
    if(i==ne_his_x(6))then          
    p_dst=p_dst+dst
	p_dsP=p_dst+dsP
   open(888,file="test_dp_fun.txt")
   write(888,88)t,p(i),e(i),q_av(i),c1,c2,1,1,strain_vol,x1,x1_1,p_dst,p_dsp    
   88 format (6(f8.4,3x),i1,3x,i1,3x,5(f8.4,3x))

	end if

   return

end subroutine dp_puff_solid


subroutine p_eos_gruneisen(i,e_old,p_old)
   implicit none
   integer(4) i,j,k
   real(8) e_old,p_old
   real(8) ea
   real(8) a1,a2,a3,ph,eh

   
       stress_d_11(i)=0.0
	   stress_d_22(i)=0.0
	   stress_d_33(i)=0.0
       stress_d_12(i)=0.0
       stress_d_13(i)=0.0
       stress_d_23(i)=0.0

	   stress_d_123(:,:,i)=0

    
	   ea=e_old-(0.5*p_old+q_av(i))*dvol(i)

	   if(problem_mode==1) then   
           ea=ea+energy_xray(i)
          
       end if
       
	   a1=1+0.5*gama(ne_layer(i))*rho0(ne_layer(i))*dvol(i)
       a2=1-rho0(ne_layer(i))*vol(i)
	   a3=(1-s(ne_layer(i))*a2)**2

	   ph=rho0(ne_layer(i))*c0(ne_layer(i))*c0(ne_layer(i))*a2/a3
	   eh=0.5*ph*(1.0/rho0(ne_layer(i))-vol(i))

	   p(i)=(ph+rho0(ne_layer(i))*gama(ne_layer(i))*(ea-eh))/a1
	   e(i)=ea-0.5*p(i)*dvol(i)       !单位质量物质的内能
	  
	   return
end subroutine p_eos_gruneisen

!------------------------------------------------------------------

subroutine p_eos_puff(i,e_old,p_old)
   implicit none 
   integer(4) i,j,k
   real(8) ea
   real(8) e_old,p_old
   real(8) a0,a1,a2,a3

   stress_d_11(i)=0.0
   stress_d_22(i)=0.0
   stress_d_33(i)=0.0
   stress_d_12(i)=0.0
   stress_d_13(i)=0.0
   stress_d_23(i)=0.0

   stress_d_123(:,:,i)=0

       
    
   ea=e_old-(0.5*p_old+q_av(i))*dvol(i)

   if(problem_mode==1) then   
       ea=ea+energy_xray(i)
       
   end if

  

   a0=1.0-rho0(ne_layer(i))*vol(i)
   a1=h0(ne_layer(i))+(gama(ne_layer(i))-h0(ne_layer(i)))/sqrt(rho0(ne_layer(i))*vol(i))
   a2=es(ne_layer(i))
   a3=cn(ne_layer(i))*rho0(ne_layer(i))*vol(i)*a0

   p(i)=a1*(ea-a2*(1.0-dexp(a3)))/(vol(i)+0.5*a1*dvol(i))
   e(i)=ea-0.5*p(i)*dvol(i) 
   
   return
end subroutine p_eos_puff

!------------------------------------------------------------------

subroutine inner_element_frct(i)
   implicit none
   integer(4) i
   real(8) f_main(3)
   real(8) f_main_vector(3,3)
   real(8)::tol
   tol=0.001

  if(mat_model(ne_layer(i))==0) then

	call eig_calculate(stress_123(:,:,i),3,f_main,tol)
		   if(max(f_main(1),f_main(2),f_main(3))>strs_frct(i)) then
			   fail(i)=1
           end if

   else if(mat_model(ne_layer(i))==1) then
   call eig_calculate(stress_123(:,:,i),3,f_main,tol)
		   if(max(f_main(1),f_main(2),f_main(3))>strs_frct_123(1,1,i)) then
			   fail(i)=1
           end if
	   end if           
  
end subroutine inner_element_frct

subroutine inner_element_frct_paper(i)
implicit none
integer::I
if(stress_123(1,1,i)>tensil_11 .or. stress_123(2,2,i)>tensil_22 .or. stress_123(3,3,i)>tensil_33 .or. stress_123(1,2,i)>tensil_12 .or. stress_123(2,3,i)>tensil_23 .or. stress_123(3,1,i)>tensil_31)then
fail(i)=1
end if


end subroutine inner_element_frct_paper




!----------------------------------------------------------------------

subroutine jacb(i)      !雅各比矩阵     
                                                                
implicit none
integer(4) i,j,k
integer::q

!默认零点l,h,m
dn_l(1)=-0.125*(1+h)*(1+m)
dn_l(2)=-0.125*(1-h)*(1+m)
dn_l(3)=-0.125*(1-h)*(1-m)
dn_l(4)=-0.125*(1+h)*(1-m)
dn_l(5)=0.125*(1+h)*(1+m)
dn_l(6)=0.125*(1-h)*(1+m)
dn_l(7)=0.125*(1-h)*(1-m)
dn_l(8)=0.125*(1+h)*(1-m)

dn_h(1)=0.125*(1-l)*(1+m)
dn_h(2)=-0.125*(1-l)*(1+m)
dn_h(3)=-0.125*(1-l)*(1-m)
dn_h(4)=0.125*(1-l)*(1-m)
dn_h(5)=0.125*(1+l)*(1+m)
dn_h(6)=-0.125*(1+l)*(1+m)
dn_h(7)=-0.125*(1+l)*(1-m)
dn_h(8)=0.125*(1+l)*(1-m)

dn_m(1)=0.125*(1-l)*(1+h)
dn_m(2)=0.125*(1-l)*(1-h)
dn_m(3)=-0.125*(1-l)*(1-h)
dn_m(4)=-0.125*(1-l)*(1+h)
dn_m(5)=0.125*(1+l)*(1+h)
dn_m(6)=0.125*(1+l)*(1-h)
dn_m(7)=-0.125*(1+l)*(1-h)
dn_m(8)=-0.125*(1+l)*(1+h)

do q=1,3
	do j=1,8
	k=ia(j,i)
	coord_xyz(q,j)=xyz(q,k)
	coord_xyz0(q,j)=xyz0(q,k)	                !完全拉格朗日,i号单元对应的8个节点坐标
    end do
end do

  




do j=1,8
	dn_lhm(j,1)=dn_l(j)
	dn_lhm(j,2)=dn_m(j)
	dn_lhm(j,3)=dn_h(j)
end do

      jj=0.0
      jj0=0.0

   jj=matmul(coord_xyz,dn_lhm)
   jj0=matmul(coord_xyz0,dn_lhm)

call solve_inverse(jj,jj_inverse,3)
call solve_inverse(jj0,jj0_inverse,3)


dn_xyz=matmul(dn_lhm,jj_inverse)
!反求dn对x,y,z


return
end subroutine jacb

    

!计算旋转矩阵R
subroutine rotate_matrix(i)!确定旋转的核心问题还是计算三维条件下的变形梯度矩阵先deform_gradient_matrix,that's the point!!!  
integer::i
real(8)::upper(3,3)
call deform_gradient_matrix(i)!根据变形梯度矩阵分解求出旋转矩阵的节奏 

 call QR_decompose(deform,rotate_R,upper,3,3) 

end subroutine rotate_matrix

subroutine original_rotation(i)              !目前看只能用于圆环
integer::i
rotate_original=0
rotate_original(3,3)=1
rotate_original(1,1)=x0(i)/sqrt(x0(i)*x0(i)+y0(i)*y0(i))
rotate_original(2,2)=x0(i)/sqrt(x0(i)*x0(i)+y0(i)*y0(i))
rotate_original(1,2)=-y0(i)/sqrt(x0(i)*x0(i)+y0(i)*y0(i))
rotate_original(2,1)=y0(i)/sqrt(x0(i)*x0(i)+y0(i)*y0(i))

end subroutine original_rotation


! 求整体内力向量(共2nj个分量)

subroutine f_int_total()                                       !由每个单元的应力张量来求对应于每个节点上的内力
   implicit none
   real(8)  f_int_matrix(8,3)
   integer i,j,k
   f_int=0.0
   f_int_matrix=0
   do i=1,ne
       
       if(fail(i)==0) then
		   call jacb(i)
			   f_int_matrix=volume(i)*matmul(dn_xyz,stress_xyz(:,:,i))     !8*J(0,0,0)就是体积
			   do j=1,8
					 k=ia(j,i)
					 f_int(1,k)=f_int(1,k)+f_int_matrix(j,1)
					 f_int(2,k)=f_int(2,k)+f_int_matrix(j,2)
					 f_int(3,k)=f_int(3,k)+f_int_matrix(j,3)
               end do
        end if
   

   
   
   end do

		
   return
end subroutine f_int_total  


!计算沙漏粘性力
subroutine f_hourglass()
implicit none
real::bbeta
real::velocity_hourglass(8,3)!速度
real::T_Hour1(8,1),T_Hour2(8,1),T_Hour3(8,1),T_Hour4(8,1)!四个Hour_glass矢量
real::T_Hour_Matrix(4,8)!四个矢量组成的矩阵
real::H_projection(4,3)
real::F_projection(4,8,3)
real::f_hour(8,3)!最终的力
integer i,j,k



T_Hour1(1,1)=1
T_Hour1(2,1)=-1
T_Hour1(3,1)=1
T_Hour1(4,1)=-1
T_Hour1(5,1)=1
T_Hour1(6,1)=-1
T_Hour1(7,1)=1
T_Hour1(8,1)=-1

T_Hour2(1,1)=1
T_Hour2(2,1)=1
T_Hour2(3,1)=-1
T_Hour2(4,1)=-1
T_Hour2(5,1)=-1
T_Hour2(6,1)=-1
T_Hour2(7,1)=1
T_Hour2(8,1)=1

T_Hour3(1,1)=1
T_Hour3(2,1)=-1
T_Hour3(3,1)=-1
T_Hour3(4,1)=1
T_Hour3(5,1)=-1
T_Hour3(6,1)=1
T_Hour3(7,1)=1
T_Hour3(8,1)=-1

T_Hour4(1,1)=-1
T_Hour4(2,1)=1
T_Hour4(3,1)=-1
T_Hour4(4,1)=1
T_Hour4(5,1)=1
T_Hour4(6,1)=-1
T_Hour4(7,1)=1
T_Hour4(8,1)=-1

	do i=1,8
		T_Hour_matrix(1,i)=T_Hour1(i,1)
		T_Hour_matrix(2,i)=T_Hour2(i,1)
		T_Hour_matrix(3,i)=T_Hour3(i,1)
		T_Hour_matrix(4,i)=T_Hour4(i,1)
	end do
do i=1,ne

		do j=1,8
				  k=ia(j,i)
				  velocity_hourglass(j,1)=v(1,k)
				  velocity_hourglass(j,2)=v(2,k)
				  velocity_hourglass(j,3)=v(3,k)		
		end do
	H_projection=matmul(T_Hour_matrix,velocity_hourglass) !每一次从新计算投影系数
    do j=1,3     !loop x,y,z
		do k=1,4 !loop T_Hour_matrix t1-t4    i:Hour_vector_number; j:xyz_number
			F_projection(k,:,j)=H_projection(k,j)*T_Hour_matrix(k,:)
		end do
	end do

	do j=1,3     !loop x,y,z
		do k=1,4 !loop T_Hour_matrix t1-t4    i:Hour_vector_number; j:xyz_number
			f_hour(:,j)=f_hour(:,j)+F_projection(k,:,j)
		end do
	end do

	 bbeta=qhg*rho(i)*(volume(i)**(2/3))*c0(ne_layer(i))/4.0
	 bbeta=0.00001*bbeta
     f_hour=bbeta*f_hour
		do j=1,8
			k=ia(j,i)
			f_hg(1,k)=f_hg(1,k)+f_hour(j,1)
			f_hg(2,k)=f_hg(2,k)+f_hour(j,2)
			f_hg(3,k)=f_hg(3,k)+f_hour(j,3)		
		end do
		if(gas(i)==1) then
		do j=1,8
			k=ia(j,i)
			f_hg(1,k)=0.0
			f_hg(2,k)=0.0
			f_hg(3,k)=0.0	
		end do
		end if

end do

end subroutine f_hourglass

 

!由于不考虑重力且节点上没有外力，因此外部节点力为0
subroutine f_ext_total()
implicit none
integer j
   do j=1,nj
       f_ext(:,j)=0.0
   end do

end subroutine f_ext_total



!计算变形梯度矩阵
subroutine deform_gradient_matrix(i)!变形梯度矩阵为jj*jj0_inverse------ps:可能涉及到一个行列问题
implicit none
integer::i
call jacb(i)
   deform=matmul(jj,jj0_inverse)

  !我尼玛干嘛不直接imsl矩阵乘法，蛋疼啊我擦 

end subroutine deform_gradient_matrix   

!计算能量：w_int,w_kin,w_ext,w_orgn(初始能量——动能)
subroutine volume_caculate(i)
integer::i,j,k
 core=0
do k=1,3
	do j=1,8
            core(k,1)=core(k,1)+xyz(k,ia(j,i))!求体心坐标
    end do
end do

    core=core/8


	volume(i)=0
	  
	do j=1,6
    call cacul_volume_triple(i,j)
    !计算i单元6个广义面（实际12个面的面积和法向单位矢量）

	call cacul_area_and_vector_surface(i,j)
	volume(i)=volume(i)+volume_s(j)    
	end do


end subroutine volume_caculate 

subroutine cacul_volume_triple(i,j)

   
	real::xyz_local(3,8)    
    real::prism_1(4,4),prism_2(4,4)
    real::vol_1,vol_2
    integer::i,j
	integer::k,t
	integer::order(6,4)
	data order /1,5,1,3,1,2,2,6,2,4,4,3,3,7,6,8,8,7,4,8,5,7,5,6/ 
	do k=1,3
		do t=1,8
	           xyz_local(k,t)=xyz(k,ia(t,i))
		end do
	end do

!到此为止i用完了
        prism_1=1
        prism_2=1
    do k=1,3
        prism_1(k,4)=core(k,1)
        prism_2(k,4)=core(k,1)!prism矩阵第一列为core的x,y,z  
    end do
    
       do t=1,3 !colume
             do k=1,3 !row
                 prism_1(k,t)=xyz_local(k,order(j,t))
                 prism_2(k,t)=prism_1(k,t)  
              end do 
       end do

        do k=1,3
        prism_2(k,2)=xyz_local(k,order(j,4))
        end do 
        
        vol_1=det(prism_1,4)
        vol_2=det(prism_2,4)
		
        volume_s(j)=(abs(vol_1)+abs(vol_2))/6           
      
end subroutine cacul_volume_triple   




subroutine cacul_area_and_vector_surface(i,j)
    integer::i,j                  !源函数引入量，单元编号和面编号
    real::xyz_local(3,8)    
    real::area(2)                 !几号面对应的两个三角面面积
    real::vector_exterior(3,2)    !几号面对应的两个三角面的外法线矢量
    real::vector_judgement(3,1)   !判定内法线还是外法线,用core和A点就好，两个面在A点共用 
    real::vector_AB(3,1),vector_AC(3,1),vector_AD(3,1)  !从A点（左下角）逆时针ABCD
    real::matrix_ABC(3,3),matrix_ACD(3,3)
    integer::k,t                    !循环变量
	integer::order(6,4)
	data order /1,5,1,3,1,2,2,6,2,4,4,3,3,7,6,8,8,7,4,8,5,7,5,6/ 

	do k=1,3
		do t=1,8
	           xyz_local(k,t)=xyz(k,ia(t,i))
		end do
	end do
!读入了编号

	vector_AB=0
	vector_AC=0
	vector_AD=0

	do k=1,3
	   vector_AB(k,1)=xyz_local(k,order(j,2))-xyz_local(k,order(j,1)) !B-A
	   vector_AC(k,1)=xyz_local(k,order(j,3))-xyz_local(k,order(j,1)) !C-A
	   vector_AD(k,1)=xyz_local(k,order(j,4))-xyz_local(k,order(j,1)) !D-A
    end do

	    matrix_ABC=1
        matrix_ACD=1
        vector_exterior=0


    do k=1,3
      matrix_ABC(2,k)=vector_AB(k,1)
      matrix_ABC(3,k)=vector_AC(k,1)
      matrix_ACD(2,k)=vector_AC(k,1)
      matrix_ACD(3,k)=vector_AD(k,1)
    end do
    
	area(1)=det(matrix_ABC,3)/2   !面积就算完了
	area(2)=det(matrix_ACD,3)/2

	vector_exterior(1,1)=matrix_ABC(2,2)*matrix_ABC(3,3)-matrix_ABC(2,3)*matrix_ABC(3,2)
    vector_exterior(2,1)=matrix_ABC(2,3)*matrix_ABC(3,1)-matrix_ABC(2,1)*matrix_ABC(3,3)
	vector_exterior(3,1)=matrix_ABC(2,1)*matrix_ABC(3,2)-matrix_ABC(2,2)*matrix_ABC(3,1)
   
	vector_exterior(1,2)=matrix_ACD(2,2)*matrix_ACD(3,3)-matrix_ACD(2,3)*matrix_ACD(3,2)
    vector_exterior(2,2)=matrix_ACD(2,3)*matrix_ACD(3,1)-matrix_ACD(2,1)*matrix_ACD(3,3)
	vector_exterior(3,2)=matrix_ACD(2,1)*matrix_ACD(3,2)-matrix_ACD(2,2)*matrix_ACD(3,1)   

        !判定方向为外并单位化

    do k=1,3
        vector_judgement(k,1)=xyz_local(k,order(j,1))-core(k,1)!
    end do
       !计算符号并单位化


		if(dot_product(vector_exterior(:,1),vector_judgement(:,1))<0)then        
			 vector_exterior(:,1)=-vector_exterior(:,1)
		else if(dot_product(vector_exterior(:,1),vector_judgement(:,1))>=0)then  !
			 vector_exterior=vector_exterior
		end if
			 vector_exterior(:,1)=vector_exterior(:,1)/sqrt(dot_product(vector_exterior(:,1),vector_exterior(:,1)))


end subroutine cacul_area_and_vector_surface








!------------------------------------------程序中力学运动与反冲运动分界-----------------------------------------------

subroutine compute_blowoff()
implicit none
integer(4) i,j
real(8) blowoff_axis  !对于平板问题
  blowoff_axis=0.0   !对于平板问题
select case(shape)
case(0)
    do i=1,n_cmup
        j=cmup(i)
        if(gas(j)==1 .and. v_cen_x(j)<0) then
		    write(4,*)t,j,v_cen_x(j)
	        blowoff_axis=blowoff_axis+rho0(ne_layer(j))*volume0(j)*abs(v_cen_x(j))
	    end if
	end do
case(1)
		
case default
    stop "shape error."
end select
write(1020,"(f12.7,3x,f12.9)")t,blowoff_axis



end subroutine compute_blowoff


!--------------------------------------X射线辐照部分-------------------------------------------
subroutine xray_init()
   implicit none
   integer(4) i
   real(8) hh

   energy_xray=0.0
   energy_der=0.0
   energy_xray_total=0.0
   wl=0.0
   wer=0.0

  !initiate rad_energy array
   select case(rad_time_mode) !0为矩形谱，1为三角形谱
   case(0)
       rad_energy_num=int(rad_time_width/dt,4)
	   rad_energy=rad_energy_total/rad_time_width !功率密度(能量密度/总时间)  

   case default
       stop "rad_time_mode error."
   end select
   !compute blackbody radiation
   call blackbody_energy()
   !compute mu
   call compute_mu()
   !compute_energy_der
   if(shape==0) then
       call compute_energy_der_plane()
   else if(shape==1) then
	   call compute_energy_der_cylinder()
   else if(shape==2) then
	   call compute_energy_der_semisphere()	
   return

   end if
end subroutine xray_init



subroutine compute_energy_xray()
   implicit none

   select case(rad_time_mode)
   case(0)
       if(shape==0) then
           call compute_rectangle_plane()
		else if(shape==1)then
		   call compute_rectangle_cylinder()
		else if(shape==2)then
		   call compute_rectangle_semisphere()
	   end if   
   case default
       stop "rad_time_mode error"
   end select

   return
end subroutine compute_energy_xray

subroutine mean_coeffcient()


end subroutine mean_coeffcient





subroutine compute_rectangle_plane()!只有单层，远远不够
   implicit none
   integer(4) i,k
  
   k=sn
   if(k<=rad_energy_num) then
       do i=1,ne
	       energy_xray(i)=energy_der(i)*rad_energy(k)*dt*dy*dz/(rho0(ne_layer(i))*volume0(i))  !dt时间内某单元中单位质量的能量沉积
           energy_xray_total(i)=energy_xray_total(i)+energy_xray(i)
       end do
   end if
   
   return

end subroutine compute_rectangle_plane   
subroutine compute_rectangle_cylinder()
!特别要注意采用简化法和六面法在这一步有区别
   implicit none
   integer(4) i,k
  
   k=sn
   if(k<=rad_energy_num) then
       do i=1,ne
	       energy_xray(i)=energy_der(i)*rad_energy(k)*dt/(rho0(ne_layer(i))*volume0(i))  !dt时间内某单元中单位质量的能量沉积
           energy_xray_total(i)=energy_xray_total(i)+energy_xray(i)
       end do
   end if
   return

end subroutine compute_rectangle_cylinder
subroutine compute_rectangle_semisphere()
!特别要注意采用简化法和六面法在这一步有区别
end subroutine compute_rectangle_semisphere


subroutine blackbody_energy()
   implicit none 
   real(8) wl_width,dwl
   real(8) wl_l,wl_r,wl_mid,dh
   integer(4) i
   real(8) sum,w

   wl_width=wl_long-wl_short
   dwl=wl_width/ken   !!将波长的宽度范围wl_width分为ken份，每一小段的宽度即为dwl，在vars里定义了5000份
   sum=0.0
   do i=1,ken
       wl(i)=wl_short+dwl*(i-0.5)  !wl(i)为第i段中间位置对应的波长
       wl_l=wl(i)-0.5*dwl
	   wl_r=wl(i)+0.5*dwl
	   wl_mid=0.5*(wl_l+wl_r)
	   dh=wl_r-wl_l
	   wer(i)=(f_blackbody(wl_l)+4.0*f_blackbody(wl_mid)+f_blackbody(wl_r))*dh/6.0
       sum=sum+wer(i)
	   write(4,*)wl(i),wer(i)/dh,f_blackbody(wl_mid)
   end do

 open(133,file="x_ray_portion.txt") 
 write(133,"('i=',10x,'wavelength=',10x,'wer=',10x,'sum=',10x,'w=')")


   do i=1,ken
       wer(i)=wer(i)/sum    !第i组光子功率占总入射功率的百分比,实际上在这里才完成了单位化
	   w=w+wer(i)

	  
write(133,19)i,wl(i),wer(i),sum,w

19 format (I5,7x,f8.6,7x,f8.6,5x,e8.2,3x,f8.6)

   end do
      
   return



end subroutine blackbody_energy

!----------------------------------------------------------------------
function f_blackbody(wavelen)
implicit none
real(8) wavelen
real(8) f_blackbody  !复合谱的强度

real(8) :: c1=3.7417749e-16
real(8) :: c2=0.01438769
real(8) :: k=1.380658e-23
real(8) :: e=1.60217733e-19

real(8)::f_blackbody_each(20)  !各黑体组分的普朗克函数
real(8)::a_planck(20)    !各黑体谱辐射强度的比例系数
real(8)::t_k(20)       !各黑体组分的黑体温度，单位为k(开)

real(8) wl_m
integer(4) i

   wl_m=wavelen*1.0e-10 !化为国际单位的波长m  
   f_blackbody=0.0
   a_planck(1)=1.0
  
   do i=1,n_blackbody   
       t_k(i)=t_kev(i)*e*1.0e3/k  !将能量转化为温度K
       a_planck(i)=(b_energy(i)*t_k(1)**4)/(b_energy(1)*t_k(i)**4)
	   f_blackbody_each(i)=c1/(wl_m**5.0*dexp(c2/(wl_m*t_k(i))-1.0)) !能谱函数
       f_blackbody=f_blackbody+a_planck(i)*f_blackbody_each(i)

   end do   
   return
end function f_blackbody
!复合谱是有问题的


!----------------------------------------------------------------------

subroutine compute_mu() !计算质量衰减系数
   implicit none
   integer(4) i,j,k
   real(8) mu_oe,mu_kn,mu_a,mu_b
   real(8) wavelen
   real(8) mean_coeff
   integer(4) zz
 mean_coeff=0.0
 open(134,file="mean_coefficient.txt")
 write(134,"('i=',10x,'wavelength=',10x,'coefficient',10x,'wer',10x,'mean_coefficient')") 
   do k=1,ken
       wavelen=wl(k)
	   do i=1,lmax
	       mu_b=0.0
	       do j=1,xlayer(i)%element_num
               zz=xlayer(i)%element_z(j)  !原子序数 
		       call compute_mu_oe(mu_oe,zz,wavelen)
		       call compute_mu_kn(mu_kn,zz,wavelen)
		       mu_a=mu_oe+mu_kn
		       mu_b=mu_b+mu_a*xlayer(i)%element_wp(j)  !对某种波长光子的质量衰减系数为所有元素质量衰减系数的质量加权平均值
          end do 

           	xlayer(i)%mu(k)=mu_b
	    end do	
	end do

  do k=1,ken
		do i=1,lmax
			mean_coeff=mean_coeff+xlayer(i)%mu(k)*wer(k)
			write(134,20)k,wl(k),xlayer(i)%mu(k),wer(k),mean_coeff
			 20 format (I5,7x,f10.3,7x,f10.3,5x,f8.6,3x,f10.2)
		end do
end do


	return
end subroutine compute_mu

!-----------------------------------------------------------------

subroutine compute_mu_oe(mu_oe,zz,wavelen)
   use mu_oe_data
   implicit none 
   real(8) mu_oe
   integer(4) zz
   real(8) wavelen
   real(8) :: planck_h=6.6260755e-34
   real(8) :: me=0.91093897e-30		
   real(8) :: light_c=299792458.0
   real(8) :: k=1.380658e-23
   real(8) :: e=1.60217733e-19
		
   real(8) wl_m, grade_kev
   real(8) alpha, xx0
		
   integer(4) i, i1
   real(8) temp
   
   wl_m=wavelen*1.0e-10
   grade_kev=(planck_h*light_c/wl_m)/(e*1.0e3)  !某波长光子的能量，(h*u的能量量级)换算成单位为Kev
                         
   alpha=planck_h/(me*light_c*wl_m)
   xx0=dlog(511.0*alpha)

   do i=13,1,-1  !与电子壳层相关
       if(grade_kev >= te(i, zz)) then		!zz元素的x射线第i壳层电子的光子吸收限能量，i=13为k层		
				i1=i
				do while(te(i1, zz) < 0.001)
					i1=i1+1
				end do

				temp &
				=whmc(1, i1, zz)+whmc(2, i1, zz)*xx0+whmc(3, i1, zz)*xx0**2.0 &   !whmc为拟合系数 
				+whmc(4, i1, zz)*xx0**3.0

				mu_oe=dexp(temp)*0.602252/aa(zz)
				return			
			end if
		end do

		i1=1
		do while(te(i1, zz) < 0.001)
			i1=i1+1
		end do

		temp &
		=whmc(1, i1, zz)+whmc(2, i1, zz)*xx0+whmc(3, i1, zz)*xx0**2.0 &
		+whmc(4, i1, zz)*xx0**3.0

		mu_oe=dexp(temp)*0.602252/aa(zz)     !光电效应质量衰减系数
		return
	end subroutine compute_mu_oe
!------------------------------------------------------------------

subroutine compute_mu_kn(mu_kn,zz,wavelen)
   use mu_oe_data
   implicit none
   real(8) mu_kn
   integer(4) zz
   real(8) wavelen

   real(8) :: planck_h=6.6260755e-34
   real(8) :: me=0.91093897e-30		
   real(8) :: light_c=299792458.0
   real(8) :: n0=6.0221367e23
   real(8) :: pi=3.1415926
   real(8) :: r0=2.817938e-15

   real(8) wl_m, alpha
   real(8) m0, x1, x2, x3, sigma_e

   wl_m=wavelen*1.0e-10
   alpha=planck_h/(me*light_c*wl_m)

   m0=2.0*pi*(r0*1.0e2)**2.0
   x1=((1.0+alpha)/alpha**2.0) &
   *(2.0*(1.0+alpha)/(1.0+2.0*alpha)-dlog(1.0+2.0*alpha)/alpha)
   x2=dlog(1.0+2.0*alpha)/(2.0*alpha)
   x3=(1.0+3.0*alpha)/(1.0+2.0*alpha)**2.0

   sigma_e=m0*(x1+x2+x3)

   mu_kn=n0*zz*sigma_e/aa(zz)

   return
end subroutine compute_mu_kn
!------------------------------------------------------------------
!单位面积单位时间单位能量对于这种规则的正平面单元的能量沉积
subroutine compute_energy_der_plane()
   implicit none 
   integer(4) i,j,k
   real(8) exp_layer,llx
   real(8) fin,fout,gr

 open(134,file="gr.txt") 
 write(134,"('i=',10x,'gr=')")

   do i=1,ne
       fin=0.0   !某单元左侧光通量
	   fout=0.0   !某单元右侧光通量 （对于入射x射线垂直于单元左右表面的情况）
	   do k=1,ken
	       if(ne_layer(i)==1) then
	           gr=xlayer(ne_layer(i))%mu(k)*rho0(ne_layer(i))
	           fin=fin+wer(k)*dexp(-gr*abs((x0(i)-0.5*dx)))           
		       fout=fout+wer(k)*dexp(-gr*(x0(i)+0.5*dx))

if(i==1)then
write(134,29)k,gr
29 format (I5,7x,f12.4)
end if
   
		   else
		       exp_layer=0.0
			   llx=0.0
		       do j=1,ne_layer(i)-1
                   gr=xlayer(j)%mu(k)*rho0(j)
				   exp_layer=exp_layer+(-gr*lx(j))
				   llx=llx+lx(j)
			   end do
			   gr=xlayer(ne_layer(i))%mu(k)*rho0(ne_layer(i))
			   fin=fin+wer(k)*dexp(exp_layer-gr*abs((x0(i)-0.5*dx-llx)))
               fout=fout+wer(k)*dexp(exp_layer-gr*(x0(i)+0.5*dx-llx))
           end if
	   end do
	   energy_der(i)=fin-fout  !x,x+dx处的能量差值(单位能通量时)

	end do

	return
end subroutine compute_energy_der_plane

subroutine compute_energy_der_cylinder()
 implicit none 
   integer:: i,j                            !单元编号，面编号
   integer::k,t !循环量
   real::xyz_local(3,8)
   real::core(3,1)
   integer::order(6,4)                      !记录器
   real::area                    !几号面对应的两个三角面面积
   real::area_center(3,1)
   real::L_1,L_2,L_3,p
   real::vector_exterior(3,1)    !几号面对应的两个三角面的外法线矢量
   real::vector_judgement(3,1)   !判定内法线还是外法线,用core和A点就好，两个面在A点共用 
   real::vector_std(3,1)  
   real::vector_AB(3,1),vector_AC(3,1),vector_AD(3,1)  !从A点（左下角）逆时针ABCD
   real::matrix_ABC(3,3),matrix_ACD(3,3)
   real::deposite_surface_dken
   real::d_surface
   real::fin_surface(6)
   real::gr
   
   data order /1,5,1,3,1,2,2,6,2,4,4,3,3,7,6,8,8,7,4,8,5,7,5,6/ 
   vector_std=0
   vector_std(1,1)=1
   do i=1,ne
      if(x0(i)<0) then 
           !i号单元内的能量
           fin_surface=0.0   !六个面的能量
           xyz_local=0 
           do k=1,3
				do t=1,8
	           xyz_local(k,t)=xyz(k,ia(t,i))
				end do
		   end do
           
                core=0
	    do k=1,3
           do t=1,8
               core(k,1)=core(k,1)+xyz_local(k,t)
           end do
        end do
             
        core=core/8

       do j=1,6
        vector_AB=0
	    vector_AC=0
	    vector_AD=0  
			do k=1,3
	   			vector_AB(k,1)=xyz_local(k,order(j,2))-xyz_local(k,order(j,1)) !B->A
	   			vector_AC(k,1)=xyz_local(k,order(j,3))-xyz_local(k,order(j,1)) !C->A
	  			vector_AD(k,1)=xyz_local(k,order(j,4))-xyz_local(k,order(j,1)) !D->A
       		end do      
	     matrix_ABC=1
         matrix_ACD=1
         vector_exterior=0


       do k=1,3
             matrix_ABC(2,k)=vector_AB(k,1)
             matrix_ABC(3,k)=vector_AC(k,1)
             matrix_ACD(2,k)=vector_AC(k,1)
             matrix_ACD(3,k)=vector_AD(k,1)
       end do
	   L_1=0
	   L_2=0
	   L_3=0
          call Lengh_caculate(vector_AB,L_1)
          call Lengh_caculate(vector_AC,L_2)
          call Lengh_caculate(vector_AB-vector_AC,L_3)
		  p=(L_1+L_2+L_3)/2
	      area=sqrt(p*(p-L_1)*(p-L_2)*(p-L_3))*2

       do k=1,3
         area_center(k,1)=xyz_local(k,order(j,1))+xyz_local(k,order(j,2))+xyz_local(k,order(j,3))+xyz_local(k,order(j,4))
       end do

        area_center=area_center/4

      	vector_exterior(1,1)=matrix_ABC(2,2)*matrix_ABC(3,3)-matrix_ABC(2,3)*matrix_ABC(3,2)
        vector_exterior(2,1)=matrix_ABC(2,3)*matrix_ABC(3,1)-matrix_ABC(2,1)*matrix_ABC(3,3)
    	vector_exterior(3,1)=matrix_ABC(2,1)*matrix_ABC(3,2)-matrix_ABC(2,2)*matrix_ABC(3,1)
  	

        !判定方向为外并单位化

       do k=1,3
          vector_judgement(k,1)=xyz_local(k,order(j,1))-core(k,1)!
       end do
       !计算符号并单位化

       if(dot_product(vector_exterior(:,1),vector_judgement(:,1))<0)then        
	  	 vector_exterior(:,1)=-vector_exterior(:,1)
	   else if(dot_product(vector_exterior(:,1),vector_judgement(:,1))>=0)then  !
	  	 vector_exterior=vector_exterior
	   end if
         vector_exterior(:,1)=vector_exterior(:,1)/sqrt(dot_product(vector_exterior(:,1),vector_exterior(:,1)))
       d_surface=-sqrt(abs(r_out*r_out-area_center(2,1)*area_center(2,1)))-area_center(1,1)


     
		   do k=1,ken
		   	   deposite_surface_dken=0
	           if(ne_layer(i)==1) then    !第一层材料
	               gr=xlayer(ne_layer(i))%mu(k)*rho0(ne_layer(i))                   		      
				   deposite_surface_dken=exp(-gr*abs(d_surface))*area*dot_product(vector_exterior(:,1),vector_std(:,1)) 
 
				   !每一份光子，ken循环，所以要累加			
                                   !我要算6个面的然后求和即可，既有光子层又有6个面层                         
	               fin_surface(j)=fin_surface(j)+wer(k)*deposite_surface_dken
			   end if		   
	       end do
           end do
            do j=1,6
	              energy_der(i)=energy_der(i)+fin_surface(j)  !四边形网格内的能量沉积(对于单位能通量)\
            end do
              energy_der(i)=-energy_der(i) 
           write(3,*)i
        else
     	energy_der(i)=0        
       end if


  end do
   write(*,"('energy in grid is computed ')") 
  print*,energy_der(1165) !1147,1165,1185
  print*,energy_der(1544) !1595,1544,1548

return


end subroutine compute_energy_der_cylinder

	subroutine Lengh_caculate(vector,L)
		implicit none
		real::vector(3,1)
		real::L
        L=sqrt(vector(1,1)*vector(1,1)+vector(2,1)*vector(2,1)+vector(3,1)*vector(3,1))
		return
	end subroutine Lengh_caculate







subroutine compute_energy_der_semisphere()

end subroutine compute_energy_der_semisphere







subroutine energy_compute()
  implicit none
   integer i,j,k

   ee=0.0
   w_kin=0.0
   w_orgn=0.0

   do i=1,ne
       ee=ee+rho0(ne_layer(i))*volume0(i)*e(i)
   end do

   do j=1,nj
	   do k=1,3
       w_ext=w_ext+0.5*dt*v_mid(k,j)*(f_ext_old(k,j)+f_ext(k,j))
       w_kin=w_kin+0.5*mass(j)*(v(k,j)**2)
	   end do
   end do

!   do i=1,ne
!       w_kin=w_kin+0.5*rho0*area0(i)*((v_cen_x(i))**2+v_cen_y(i)**2)
!   end do

   if(problem_mode==0) then
       w_orgn=0.5*rho0(1)*(0.4*length_y*length_z)*(vf**2)	 !第一层材料的动能
   else if(problem_mode==1) then
       do i=1,ne
	       w_orgn=w_orgn+energy_xray_total(i)*rho0(ne_layer(i))*volume0(i)
       end do
   end if 
   
   w_err=(w_orgn+w_ext-(ee+w_kin))/w_orgn*100
  
end subroutine energy_compute  
 


!----------------------------------------------------------------------------------------------x_ray




!-----------------------------------------------------------------------------------------------out_put_subs
subroutine out_put()          
implicit none
integer(4) i,j,k
integer(4) n_ne    !与每个节点相关的单元数
integer(8)::no_fail_num
character(len=25) layer_output(lmax) 
   open(150,file="w_eff.txt")   
   write(150,10)t,w_err,w_orgn,ee,w_kin
   10 format(3x,"t=",f9.6,8x,"w_err=",f10.6,8x,"w_orgn=",f10.6,8x,"ee=",f10.6,8x,"w_kin=",f10.6)
   
   write(*,"('t=',f10.6)")t
     	no_fail_num=ne
		do i=1,ne
			if(fail(i)==1) then
				no_fail_num=no_fail_num-1
			end if
		end do
open(901,file="failure_condtion.txt")
write(901,91)t,ne-no_fail_num,(ne-no_fail_num)/ne	
91 format(f8.6,5x,8I,5x,f8.6)
!按总体节点输出

   do i=1,num_time_record

       if(abs(sn*dt-time_record(i))<=dt/3.0) then 
	   
           open(200+i,file=time_output(i))
		   write(200+i,"(A)"),'TITLE="EXAMPLE"'
           write(200+i,"(A)"),'VARIABLES="X","Y","Z","stress_node_xx","stress_node_yy","stress_node_ZZ","stress_node_XY","stress_node_XZ","stress_node_YZ","P","E"'
		   write(200+i,"(A)"),'ZONE T="p_1"'
		   write(200+i,"('N=',i5,3x,'E=',i5)")nj,no_fail_num
		   write(200+i,"(A)"),'F=FEPOINT,ET="BRICK"'      
		end if

   end do
 

   do i=1,ne     
       index_ele_write(i)=0
	   if(shape==0 .and. fail(i)==0) then
	       index_ele_write(i)=1    
       else if (shape==1 .and. fail(i)==0) then      !改半圆就不是个问题
	       index_ele_write(i)=1  !单元输出控制指数，index_node_write(j)=1该单元处于上半平面需要输出相应结果
	   else
	       index_ele_write(i)=0
	   end if     
   end do

   do j=1,nj
       stress_node_xx(j)=0
	   stress_node_yy(j)=0
	   stress_node_zz(j)=0
	   stress_node_xy(j)=0
	   stress_node_xz(j)=0
	   stress_node_yz(j)=0
	   p_node(j)=0
	   e_node(j)=0
       n_ne=0	   
       do i=1,ne
	       do k=1,8
		       if (ia(k,i)==j) then
			       n_ne=n_ne+1
				   stress_node_xx(j)=stress_node_xx(j)+stress_xyz(1,1,i)  !线性插值
                   stress_node_yy(j)=stress_node_yy(j)+stress_xyz(2,2,i) 
				   stress_node_zz(j)=stress_node_zz(j)+stress_xyz(3,3,i) 
				   stress_node_xy(j)=stress_node_xy(j)+stress_xyz(1,2,i)
				   stress_node_xz(j)=stress_node_xz(j)+stress_xyz(1,3,i)
				   stress_node_yz(j)=stress_node_yz(j)+stress_xyz(2,3,i)
				   p_node(j)=p_node(j)+p(i) 
				   e_node(j)=e_node(j)+e(i)                  
			   end if
		   end do	  
	   end do


	   stress_node_xx(j)=stress_node_xx(j)/n_ne      
       stress_node_yy(j)=stress_node_yy(j)/n_ne
	   stress_node_zz(j)=stress_node_zz(j)/n_ne
	   stress_node_xy(j)=stress_node_xy(j)/n_ne
	   stress_node_xz(j)=stress_node_xz(j)/n_ne
	   stress_node_yz(j)=stress_node_yz(j)/n_ne 
	   p_node(j)=p_node(j)/n_ne
	   e_node(j)=e_node(j)/n_ne  
       
       do i=1,num_time_record
		   if(abs(sn*dt-time_record(i))<=dt/3.0) then
	           write(200+i,20)xyz(1,j),xyz(2,j),xyz(3,j),stress_node_xx(j)*100,stress_node_yy(j)*100,stress_node_zz(j)*100,stress_node_xy(j)*100,stress_node_xz(j)*100,stress_node_yz(j)*100,p_node(j)*100,e_node(j)
		       20 format(11(f10.6,3x))
		   end if
	   end do
   
   end do
  
   do i=1,ne
       if(index_ele_write(i)==1) then
	        do j=1,num_time_record
		       if(abs(sn*dt-time_record(j))<=dt/3.0) then
	               write(200+j,25)ia(1,i),ia(2,i),ia(3,i),ia(4,i),ia(5,i),ia(6,i),ia(7,i),ia(8,i)
                   25 format(8(i6,3x))
			   end if
			end do
	   end if
   end do
   return
end subroutine out_put

subroutine out_put_element()
implicit none
integer i
   if(shape==0) then
       open(400,file="out_ele_up.txt")
       write(400,"('t=',f10.6)")t
       write(400,"(4x,'x=',12x,'y=',12x,'z=',12x,'v_cen_x',7x,'v_cen_y',5x,'v_cen_z',5x,  &
           'stress_xx',4x,'stress_yy',4x,'stress_zz',4x,'stress_xy',4x,'stress_xz',4x,'stress_yz',&
		    9x,'p=',9x,'q=',9x,'e=',5x,   &
		    'f1_yield',8x,'strain_xx',4x,'strain_yy',4x,'strain_zz',4x,'strain_xy',4x,'strain_xy',4x,'strain_yz',5x,&
			'rate',5x, 'plastic_index',9x,'yield_mode',5x,'gas',5x,'fail')") 

       do i=1,n_cmup         
          write(400,40) x(cmup(i)),y(cmup(i)),z(cmup(i)),v_cen_x(cmup(i)),v_cen_y(cmup(i)),v_cen_z(cmup(i)),        &
                  stress_xyz(1,1,cmup(i)),stress_xyz(2,2,cmup(i)),stress_xyz(3,3,cmup(i)),stress_xyz(1,2,cmup(i)),stress_xyz(1,3,cmup(i)),stress_xyz(2,3,cmup(i)), &
				  p(cmup(i)),q(cmup(i)),e(cmup(i)),    &
				  f1_yield(cmup(i)),strain_xyz(1,1,cmup(i)),strain_xyz(2,2,cmup(i)),strain_xyz(3,3,cmup(i)),strain_xyz(1,2,cmup(i)),strain_xyz(1,3,cmup(i)),strain_xyz(2,3,cmup(i)),  &
				  rate(cmup(i)),plastic_index(cmup(i)),yield_mode(ne_layer(cmup(i))),gas(cmup(i)),fail(cmup(i))
          40 format(23(f12.7,3x),4(i3,2x))
!          write(401,"(7(f12.7,3x))")strain_xx(cmup(i)),strain_yy(cmup(i)),strain_xx(cmup(i))+strain_yy(cmup(i)),        &
!		          strain_11(cmup(i)),strain_22(cmup(i)),strain_11(cmup(i))+strain_22(cmup(i)),(vol(cmup(i))-vol0(cmup(i)))/vol(cmup(i))
       end do
	end if
end subroutine out_put_element

!------------------------------------------------------------------------------------------------------文件添名字
subroutine  output_filename_definition
implicit none
integer(4) i
character(len=2) temp


   do i=1,num_time_record
       select case(i)
	   case(0 : 9)
	       call convert_char(temp, '0', i-0)
	   case(10 : 19)
		   call convert_char(temp, '1', i-10)	  
       case(20)
		   call convert_char(temp, '2', i-20)
	   case default
		   stop "Global variable num_time_record> max_num_time_record_20,    &
				      program terminated."
	   end select

	   time_output(i)="time_"//temp//"_output.txt" 
   end do

   do i=1,num_elefix_x
       select case(i)
	   case(0 : 9)
	       call convert_char(temp, '0', i-0)
	   case(10 : 19)
		   call convert_char(temp, '1', i-10)	  
       case(20)
		   call convert_char(temp, '2', i-20)
	   case default
		   stop "Global variable num_time_record> max_num_time_record_20,    &
				      program terminated."
	   end select

	   fne_his_x(i)="fix_ex_"//temp//".dat" 
   end do

   do i=1,num_elefix_y
       select case(i)
	   case(0 : 9)
	       call convert_char(temp, '0', i-0)
	   case(10 : 19)
		   call convert_char(temp, '1', i-10)	  
       case(20)
		   call convert_char(temp, '2', i-20)
	   case default
		   stop "Global variable num_time_record> max_num_time_record_20,    &
				      program terminated."
	   end select

	   fne_his_y(i)="fix_ey_"//temp//".dat" 
   end do


   do i=1,num_elefix_z
       select case(i)
	   case(0 : 9)
	       call convert_char(temp, '0', i-0)
	   case(10 : 19)
		   call convert_char(temp, '1', i-10)	  
       case(20)
		   call convert_char(temp, '2', i-20)
	   case default
		   stop "Global variable num_time_record> max_num_time_record_20,    &
				      program terminated."
	   end select

	   fne_his_z(i)="fix_ez_"//temp//".dat" 
   end do
 end subroutine output_filename_definition
 
 
 
 
 	subroutine convert_char(temp, ch, n)
		implicit none
		character(len=2) temp
		character(len=1) ch
		integer(4) n
		
		temp=ch//achar(iachar('0')+n)

		return
	end subroutine convert_char
!-------------------------------------------------------------------------------------------------------------------------------
subroutine out_put_energy
implicit none
integer(4) i,k

k=sn
if(k<=rad_energy_num) then  
    
    if(shape==0) then
	    open(1010,file="energy_ele_up.txt")
		write(1010,"('t=',f10.6)")t
		write(1010,"('x=',10x,'y=',10x,'z=',10x,'e=',7x,'e_deposited=')") 

	    do i=1,n_cmup         
            write(1010,"(5(f12.7,3x))") x(cmup(i)),y(cmup(i)),z(cmup(i)),e(cmup(i)),energy_xray_total(cmup(i))
        end do   
	end if    
end if
   
return
end subroutine out_put_energy

subroutine out_put_history()
implicit none
integer(4) i,j
real(8) t_axix(500000)
real(8) peak_stress_xx,peak_x,peak_y,peak_z
real(8) peak_axis(500000)
real(8) peak_axis_stress_xx,peak_axis_x,peak_axis_y,peak_axis_z
integer(4) peak_i,peak_axis_i


do i=1,num_elefix_x
   j=ne_his_x(i)
   
   open(600+i,file=fne_his_x(i))
   write(600+i,60)t,x(j),y(j),z(j),v_cen_x(j),v_cen_y(j),v_cen_z(j),strain_xyz(1,1,j),strain_xyz(2,2,j),strain_xyz(3,3,j),strain_xyz(1,2,j),strain_xyz(1,3,j),strain_xyz(2,3,j),       &
                         stress_xyz(1,1,j),stress_xyz(2,2,j),stress_xyz(3,3,j),stress_xyz(1,2,j),stress_xyz(1,3,j),stress_xyz(2,3,j),p(j),e(j)
   60 format(21(f12.7,2x))

end do

do i=1,num_elefix_y
   j=ne_his_y(i)
   open(700+i,file=fne_his_y(i)) 
   write(700+i,60)t,x(j),y(j),z(j),v_cen_x(j),v_cen_y(j),v_cen_z(j),strain_xyz(1,1,j),strain_xyz(2,2,j),strain_xyz(3,3,j),strain_xyz(1,2,j),strain_xyz(1,3,j),strain_xyz(2,3,j),       &
                         stress_xyz(1,1,j),stress_xyz(2,2,j),stress_xyz(3,3,j),stress_xyz(1,2,j),stress_xyz(1,3,j),stress_xyz(2,3,j),p(j),e(j)
   70 format(21(f12.7,2x))
 
end do

 do i=1,num_elefix_z
   j=ne_his_z(i)
   open(800+i,file=fne_his_z(i))  
   write(800+i,60)t,x(j),y(j),z(j),v_cen_x(j),v_cen_y(j),v_cen_z(j),strain_xyz(1,1,j),strain_xyz(2,2,j),strain_xyz(3,3,j),strain_xyz(1,2,j),strain_xyz(1,3,j),strain_xyz(2,3,j),       &
                         stress_xyz(1,1,j),stress_xyz(2,2,j),stress_xyz(3,3,j),stress_xyz(1,2,j),stress_xyz(1,3,j),stress_xyz(2,3,j),p(j),e(j)
   80 format(21(f12.7,2x))
                    
end do


  ! if(problem_mode==1) then
       if(shape==0) then
           do i=1,n_cmup
	           j=cmup(i)
	           if(stress_xx(j)<peak_axis(j)) then
		           peak_axis(j)=stress_xyz(1,1,j)
				   t_axix(j)=t
		       end if

		       if(abs(t-t_total)<dt/3.0) then
			       write(950,"(2x,i7,4x,f8.6,4x,f8.6,4x,f8.6,4x,f8.6)")j,x0(j),y0(j),z0(j),peak_axis(j),t_axix(j)
		       end if
	       end do
       end if
 !  end if

 !  if(problem_mode==1) then
       if(shape==0) then       
           peak_axis_stress_xx=stress_xyz(1,1,cmup(1))
	       peak_axis_i=cmup(1)
	       do i=2,n_cmup
	           if(stress_xyz(1,1,cmup(i))<peak_axis_stress_xx) then
		           peak_axis_stress_xx=stress_xyz(1,1,cmup(i))
			       peak_axis_i=cmup(i)
			       peak_axis_x=x0(cmup(i))
			       peak_axis_y=y0(cmup(i))
				   peak_axis_z=z0(cmup(i))			      
		       end if
	       end do
		   write(900,"(f8.6,4x,i7,4x,4(f12.7,3x))")t,peak_axis_i,peak_axis_x,peak_axis_y,peak_axis_z,peak_axis_stress_xx	      
	   end if 	  	       
 !  end if   

 !  if(problem_mode==1) then
       
       peak_stress_xx=stress_xyz(1,1,cmup(1))
	   peak_i=1
	   do i=2,ne
	       if(stress_xyz(1,1,i)<peak_stress_xx) then
		       peak_stress_xx=stress_xx(i)
			   peak_i=i
			   peak_x=x0(i)
			   peak_y=y0(i)
			   peak_z=z0(i)
		   end if
	   end do
		   write(1000,"(f8.6,4x,i7,4x,4(f12.7,3x))")t,peak_i,peak_x,peak_y,peak_z,peak_stress_xx
  ! end if

   return


end subroutine out_put_history


subroutine allocate_mode()

	print*,"number of element is",ne
	print*,"number of node is",nj

allocate(fail(ne))
allocate(ne_layer(ne))
allocate(ia(9,ne))
allocate(index_ele_write(ne))

allocate(rho(ne),rho_old(ne))
allocate(volume0(ne),volume(ne),volume_old(ne),dvolume(ne),volume_average(ne))

allocate(vol0(ne),vol(ne),vol_old(ne),dvol(ne),vol_average(ne))


allocate(e(ne),de(ne))
allocate(p(ne),dp(ne))

allocate(strs_yld_11(ne),strs_yld_22(ne),strs_yld_33(ne),strs_yld_12(ne),strs_yld_13(ne),strs_yld_23(ne))
allocate(strs_yld_123(3,3,ne))

allocate(rate(ne))
allocate(strain_rate_eff(ne))

allocate(energy_der(ne))
allocate(energy_xray_total(ne))
allocate(energy_xray(ne))

allocate(c11(ne),c22(ne),c33(ne),c12(ne),c21(ne),c13(ne),c31(ne),c23(ne),c32(ne),c44(ne),c55(ne),c66(ne))

!allocate(c123(6,6,ne))
allocate(ex(ne),ey(ne),ez(ne),gxy(ne),gxz(ne),gyz(ne),pxy(ne),pxz(ne),pyz(ne))

allocate(strs_yld(ne))!各向同性材料的屈服强度,前者为以层表示，后者以单元表示
allocate(f1_yield(ne))

allocate(strs_frct(ne))
allocate(strs_frct_123(3,3,ne))

!--------------------------------
allocate(strs_frct_11(ne),strs_frct_22(ne),strs_frct_33(ne),strs_frct_12(ne),strs_frct_13(ne),strs_frct_23(ne))
allocate(plastic_index(ne))
allocate(gas(ne))
allocate(stress_xx(ne),stress_yy(ne),stress_zz(ne),stress_xy(ne),stress_xz(ne),stress_yz(ne) )
allocate(stress_xyz(3,3,ne))
allocate(stress_11(ne),stress_22(ne),stress_33(ne),stress_12(ne),stress_13(ne),stress_23(ne)) 
allocate(stress_123(3,3,ne))
allocate(stress_d_11(ne),stress_d_22(ne),stress_d_33(ne),stress_d_12(ne),stress_d_13(ne),stress_d_23(ne))
allocate(stress_d_123(3,3,ne))
allocate(q(ne),q_old(ne),q_av(ne)) 


allocate(strain_xx(ne),strain_yy(ne),strain_zz(ne),strain_xy(ne),strain_xz(ne),strain_yz(ne))
allocate(strain_xyz(3,3,ne))

allocate(strain_xx_rate(ne),strain_yy_rate(ne),strain_zz_rate(ne),strain_xy_rate(ne),strain_xz_rate(ne),strain_yz_rate(ne))   !应变率分量
allocate(strain_xyz_rate(3,3,ne))
allocate(strain_11(ne),strain_22(ne),strain_33(ne),strain_12(ne),strain_13(ne),strain_23(ne))      
allocate(strain_123(3,3,ne))
allocate(strain_11_rate(ne),strain_22_rate(ne),strain_33_rate(ne),strain_12_rate(ne),strain_13_rate(ne),strain_23_rate(ne))
allocate(strain_p_11(ne),strain_p_22(ne),strain_p_33(ne),strain_p_12(ne),strain_p_13(ne),strain_p_23(ne))
allocate(strain_p_123(3,3,ne))
allocate(strain_p_eff(ne),strain_p_rate_eff(ne))
!-------------------------

allocate(xyz0(3,nj),xyz(3,nj))
allocate(stress_node_xx(nj),stress_node_yy(nj),stress_node_zz(nj),stress_node_xy(nj),stress_node_xz(nj),stress_node_yz(nj),p_node(nj),e_node(nj)) 
allocate(u(3,nj),du(3,nj),v(3,nj),a(3,nj),v_mid(3,nj))  
allocate(mass(nj))     !在:个自由度上的质量矩阵（对角化）
allocate(x0(nj),y0(nj),z0(nj),x(nj),y(nj),z(nj))   !初始时每单元中心的整体坐标
allocate(v_cen_x(nj),v_cen_y(nj),v_cen_z(nj))
allocate(f_ext(3,nj),f_int(3,nj),f_ext_old(3,nj),f_int_old(3,nj))
allocate(hggama1(nj),hggama2(nj),hggama3(nj),hggama4(nj))  !每个单元沙漏基矢量γ的四个分量（三维情况下四个是否够？？）
allocate(f_hg(3,nj))

print*,"allocation is complete"




end subroutine allocate_mode



end module subs


