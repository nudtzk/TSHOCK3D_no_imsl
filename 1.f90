module vars                            !应该除了沙漏部分需要看看其他都定义完成
implicit none
!namelist:compute_control-----------------------------------------------------------------------------------------------------------
character(len=50),save:: problem
integer(4), save :: problem_mode
integer(4),save:: shape                 !0代表平板形状，1代表圆柱壳体形状
integer(4),save:: eos_mode  
real(8),save:: flyer_length
integer(4),save::lmax,lmax_iso,lmax_aniso     !各向同性材料与各向异性材料的层数  
real(8),save:: length_1,length_2,length_3
real(8),save::dl_1,dl_2,dl_3  
real(8),save::q1,q2,qhg    !粘性系数!沙漏粘性系数
real(8),save:: t,dt,t_total   !时间参数
integer(4),save::nscr
integer(8),save::ne,nj
real(8),save::vf
integer(4),save:: num_elefix_x,num_elefix_y,num_elefix_z  !选取x/y/z坐标相同的一组点最多各取20个
integer(4),save:: ne_his_x(20),ne_his_y(20),ne_his_z(20)  !对应参数
integer(4),save::num_time_record  !指定输出结果的时刻总数,最多=20
real(8),save::time_record(20)     !指定需要输出结果的具体时刻值
!----------------------------------------------------------------------------------------------------------------------------------------

!namelist layer_setting--------------------------------------------------------
integer(4),save::mat_model(30)     !各层材料类型指示，0->各向同性，1->各向异性
integer(4),save::yield_mode(30)
real(8),save::length_layer(30)
!------------------------------------------------------------------------------


!namelist model_0------------------------------------------------------------------------------------------------------------
integer(4),save::layer_no_iso(30)
real(8),save::rho0_iso(30),c0_iso(30),s_iso(30),gama_iso(30)
real(8),save::h0_iso(30),es_iso(30)
real(8),save::ex_iso(30),ey_iso(30),ez_iso(30),pxy_iso(30),pxz_iso(30),pyz_iso(30),gxy_iso(30),gxz_iso(30),gyz_iso(30)
real(8),save::strs_yld_0(30),strs_frct_0(30)
!----------------------------------------------------------------------------------------------------------------------------

!namelist model_0_yield_1----------------------------------------------
real(8),save::ax_iso(30),nx_iso(30),beta0_iso(30),strain_rate0_iso(30)
!----------------------------------------------------------------------

!namelist model_1---------------------------------------------------
integer(4),save::layer_no_aniso(30) !各向异性材料中第i层材料在整体层数中的位置，主要用于数据输入
real(8),save::rho0_aniso(30),c0_aniso(30),s_aniso(30),gama_aniso(30)
real(8),save::h0_aniso(30),es_aniso(30)
real(8),save::ex_aniso(30),ey_aniso(30),ez_aniso(30),pxy_aniso(30),pxz_aniso(30),pyz_aniso(30),gxy_aniso(30),gxz_aniso(30),gyz_aniso(30)
real(8),save::strs_yld_110(30),strs_yld_220(30),strs_yld_330(30),strs_yld_120(30),strs_yld_130(30),strs_yld_230(30)  
real(8),save::strs_frct_110(30),strs_frct_220(30),strs_frct_330(30),strs_frct_120(30),strs_frct_130(30),strs_frct_230(30)
!-------------------------------------------------------------------

!namelist model_1_yield_1-------------------------------------------
real(8),save::ax_aniso(30),ay_aniso(30),az_aniso(30),nx_aniso(30),ny_aniso(30),nz_aniso(30)
 real(8),save::axy_aniso(30),nxy_aniso(30),axz_aniso(30),nxz_aniso(30),ayz_aniso(30),nyz_aniso(30)
real(8),save::beta0_aniso(30),strain_rate0_aniso(30)
!-------------------------------------------------------------------

!namelist xray-------------------------------------------
integer(4),save::n_blackbody    !黑体谱的个数，最多取20组
real(8),save::t_kev(20)         !各黑体组分的黑体温度，单位为keV
real(8),save::b_energy(20)      !各黑体组分所占总能量的能量份额
real(8),save::wl_short,wl_long
integer(4),save::rad_time_mode  !时间谱模式，0为矩形谱，1为三角形谱
real(8),save::rad_time_width    !时间谱的宽度，三角形谱为半高宽
real(8),save::rad_energy_total  !入射通量密度
integer(4),parameter::ken=5000
integer(4),parameter:: max_element_num=10
integer(4),parameter::max_layer=30
type::layer_xray
    integer(4) element_num    !各层的元素总数，<=10
	integer(4) element_z(max_element_num)    !元素的原子序数，元素种类不超过10个
	real(8) element_wp(max_element_num)     !各元素的质量百分比
	real(8) mu(ken)
end type layer_xray
type(layer_xray), save :: xlayer(max_layer)

!------------------------------------------------------

!namelist xray_output-----------------------------------------------------------------------------------------------------------
integer(4),save::cm_r1(5000),cm_r2(5000)  
integer(4),save::cm_theta_0(2000),cm_theta_10(2000),cm_theta_20(2000),cm_theta_30(2000),cm_theta_40(2000),cm_theta_45(2000),   &
                 cm_theta_50(2000),cm_theta_60(2000),cm_theta_70(2000),cm_theta_80(2000),cm_theta_90(2000)
!控制输出部分参数-------------------------------------------------------------------------------------
character(len=18) time_output(20)                                    !基本形式为time_01_output.txt
character(len=13) fne_his_x(20),fne_his_y(20),fne_his_z(20)         !基本形式为fix_ex_01.dat,fix_ey_01.dat
integer(4),save::st,sn  !关于计算步数的控制变量
!基本没用参数，除了dx在算平板会简单----------------------------------------------------------------------------------
real(8),save::length_x,length_y,length_z,lx(30),dx,dy,dz          !对应平面x方向总长度，y方向总长度，z方向总长度各层的厚度，x方向y方向z方向的网格尺寸
real(8),save:: r_in,r_out,r_height,lr(30),dr,d_theta,d_height     !分别为圆环的内半径、外半径，各层沿半径方向的厚度，半径方向的网格尺寸,角度划分尺寸
!-------------------------------------------------------------------------------------
integer(4),save::n_cmup  !中轴线上方的单元数量
integer(4),save::cmup(10000)    !沿x方向的中轴上方所包含的单元(平板模型)

real(8),save::ex_0(30),ey_0(30),ez_0(30),gxy_0(30),gxz_0(30),gyz_0(30)!,pxy(30),pxz(30),pyz(30)!三个弹性模量，剪切模量，泊松比
real(8),save::c11_0(30),c22_0(30),c33_0(30),c12_0(30),c13_0(30),c23_0(30),c44_0(30),c55_0(30),c66_0(30)         !初始各层本构系数，由刚度矩阵转换得到,包括3*3的三个正应力关系和3个剪切关系
real(8),save::c123_0(6,6,20)
real(8),save::rho0(30)  !各层的初始密度

real(8),save::shape_n(8),shape_n0(8) !形状函数,中心点的形状函数
real(8),save::dn_l(8),dn_h(8),dn_m(8),dn0_l(8),dn0_h(8),dn0_m(8)  !形状函数、中心点的形状函数对自然坐标l,h,m的导数
real(8),save::dn_lhm(8,3)
real(8),save::dn_xyz(8,3)
real(8),save::coord_xyz(3,8)
real(8),save::coord_xyz0(3,8)
real(8),save::dn_x(8),dn_y(8),dn_z(8)           !形状函数对整体坐标x,y,z的导数
real(8),save::jj(3,3),jj_inverse(3,3)       !(x,y,z)对(l,h,m)偏导数的雅各比矩阵,矩阵行列式，逆矩阵
real(8),save::jj0(3,3),jj0_inverse(3,3)    !(x0,y0,z0)对(l,h,m)偏导数的雅各比矩阵,矩阵行列式，逆矩阵
real(8),save::f11,f12,f13,f21,f22,f23,f31,f32,f33                  !(x,y,z)对(x0,y0,z0)的偏导数的雅各比矩阵，变形梯度矩阵
real(8),save::deform(3,3)
real(8),save::rotate_R(3,3)
real(8),save::rotate_original(3,3)

!-----------------------------------------------
real(8),save::c0(30),s(30),gama(30)  !物态方程参数
real(8),save::h0(30),es(30)      !PUFF物态方程参数及升华能
real(8),save::ax(30),ay(30),az(30),nx(30),ny(30),nz(30),axy(30),nxy(30),axz(30),nxz(30),ayz(30),nyz(30)
real(8),save::cn(30)
real(8),save::beta0(30),strain_rate0(30)        !应变率有效因子的参数0.02183,初始参考应变率 取值为0.001/s,应变率因子
!---------------------------------
real(8),save::diag_vector(3)
real(8),save::diag_matrix(3,3)
real(8),save::diag_sum
real(8),save::p_dst,p_dsp,ee         !偏应变贡献压力，塑性应变贡献压力,能量项




integer(4),allocatable,save::fail(:)  !材料单元当前失效判定，0代表材料未失效，1代表材料失效
integer(4),allocatable,save::ne_layer(:)  !各单元所属层数
integer(4),allocatable,save::ia(:,:)   !j=1-8表示组成某单元的8个节点编号,j=9代表该单元所在层数的属性
integer(4),allocatable,save::index_ele_write(:)         !对于圆柱壳体模型，index_ele_write(i)==1代表该单元处于上半平面，需要输出结果

real(8),allocatable,save::rho(:),rho_old(:),volume0(:),volume(:),volume_old(:),dvolume(:),volume_average(:)
real(8),allocatable,save::vol0(:),vol(:),vol_old(:),dvol(:),vol_average(:)  !单元密度，现时面积与上一时刻的面积
real(8),allocatable,save::e(:),de(:),p(:),dp(:)         !各单元的内能,平均应力,总内能

real(8),allocatable,save::strs_yld_11(:),strs_yld_22(:),strs_yld_33(:),strs_yld_12(:),strs_yld_13(:),strs_yld_23(:)
real(8),allocatable,save::strs_yld_123(:,:,:)

real(8),allocatable,save::rate(:),strain_rate_eff(:)

real(8),allocatable,save::energy_der(:),energy_xray_total(:),energy_xray(:)



real(8),allocatable,save::c11(:),c22(:),c33(:),c12(:),c21(:),c13(:),c31(:),c23(:),c32(:),c44(:),c55(:),c66(:)         !各单元本构系数
real(8),save::c123(6,6,360000)   !这个有问题
!real(8),allocatable,save::c123(:,:,:) 

real(8),allocatable,save::ex(:),ey(:),ez(:),gxy(:),gxz(:),gyz(:),pxy(:),pxz(:),pyz(:)
real(8),allocatable,save::strs_yld(:)!各向同性材料的屈服强度,前者为以层表示，后者以单元表示
real(8),allocatable,save::f1_yield(:)
real(8),allocatable,save::strs_frct(:)
real(8),allocatable,save::strs_frct_123(:,:,:)
real(8),allocatable,save::strs_frct_11(:),strs_frct_22(:),strs_frct_33(:),strs_frct_12(:),strs_frct_13(:),strs_frct_23(:)
integer(4),allocatable,save::plastic_index(:)   !塑性变形指示因子，为0时未发生塑性变形，为1时发生塑性变形
integer(4),allocatable,save::gas(:)  !0代表物质不是气体，1代表物质为气体

real(8),allocatable,save::stress_xx(:),stress_yy(:),stress_zz(:),stress_xy(:),stress_xz(:),stress_yz(:) !每个单元（l=0,s=0,m=0）处的应力
real(8),allocatable,save::stress_xyz(:,:,:)
real(8),allocatable,save::stress_11(:),stress_22(:),stress_33(:),stress_12(:),stress_13(:),stress_23(:) !在主轴坐标中，每个单元（l=0,s=0,m=0）处的应力
real(8),allocatable,save::stress_123(:,:,:)
real(8),allocatable,save::stress_d_11(:),stress_d_22(:),stress_d_33(:),stress_d_12(:),stress_d_13(:),stress_d_23(:)   !每个单元（l=0,s=0,m=0）处当前时刻的偏应力
real(8),allocatable,save::stress_d_123(:,:,:)



!球—偏分解专用变量

real(8),allocatable,save::q(:),q_old(:),q_av(:)   !人工粘性项（当前时刻和上一时刻的值及平均值）
!==========================================================================================三维的沙漏修正
real(8),allocatable,save::hggama1(:),hggama2(:),hggama3(:),hggama4(:)  !每个单元沙漏基矢量γ的四个分量（三维情况下四个是否够？？）
real(8),allocatable,save::f_hg(:,:)     !每个节点处的x,y方向的沙漏力

real(8),allocatable,save::strain_xx(:),strain_yy(:),strain_zz(:),strain_xy(:),strain_xz(:),strain_yz(:)      !每个单元（l=0,s=0,m=0）处的应变,等效塑性应变率
real(8),allocatable,save::strain_xyz(:,:,:)

real(8),allocatable,save::strain_xx_rate(:),strain_yy_rate(:),strain_zz_rate(:),strain_xy_rate(:),strain_xz_rate(:),strain_yz_rate(:)   !应变率分量
real(8),allocatable,save::strain_xyz_rate(:,:,:)
real(8),allocatable,save::strain_11(:),strain_22(:),strain_33(:),strain_12(:),strain_13(:),strain_23(:)      !在主轴坐标中，每个单元（l=0,s=0,m=0）处的应变
real(8),allocatable,save::strain_123(:,:,:)

real(8),allocatable,save::strain_11_rate(:),strain_22_rate(:),strain_33_rate(:),strain_12_rate(:),strain_13_rate(:),strain_23_rate(:)   !主轴坐标系中应变率分量
real(8),allocatable,save::strain_p_11(:),strain_p_22(:),strain_p_33(:),strain_p_12(:),strain_p_13(:),strain_p_23(:)      !在主轴坐标中，每个单元（l=0,s=0,m=0）处的塑性应变
real(8),allocatable,save::strain_p_123(:,:,:)

real(8),allocatable,save::strain_p_eff(:),strain_p_rate_eff(:)  !等效塑性应变,等效塑性应变增量，等效塑性应变率

real(8),allocatable,save::xyz0(:,:),xyz(:,:)  !j=1表示x坐标,j=2表示y坐标,j=3表示z坐标
real(8),allocatable,save::u(:,:),du(:,:),v(:,:),a(:,:),v_mid(:,:)  !在::个自由度上的节点位移,节点速度,节点加速度,中间时刻节点速度
real(8),allocatable,save::mass(:)     !在:个自由度上的质量矩阵（对角化）
real(8),allocatable,save::x0(:),y0(:),z0(:)   !初始时每单元中心的整体坐标
real(8),allocatable,save::x(:),y(:),z(:),v_cen_x(:),v_cen_y(:),v_cen_z(:)
real(8),allocatable,save::f_ext(:,:),f_int(:,:),f_ext_old(:,:),f_int_old(:,:) 
real(8),allocatable,save::stress_node_xx(:),stress_node_yy(:),stress_node_zz(:),stress_node_xy(:),stress_node_xz(:),stress_node_yz(:),p_node(:),e_node(:) 
real(8),save::l,h,m      !每单元中心的整体坐标，速度,自然坐标





integer(4),save::rad_energy_num !时间份数=时间谱宽度/dt
real(8),save::rad_energy(500000)   !功率密度(energy_num维)
real(8),save::wl(5000)
real(8),save::wer(5000)





real(8),save::parameter_A1,parameter_A2,parameter_A3,parameter_B0,parameter_T1                  
real(8),save::a_11,a_22,a_33,a_12,a_13,a_23,a_44,a_55,a_66,sigma_well_1,sigma_well_10,mean_plastic    
real(8),save::tensil_11,tensil_22,tensil_33,tensil_12,tensil_23,tensil_31                        
real(8),save::fracture_e11,fracture_e22,fracture_e33,fracture_e12,fracture_e23,fracture_e31     
real(8),save::damage_coupling





real(8),save::w_int,w_ext,w_kin,w_orgn,w_err     !内能，外能，动能，初始能量，能量误差百分比
!体积计算相关变量--------------------------------
real,save::core(3,1)
real,save::volume_s(6)
!------------------------------------------------
!对称面--------------------------------------------
integer(4),save::symm_xy(100000),symm_xz(100000)
!----------------------------------------------------
logical,save::mat0,mat1   !mat0为真时表示各向同性材料，mat1为真时表示各向异性材料
real(8),save::h_dt   !半时间步长

end module vars   