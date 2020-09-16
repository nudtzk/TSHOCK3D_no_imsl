module vars                            !Ӧ�ó���ɳ©������Ҫ�����������������
implicit none
!namelist:compute_control-----------------------------------------------------------------------------------------------------------
character(len=50),save:: problem
integer(4), save :: problem_mode
integer(4),save:: shape                 !0����ƽ����״��1����Բ��������״
integer(4),save:: eos_mode  
real(8),save:: flyer_length
integer(4),save::lmax,lmax_iso,lmax_aniso     !����ͬ�Բ�����������Բ��ϵĲ���  
real(8),save:: length_1,length_2,length_3
real(8),save::dl_1,dl_2,dl_3  
real(8),save::q1,q2,qhg    !ճ��ϵ��!ɳ©ճ��ϵ��
real(8),save:: t,dt,t_total   !ʱ�����
integer(4),save::nscr
integer(8),save::ne,nj
real(8),save::vf
integer(4),save:: num_elefix_x,num_elefix_y,num_elefix_z  !ѡȡx/y/z������ͬ��һ�������ȡ20��
integer(4),save:: ne_his_x(20),ne_his_y(20),ne_his_z(20)  !��Ӧ����
integer(4),save::num_time_record  !ָ����������ʱ������,���=20
real(8),save::time_record(20)     !ָ����Ҫ�������ľ���ʱ��ֵ
!----------------------------------------------------------------------------------------------------------------------------------------

!namelist layer_setting--------------------------------------------------------
integer(4),save::mat_model(30)     !�����������ָʾ��0->����ͬ�ԣ�1->��������
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
integer(4),save::layer_no_aniso(30) !�������Բ����е�i���������������е�λ�ã���Ҫ������������
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
integer(4),save::n_blackbody    !�����׵ĸ��������ȡ20��
real(8),save::t_kev(20)         !��������ֵĺ����¶ȣ���λΪkeV
real(8),save::b_energy(20)      !�����������ռ�������������ݶ�
real(8),save::wl_short,wl_long
integer(4),save::rad_time_mode  !ʱ����ģʽ��0Ϊ�����ף�1Ϊ��������
real(8),save::rad_time_width    !ʱ���׵Ŀ�ȣ���������Ϊ��߿�
real(8),save::rad_energy_total  !����ͨ���ܶ�
integer(4),parameter::ken=5000
integer(4),parameter:: max_element_num=10
integer(4),parameter::max_layer=30
type::layer_xray
    integer(4) element_num    !�����Ԫ��������<=10
	integer(4) element_z(max_element_num)    !Ԫ�ص�ԭ��������Ԫ�����಻����10��
	real(8) element_wp(max_element_num)     !��Ԫ�ص������ٷֱ�
	real(8) mu(ken)
end type layer_xray
type(layer_xray), save :: xlayer(max_layer)

!------------------------------------------------------

!namelist xray_output-----------------------------------------------------------------------------------------------------------
integer(4),save::cm_r1(5000),cm_r2(5000)  
integer(4),save::cm_theta_0(2000),cm_theta_10(2000),cm_theta_20(2000),cm_theta_30(2000),cm_theta_40(2000),cm_theta_45(2000),   &
                 cm_theta_50(2000),cm_theta_60(2000),cm_theta_70(2000),cm_theta_80(2000),cm_theta_90(2000)
!����������ֲ���-------------------------------------------------------------------------------------
character(len=18) time_output(20)                                    !������ʽΪtime_01_output.txt
character(len=13) fne_his_x(20),fne_his_y(20),fne_his_z(20)         !������ʽΪfix_ex_01.dat,fix_ey_01.dat
integer(4),save::st,sn  !���ڼ��㲽���Ŀ��Ʊ���
!����û�ò���������dx����ƽ����----------------------------------------------------------------------------------
real(8),save::length_x,length_y,length_z,lx(30),dx,dy,dz          !��Ӧƽ��x�����ܳ��ȣ�y�����ܳ��ȣ�z�����ܳ��ȸ���ĺ�ȣ�x����y����z���������ߴ�
real(8),save:: r_in,r_out,r_height,lr(30),dr,d_theta,d_height     !�ֱ�ΪԲ�����ڰ뾶����뾶�������ذ뾶����ĺ�ȣ��뾶���������ߴ�,�ǶȻ��ֳߴ�
!-------------------------------------------------------------------------------------
integer(4),save::n_cmup  !�������Ϸ��ĵ�Ԫ����
integer(4),save::cmup(10000)    !��x����������Ϸ��������ĵ�Ԫ(ƽ��ģ��)

real(8),save::ex_0(30),ey_0(30),ez_0(30),gxy_0(30),gxz_0(30),gyz_0(30)!,pxy(30),pxz(30),pyz(30)!��������ģ��������ģ�������ɱ�
real(8),save::c11_0(30),c22_0(30),c33_0(30),c12_0(30),c13_0(30),c23_0(30),c44_0(30),c55_0(30),c66_0(30)         !��ʼ���㱾��ϵ�����ɸնȾ���ת���õ�,����3*3��������Ӧ����ϵ��3�����й�ϵ
real(8),save::c123_0(6,6,20)
real(8),save::rho0(30)  !����ĳ�ʼ�ܶ�

real(8),save::shape_n(8),shape_n0(8) !��״����,���ĵ����״����
real(8),save::dn_l(8),dn_h(8),dn_m(8),dn0_l(8),dn0_h(8),dn0_m(8)  !��״���������ĵ����״��������Ȼ����l,h,m�ĵ���
real(8),save::dn_lhm(8,3)
real(8),save::dn_xyz(8,3)
real(8),save::coord_xyz(3,8)
real(8),save::coord_xyz0(3,8)
real(8),save::dn_x(8),dn_y(8),dn_z(8)           !��״��������������x,y,z�ĵ���
real(8),save::jj(3,3),jj_inverse(3,3)       !(x,y,z)��(l,h,m)ƫ�������Ÿ��Ⱦ���,��������ʽ�������
real(8),save::jj0(3,3),jj0_inverse(3,3)    !(x0,y0,z0)��(l,h,m)ƫ�������Ÿ��Ⱦ���,��������ʽ�������
real(8),save::f11,f12,f13,f21,f22,f23,f31,f32,f33                  !(x,y,z)��(x0,y0,z0)��ƫ�������Ÿ��Ⱦ��󣬱����ݶȾ���
real(8),save::deform(3,3)
real(8),save::rotate_R(3,3)
real(8),save::rotate_original(3,3)

!-----------------------------------------------
real(8),save::c0(30),s(30),gama(30)  !��̬���̲���
real(8),save::h0(30),es(30)      !PUFF��̬���̲�����������
real(8),save::ax(30),ay(30),az(30),nx(30),ny(30),nz(30),axy(30),nxy(30),axz(30),nxz(30),ayz(30),nyz(30)
real(8),save::cn(30)
real(8),save::beta0(30),strain_rate0(30)        !Ӧ������Ч���ӵĲ���0.02183,��ʼ�ο�Ӧ���� ȡֵΪ0.001/s,Ӧ��������
!---------------------------------
real(8),save::diag_vector(3)
real(8),save::diag_matrix(3,3)
real(8),save::diag_sum
real(8),save::p_dst,p_dsp,ee         !ƫӦ�乱��ѹ��������Ӧ�乱��ѹ��,������




integer(4),allocatable,save::fail(:)  !���ϵ�Ԫ��ǰʧЧ�ж���0�������δʧЧ��1�������ʧЧ
integer(4),allocatable,save::ne_layer(:)  !����Ԫ��������
integer(4),allocatable,save::ia(:,:)   !j=1-8��ʾ���ĳ��Ԫ��8���ڵ���,j=9����õ�Ԫ���ڲ���������
integer(4),allocatable,save::index_ele_write(:)         !����Բ������ģ�ͣ�index_ele_write(i)==1����õ�Ԫ�����ϰ�ƽ�棬��Ҫ������

real(8),allocatable,save::rho(:),rho_old(:),volume0(:),volume(:),volume_old(:),dvolume(:),volume_average(:)
real(8),allocatable,save::vol0(:),vol(:),vol_old(:),dvol(:),vol_average(:)  !��Ԫ�ܶȣ���ʱ�������һʱ�̵����
real(8),allocatable,save::e(:),de(:),p(:),dp(:)         !����Ԫ������,ƽ��Ӧ��,������

real(8),allocatable,save::strs_yld_11(:),strs_yld_22(:),strs_yld_33(:),strs_yld_12(:),strs_yld_13(:),strs_yld_23(:)
real(8),allocatable,save::strs_yld_123(:,:,:)

real(8),allocatable,save::rate(:),strain_rate_eff(:)

real(8),allocatable,save::energy_der(:),energy_xray_total(:),energy_xray(:)



real(8),allocatable,save::c11(:),c22(:),c33(:),c12(:),c21(:),c13(:),c31(:),c23(:),c32(:),c44(:),c55(:),c66(:)         !����Ԫ����ϵ��
real(8),save::c123(6,6,360000)   !���������
!real(8),allocatable,save::c123(:,:,:) 

real(8),allocatable,save::ex(:),ey(:),ez(:),gxy(:),gxz(:),gyz(:),pxy(:),pxz(:),pyz(:)
real(8),allocatable,save::strs_yld(:)!����ͬ�Բ��ϵ�����ǿ��,ǰ��Ϊ�Բ��ʾ�������Ե�Ԫ��ʾ
real(8),allocatable,save::f1_yield(:)
real(8),allocatable,save::strs_frct(:)
real(8),allocatable,save::strs_frct_123(:,:,:)
real(8),allocatable,save::strs_frct_11(:),strs_frct_22(:),strs_frct_33(:),strs_frct_12(:),strs_frct_13(:),strs_frct_23(:)
integer(4),allocatable,save::plastic_index(:)   !���Ա���ָʾ���ӣ�Ϊ0ʱδ�������Ա��Σ�Ϊ1ʱ�������Ա���
integer(4),allocatable,save::gas(:)  !0�������ʲ������壬1��������Ϊ����

real(8),allocatable,save::stress_xx(:),stress_yy(:),stress_zz(:),stress_xy(:),stress_xz(:),stress_yz(:) !ÿ����Ԫ��l=0,s=0,m=0������Ӧ��
real(8),allocatable,save::stress_xyz(:,:,:)
real(8),allocatable,save::stress_11(:),stress_22(:),stress_33(:),stress_12(:),stress_13(:),stress_23(:) !�����������У�ÿ����Ԫ��l=0,s=0,m=0������Ӧ��
real(8),allocatable,save::stress_123(:,:,:)
real(8),allocatable,save::stress_d_11(:),stress_d_22(:),stress_d_33(:),stress_d_12(:),stress_d_13(:),stress_d_23(:)   !ÿ����Ԫ��l=0,s=0,m=0������ǰʱ�̵�ƫӦ��
real(8),allocatable,save::stress_d_123(:,:,:)



!��ƫ�ֽ�ר�ñ���

real(8),allocatable,save::q(:),q_old(:),q_av(:)   !�˹�ճ�����ǰʱ�̺���һʱ�̵�ֵ��ƽ��ֵ��
!==========================================================================================��ά��ɳ©����
real(8),allocatable,save::hggama1(:),hggama2(:),hggama3(:),hggama4(:)  !ÿ����Ԫɳ©��ʸ���õ��ĸ���������ά������ĸ��Ƿ񹻣�����
real(8),allocatable,save::f_hg(:,:)     !ÿ���ڵ㴦��x,y�����ɳ©��

real(8),allocatable,save::strain_xx(:),strain_yy(:),strain_zz(:),strain_xy(:),strain_xz(:),strain_yz(:)      !ÿ����Ԫ��l=0,s=0,m=0������Ӧ��,��Ч����Ӧ����
real(8),allocatable,save::strain_xyz(:,:,:)

real(8),allocatable,save::strain_xx_rate(:),strain_yy_rate(:),strain_zz_rate(:),strain_xy_rate(:),strain_xz_rate(:),strain_yz_rate(:)   !Ӧ���ʷ���
real(8),allocatable,save::strain_xyz_rate(:,:,:)
real(8),allocatable,save::strain_11(:),strain_22(:),strain_33(:),strain_12(:),strain_13(:),strain_23(:)      !�����������У�ÿ����Ԫ��l=0,s=0,m=0������Ӧ��
real(8),allocatable,save::strain_123(:,:,:)

real(8),allocatable,save::strain_11_rate(:),strain_22_rate(:),strain_33_rate(:),strain_12_rate(:),strain_13_rate(:),strain_23_rate(:)   !��������ϵ��Ӧ���ʷ���
real(8),allocatable,save::strain_p_11(:),strain_p_22(:),strain_p_33(:),strain_p_12(:),strain_p_13(:),strain_p_23(:)      !�����������У�ÿ����Ԫ��l=0,s=0,m=0����������Ӧ��
real(8),allocatable,save::strain_p_123(:,:,:)

real(8),allocatable,save::strain_p_eff(:),strain_p_rate_eff(:)  !��Ч����Ӧ��,��Ч����Ӧ����������Ч����Ӧ����

real(8),allocatable,save::xyz0(:,:),xyz(:,:)  !j=1��ʾx����,j=2��ʾy����,j=3��ʾz����
real(8),allocatable,save::u(:,:),du(:,:),v(:,:),a(:,:),v_mid(:,:)  !��::�����ɶ��ϵĽڵ�λ��,�ڵ��ٶ�,�ڵ���ٶ�,�м�ʱ�̽ڵ��ٶ�
real(8),allocatable,save::mass(:)     !��:�����ɶ��ϵ��������󣨶Խǻ���
real(8),allocatable,save::x0(:),y0(:),z0(:)   !��ʼʱÿ��Ԫ���ĵ���������
real(8),allocatable,save::x(:),y(:),z(:),v_cen_x(:),v_cen_y(:),v_cen_z(:)
real(8),allocatable,save::f_ext(:,:),f_int(:,:),f_ext_old(:,:),f_int_old(:,:) 
real(8),allocatable,save::stress_node_xx(:),stress_node_yy(:),stress_node_zz(:),stress_node_xy(:),stress_node_xz(:),stress_node_yz(:),p_node(:),e_node(:) 
real(8),save::l,h,m      !ÿ��Ԫ���ĵ��������꣬�ٶ�,��Ȼ����





integer(4),save::rad_energy_num !ʱ�����=ʱ���׿��/dt
real(8),save::rad_energy(500000)   !�����ܶ�(energy_numά)
real(8),save::wl(5000)
real(8),save::wer(5000)





real(8),save::parameter_A1,parameter_A2,parameter_A3,parameter_B0,parameter_T1                  
real(8),save::a_11,a_22,a_33,a_12,a_13,a_23,a_44,a_55,a_66,sigma_well_1,sigma_well_10,mean_plastic    
real(8),save::tensil_11,tensil_22,tensil_33,tensil_12,tensil_23,tensil_31                        
real(8),save::fracture_e11,fracture_e22,fracture_e33,fracture_e12,fracture_e23,fracture_e31     
real(8),save::damage_coupling





real(8),save::w_int,w_ext,w_kin,w_orgn,w_err     !���ܣ����ܣ����ܣ���ʼ�������������ٷֱ�
!���������ر���--------------------------------
real,save::core(3,1)
real,save::volume_s(6)
!------------------------------------------------
!�Գ���--------------------------------------------
integer(4),save::symm_xy(100000),symm_xz(100000)
!----------------------------------------------------
logical,save::mat0,mat1   !mat0Ϊ��ʱ��ʾ����ͬ�Բ��ϣ�mat1Ϊ��ʱ��ʾ�������Բ���
real(8),save::h_dt   !��ʱ�䲽��

end module vars   