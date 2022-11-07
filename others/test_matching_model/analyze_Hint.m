% 大源 匹配电路模型
clear
addpath(genpath('../../packages'))

%% 粗糙尝试-使用2016岳海昆数据
f=1.02e6;
fprintf('---------- %d MHz\n',f)
%--------------------- 初始化
% 2016岳海昆 9匝线圈
ZD0.R=1.8;
ZD0.L=20.5e-6;
Zs0=50;
Cs0=1.3e-9;  % 运行时参数，非设计参数
Cp0=7.6e-9;
% 网分 直接测量 RF变 应为1MHz
t_YHK.L1 = 135.4e-6;             %变压器原边电感
t_YHK.L2 = 16.8e-6;              %变压器副边线圈电感
t_YHK.k =0.96 ;                     %互感耦合系数
t_YHK.M = t_YHK.k*sqrt(t_YHK.L1*t_YHK.L2);              %互感

c_YHK=c_RF_circuit(f);
% init matching t
c_YHK=c_YHK.init_matching_t('m',t_YHK.M,t_YHK.L1,t_YHK.L2);
% init matching C
c_YHK=c_YHK.init_matching_C( 'gamma', Cs0, Cp0);
% init Zs
c_YHK.s.Z.cm=Zs0;

%--------------------- get ZD
fprintf('ZD(网分) = %.2f Ω+ %.2f μH\n',ZD0.R,ZD0.L*1e6)
c_YHK=c_YHK.get_ZD();
fprintf('ZD(t_YHK) = %.2f Ω+ %.2f μH\n',c_YHK.D.Z.R,c_YHK.D.Z.L*1e6)
% fprintf('ZD(网分) = %.2f + %.2f Ω\n',ZD0.R,ZD0.X)

%--------------------- get Zs and compare with Multisim
c_YHK=c_YHK.init_ZD('RL',1,ZD0);
c_YHK=c_YHK.get_Zs();
fprintf('Zs = %.2f +∠ %.2f \n',c_YHK.s.Z.abs,c_YHK.s.Z.phase)
fprintf('ZD = %.2f +∠ %.2f \n',c_YHK.D.Z.abs,c_YHK.D.Z.phase)