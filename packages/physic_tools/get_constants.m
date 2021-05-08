function [ constants ] = get_constants(  )
% 返回 存储常数 的结构体
%--------物理常数--------------------------------
constants.e=1.6022e-19;%基本电荷 [C]
constants.me=9.1094e-31;%电子质量 [kg]
constants.mH=1.6726e-27;%质子质量 [kg]
constants.mH2=3.3474e-27;%质子质量 [kg]
constants.c=2.9979e8;%真空光速[m/s]
constants.eps0=8.8542e-12;%真空介电常数 [F/m]
constants.mu0=4*pi*1e-7;%真空磁导率 [H/m]
constants.kB=1.3807e-23;%玻尔兹曼常数 [ J/K]
constants.amu=1.66053904e-27; %原子质量单位[kg]

constants.q_m_ratio_e=-constants.e/constants.me; %电子荷质比 [C/kg]
constants.q_m_ratio_H=constants.e/constants.mH; 
constants.q_m_ratio_H2p=constants.e/constants.mH/2; 

constants.sigma_Cu=5.8e7;       %铜的电导率

%--------物理常数--------------------------------
end
