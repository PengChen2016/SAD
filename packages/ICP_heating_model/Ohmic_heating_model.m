function [ plasma ] = Ohmic_heating_model( plasma )
% Ohmic_heating_model
constants=get_constants();
Xsec=get_Xsec(false);
%% 电子与中性粒子动量转移碰撞
plasma.kenp=get_k(plasma.Te, Xsec.enp); %速率系数
plasma.nu_enp=plasma.ng.*plasma.kenp; % 碰撞频率
lambda_coll_enp=plasma.ve./plasma.nu_enp; %平均自由程
%% 电子与中性粒子电离碰撞
plasma.keniz=get_k(plasma.Te, Xsec.eniz);
plasma.nu_eniz=plasma.ng.*plasma.keniz;
lambda_coll_eniz=plasma.ve./plasma.nu_eniz;
%% 电子与离子动量损失碰撞
lambda_e=(constants.eps0*constants.e*plasma.Te.^(3/2))./sqrt(plasma.ne)*...
    (4*pi*constants.mH/(constants.e^3)/(constants.me+constants.mH));%库仑算子中Λ
plasma.nu_eip=plasma.ne.*log(lambda_e)./sqrt((constants.e*plasma.Te).^3*constants.me)*...
    4*sqrt(2*pi)/3*((constants.e^2/4/pi/constants.eps0)^2);
lambda_coll_eip=plasma.ve./plasma.nu_eip;

%% total
plasma.nu_m=plasma.nu_enp+plasma.nu_eniz+plasma.nu_eip;

end