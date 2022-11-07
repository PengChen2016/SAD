function [ plasma ] = get_characterstic_parameters( flag, plasma )
% abstract 等离子体特征参数计算
constants=get_constants();
plasma.wpe=get_omega_pe(plasma.ne); %电子等离子体频率
switch flag.type_Xsec(3:4)
    case 'H2'
        plasma.wpi=get_omega_pi(plasma.ne,1,1); %离子等离子体频率
    case 'Ar'
        plasma.wpi=get_omega_pi(plasma.ne,1,39.948); %离子等离子体频率
end

plasma.ve=sqrt(8*plasma.Te*constants.e/(pi*constants.me));    %电子平均速率，计算自由程用
plasma.veth=sqrt(2*plasma.Te*constants.e/constants.me);    %电子热速率
%             wce=constants.e;     %电子拉莫尔运动频率
end